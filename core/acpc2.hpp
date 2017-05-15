#pragma once

#include "cc_sampler.hpp"
#include "cluster_vertex.hpp"
#include "counts_cache.hpp"
#include "experiment_reporter.hpp"
#include "guesser.hpp"
#include "require.hpp"
#include "termcolor.hpp"
#include "types.hpp"

void pick_centers(const ugraph_t &graph,
                  const std::vector<ClusterVertex> &vinfo,
                  const size_t num_centers,
                  Xorshift1024star &rnd,
                  ConnectionCountsCache &cccache,
                  std::vector<ugraph_vertex_t> &uncovered_scratch,
                  std::vector<ugraph_vertex_t> &selected_centers) {
  selected_centers.clear();
  uncovered_scratch.clear();
  auto n = boost::num_vertices(graph);
  REQUIRE(num_centers < n, "There are not enough nodes in the graph.");
  for (ugraph_vertex_t i = 0; i < n; i++) {
    if (!vinfo[i].is_covered()) {
      uncovered_scratch.push_back(i);
    } else if (vinfo[i].is_covered() && cccache.contains(i)) {
      // mark the node for eviction from the cache
      cccache.set_accessed(i, 0);
    }
  }
  if (uncovered_scratch.size() == 0) {
    throw std::logic_error("No uncovered node to select");
  }

  std::random_shuffle(
      uncovered_scratch.begin(), uncovered_scratch.end(),
      [&](size_t i) { return (size_t)std::floor(rnd.next_double() * i); });

  size_t to_select = std::min(num_centers, uncovered_scratch.size());

  selected_centers.resize(to_select);
  std::copy(uncovered_scratch.begin(), uncovered_scratch.begin() + to_select,
            selected_centers.begin());
}

double
sum_center_connection_probabilities(const std::vector<ClusterVertex> &vinfo) {
  double sum = 0.0;
  for (const ClusterVertex &v : vinfo) {
    if (v.is_covered()) {
      sum += v.probability();
    }
  }
  return sum;
}

template <typename Sampler>
ugraph_vertex_t best_center(
    // Input parameters
    const ugraph_t &graph,
    const size_t batch_centers,
    const probability_t p_bound,
    // Input
    const std::vector<ClusterVertex> &current_clustering,
    // Output
    std::vector<probability_t> &best_probabilities,
    // Utilities
    Sampler &sampler, Xorshift1024star &rnd, ConnectionCountsCache &cccache,
    // Scratch areas
    std::vector<ugraph_vertex_t> &centers,
    std::vector<ugraph_vertex_t> &uncovered_scratch,
    std::vector<probability_t> &scratch_probabilities) {

  const size_t n = boost::num_vertices(graph);
  
  pick_centers(graph, current_clustering, batch_centers, rnd, cccache,
               uncovered_scratch, centers);

  size_t max_newly_covered = 0;
  ugraph_vertex_t best_center;
  for (ugraph_vertex_t candidate_center : centers) {
    sampler.connection_probabilities_cache(graph, candidate_center, cccache,
                                           scratch_probabilities);
    size_t newly_covered = 0;
    for (ugraph_vertex_t i=0; i<n; i++) {
      if (!current_clustering[i].is_covered() && scratch_probabilities[i] >= p_bound) {
        newly_covered++;
      }
    }
    if (newly_covered > max_newly_covered) {
      max_newly_covered = newly_covered;
      best_center = candidate_center;
      std::copy(scratch_probabilities.begin(), scratch_probabilities.end(),
                best_probabilities.begin());
    }
  }
  return best_center;
}

template<typename Sampler>
std::vector< ClusterVertex >
average_connection_probability_clustering_2(const ugraph_t & graph,
                                            Sampler & sampler,
                                            Xorshift1024star & rnd,
                                            const size_t k,
                                            const size_t batch_centers,
                                            const double rate,
                                            const probability_t p_low) {

  // Workspace areas
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > scratch_probabilities(n);
  std::vector< probability_t > best_probabilities(n);
  std::vector< ugraph_vertex_t > uncovered_scratch(n);
  std::vector< ugraph_vertex_t > centers(batch_centers);
  ConnectionCountsCache cccache(std::min(n, 2000ul));
  size_t iteration = 0;
  probability_t p_curr = 1.0;
  probability_t best_sum_p = 0.0;
  ExponentialGuesser guesser(rate, p_low);
  size_t uncovered = n;

  while (!guesser.stop()) {
    LOG_INFO(">>> Build clustering with p_curr=" << p_curr);
    cccache.cleanup();
    LOG_DEBUG(cccache.str());
    std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());
    uncovered = n;
    sampler.min_probability(graph, p_curr);

    // Build the clustering. Start the count from 1 because the
    // stopping condition is _inside_ the cycle
    for (size_t center_cnt = 1; center_cnt < k; center_cnt++) {
      ugraph_vertex_t center = best_center(graph, batch_centers, p_curr, vinfo,
                                           best_probabilities,
                                           sampler, rnd, cccache, centers,
                                           uncovered_scratch, scratch_probabilities);
      vinfo[center].make_center(center);
      uncovered--;
      
      // Cover the nodes
      for (ugraph_vertex_t i = 0; i < n; i++) {
        if (best_probabilities[i] >= p_curr) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, best_probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < best_probabilities[i]) {
            vinfo[i].cover(center, best_probabilities[i]);
          } else if (!vinfo[i].is_covered() &&
                     vinfo[i].unreliable_probability() < best_probabilities[i]) {
            vinfo[i].unreliable_cover(center, best_probabilities[i]);
          }
        }
      }
      if (center_cnt + uncovered <= k) {
        // Complete the clustering
        for (ugraph_vertex_t i = 0; i < n; i++) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].make_center(i);
            uncovered--;
          }
        }
        break;
      }
    }

    probability_t sum_p = sum_center_connection_probabilities(vinfo);
    if (sum_p > best_sum_p) {
      LOG_INFO("Better average probability: " <<
               termcolor::green << sum_p / n << termcolor::reset);
      best_sum_p = sum_p;
      for (ugraph_vertex_t i = 0; i < n; i++) {
        valid_clustering[i] = vinfo[i];
      }
    } else {
      LOG_INFO("Average probability: " << sum_p / n);
    }
    
    // Check if the corresponding k-center clustering is valid
    if (uncovered == 0) {
      guesser.below();
    } else {
      guesser.above();
    }

    EXPERIMENT_APPEND("average-probability",
                  {{"p_curr", p_curr},
                   {"average-probability", sum_p / n},
                   {"iteration", iteration}});
    
    LOG_INFO("Still " << uncovered << " nodes to cover");
    LOG_INFO("Cache hit rate: " << std::fixed << std::setprecision(2)
             << cccache.perc_hits() << "%");
    // update the probability
    p_curr = guesser.guess();
    iteration++;
  }

    size_t returned_uncovered = 0;
  for(ugraph_vertex_t i = 0; i<n; i++) {
    if (!valid_clustering[i].is_covered()) {
      returned_uncovered++;
    }
  }
  LOG_INFO("There are " << returned_uncovered << " uncovered nodes in the solution");
  if (returned_uncovered > 0) {
    // Fix the uncovered nodes by adding probability guesses. We know
    // that p_curr is sufficient to cover the whole graph. Probably the
    // centers are not in the cache.
    sampler.min_probability(graph, p_curr);
    if (returned_uncovered >= k) {
      // Sample from the centers
      LOG_INFO("Sampling to cover the uncovered nodes from centers");
      for(ugraph_vertex_t i = 0; i<n; i++) {
        if (valid_clustering[i].is_center()) {
          sampler.connection_probabilities_cache(graph, i, cccache, best_probabilities);
          for(ugraph_vertex_t j = 0; j<n; j++) {
            if (!valid_clustering[j].is_covered() ||
                best_probabilities[j] > valid_clustering[j].probability()) {
              valid_clustering[j].cover(i, best_probabilities[j]);
            }
          }
        }
      }
    } else {
      // Sample from the uncovered nodes
      LOG_INFO("Sampling to cover uncovered nodes from the nodes themselves");
      for(ugraph_vertex_t i = 0; i<n; i++) {
        if (!valid_clustering[i].is_covered()) {
          sampler.connection_probabilities_cache(graph, i, cccache, best_probabilities);
          probability_t best_p = 0.0;
          ugraph_vertex_t best_c = 0;
          for(ugraph_vertex_t j = 0; j<n; j++) {
            if (valid_clustering[j].is_center() && best_probabilities[j] > best_p) {
              best_p = best_probabilities[j];
              best_c = j;
            }
          }
          valid_clustering[i].cover(best_c, best_p);
        }
      }
    }
  }

  return valid_clustering;

}
