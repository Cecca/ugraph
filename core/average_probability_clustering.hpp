#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include "guesser.hpp"
#include "counts_cache.hpp"
#include "termcolor.hpp"

size_t count_uncovered(const std::vector< ClusterVertex > & vinfo, double p_curr) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered() || v.probability() < p_curr) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            const std::vector< ClusterVertex > & vinfo) {
  auto n = boost::num_vertices(graph);
  // Give precedence to really uncovered nodes
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered()) {
      return i;
    }
  }
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (vinfo[i].probability() < p_curr) {
      return i;
    }
  }
  throw std::logic_error("No uncovered node to select");
}

ugraph_vertex_t pick_vertex_rnd(const ugraph_t & graph,  
                                Xorshift1024star & rnd,
                                ConnectionCountsCache & cccache,
                                probability_t p_curr,
                                std::vector<ugraph_vertex_t> & uncovered_scratch,
                                const std::vector< ClusterVertex > & vinfo) {
  uncovered_scratch.clear();
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered() || vinfo[i].probability() < p_curr) {
      uncovered_scratch.push_back(i);
    } else if (vinfo[i].is_covered() && cccache.contains(i)) {
      // mark the node for eviction from the cache
      cccache.set_accessed(i, 0);
    }
  }
  if (uncovered_scratch.size() == 0) {
    throw std::logic_error("No uncovered node to select");
  }
  size_t i = (size_t) std::floor(rnd.next_double()*uncovered_scratch.size());
  return uncovered_scratch[i];
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            ConnectionCountsCache & cccache,
                            const std::vector< ClusterVertex > & vinfo) {
  size_t n = vinfo.size();
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (cccache.contains(i)) {
      if (!vinfo[i].is_covered()) {
        return i;
      } else if (vinfo[i].is_covered() && !vinfo[i].is_center()) {
        // reset counter to mark for eviction for nodes that are
        // covered but are not centers.
        cccache.set_accessed(i, 0);
      }
    }
  }
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (cccache.contains(i)) {
      if (vinfo[i].probability() < p_curr) {
        return i;
      } else if (vinfo[i].is_covered() && !vinfo[i].is_center()) {
        // reset counter to mark for eviction for nodes that are
        // covered but are not centers.
        cccache.set_accessed(i, 0);
      }
    }
  }

  
  return pick_vertex(graph, p_curr, vinfo);
}

double sum_center_connection_probabilities(const std::vector< ClusterVertex > & vinfo) {
  double sum = 0.0;
  for (const ClusterVertex & v : vinfo) {
    if (v.is_covered()) {
      sum += v.probability();
    }
  }
  return sum;
}

template<typename Sampler>
std::vector< ClusterVertex >
average_probability_cluster(const ugraph_t & graph,
                            Sampler & sampler,
                            Xorshift1024star & rnd,
                            const size_t k,
                            const double rate,
                            const probability_t p_low,
                            ExperimentReporter & experiment) {

  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > probabilities(n);
  std::vector< ugraph_vertex_t > uncovered_scratch(n);
  ConnectionCountsCache cccache(std::min(k, 500ul));
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
      assert(uncovered == count_uncovered(vinfo));
      ugraph_vertex_t center = pick_vertex_rnd(graph, rnd, cccache, p_curr,
                                               uncovered_scratch, vinfo);
      vinfo[center].make_center(center);
      uncovered--;
      sampler.connection_probabilities_cache(graph, center, cccache,
                                             probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i = 0; i < n; i++) {
        if (probabilities[i] >= p_curr) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
          } else if (!vinfo[i].is_covered() &&
                     vinfo[i].unreliable_probability() < probabilities[i]) {
            vinfo[i].unreliable_cover(center, probabilities[i]);
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

    experiment.append("average-probability",
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
          sampler.connection_probabilities_cache(graph, i, cccache, probabilities);
          for(ugraph_vertex_t j = 0; j<n; j++) {
            if (!valid_clustering[j].is_covered() ||
                probabilities[j] > valid_clustering[j].probability()) {
              valid_clustering[j].cover(i, probabilities[j]);
            }
          }
        }
      }
    } else {
      // Sample from the uncovered nodes
      LOG_INFO("Sampling to cover uncovered nodes from the nodes themselves");
      for(ugraph_vertex_t i = 0; i<n; i++) {
        if (!valid_clustering[i].is_covered()) {
          sampler.connection_probabilities_cache(graph, i, cccache, probabilities);
          probability_t best_p = 0.0;
          ugraph_vertex_t best_c = 0;
          for(ugraph_vertex_t j = 0; j<n; j++) {
            if (valid_clustering[j].is_center() && probabilities[j] > best_p) {
              best_p = probabilities[j];
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
