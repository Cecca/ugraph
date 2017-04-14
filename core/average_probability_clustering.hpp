#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include "guesser.hpp"
#include "counts_cache.hpp"

size_t count_uncovered(const std::vector< ClusterVertex > & vinfo) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered()) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const std::vector< ClusterVertex > & vinfo) {
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered()) {
      return i;
    }
  }
  throw std::logic_error("No uncovered node to select");
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            ConnectionCountsCache & cccache,
                            const std::vector< ClusterVertex > & vinfo) {
  int query_result = cccache.uncovered_node(vinfo);
  if (query_result >= 0) {
    return (ugraph_vertex_t) query_result;
  } else {
    return pick_vertex(graph, vinfo);
  }
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
                            const size_t k,
                            const double rate,
                            const probability_t p_low,
                            ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > probabilities(n);
  ConnectionCountsCache cccache(std::min(k, 500ul));
  size_t iteration = 0;
  probability_t p_curr = 1.0;
  probability_t reliable_estimate_lower_bound = 1.0;
  // TODO Implement exponential guesser with better binary search,
  // that searched from the second-to-last increasing p_min
  ExtendedExponentialGuesser guesser(rate, p_low);
  size_t uncovered = n;
  double max_sum = 0.0;

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
      ugraph_vertex_t center = pick_vertex(graph, cccache, vinfo);
      vinfo[center].make_center(center);
      uncovered--;
      sampler.connection_probabilities_cache(graph, center, cccache, probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i=0; i<n; i++) {
        if (probabilities[i] >= reliable_estimate_lower_bound) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
          }
        }
      }

      if (uncovered == 0) {
        break;
      }

    }
    double prob_sum = sum_center_connection_probabilities(vinfo);
    LOG_INFO("Average connection probability " << prob_sum / n);

    if (prob_sum >= max_sum) {
      guesser.above();
      for (ugraph_vertex_t i=0; i<n; i++) {
        valid_clustering[i] = vinfo[i];
      }
      max_sum = prob_sum;
    } else {
      guesser.below();
    }

    LOG_INFO("Cache hit rate: " << std::fixed << std::setprecision(2)
             << cccache.perc_hits() << "%");
    // update the probability
    p_curr = guesser.guess();
    reliable_estimate_lower_bound =
      (p_curr < reliable_estimate_lower_bound)? p_curr : reliable_estimate_lower_bound;
    iteration++;
  }

  // Fix uncovered nodes
  // get the first center, so to assign to it uncovered nodes
  ugraph_vertex_t first_center_idx=0;
  for (; first_center_idx<n || valid_clustering[first_center_idx].is_center(); first_center_idx++){}
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!valid_clustering[i].is_covered()) {
      valid_clustering[i].cover(first_center_idx, 0.0);
    }
  }

  return valid_clustering;
  
}
