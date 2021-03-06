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

size_t count_uncovered(const std::vector< ClusterVertex > & vinfo) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered()) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex_rnd(const ugraph_t & graph,  
                                Xorshift1024star & rnd,
                                ConnectionCountsCache & cccache,
                                std::vector<ugraph_vertex_t> & uncovered_scratch,
                                const std::vector< ClusterVertex > & vinfo) {
  uncovered_scratch.clear();
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered()) {
      uncovered_scratch.push_back(i);
    } else if (cccache.contains(i)) {
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

template<typename Sampler>
std::vector< ClusterVertex >
minimum_connection_probability_clustering(const ugraph_t & graph,
                                          Sampler & sampler,
                                          const size_t k,
                                          const size_t slack,
                                          const double rate,
                                          const probability_t p_low,
                                          Xorshift1024star & rnd) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > probabilities(n);
  std::vector< ugraph_vertex_t > uncovered_scratch(n);
  ConnectionCountsCache cccache(std::min(k, 500ul));
  size_t iteration = 0;
  probability_t p_curr = 1.0;
  ExponentialGuesser guesser(rate, p_low);
  size_t uncovered = n;
  size_t used_slack = 0;
  size_t best_clustering_iteration = 0;

  while (!guesser.stop()) {
    LOG_INFO(">>> Build clustering with p_curr=" << p_curr);
    cccache.cleanup();
    LOG_DEBUG(cccache.str());
    std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());
    uncovered = n;
    used_slack = 0;
    sampler.min_probability(graph, p_curr);

    // Build the clustering. Start the count from 1 because the
    // stopping condition is _inside_ the cycle
    for (size_t center_cnt = 1; center_cnt < k; center_cnt++) {
      assert(uncovered == count_uncovered(vinfo));
      ugraph_vertex_t center = pick_vertex_rnd(graph, rnd, cccache, uncovered_scratch, vinfo);
      vinfo[center].make_center(center);
      uncovered--;
      sampler.connection_probabilities_cache(graph, center, cccache, probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i=0; i<n; i++) {
        if (probabilities[i] >= p_curr) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
          } else if (!vinfo[i].is_covered() && vinfo[i].unreliable_probability() < probabilities[i]) {
            vinfo[i].unreliable_cover(center, probabilities[i]);
          }
        }
      }
      if (center_cnt + uncovered <= k + slack) {
        // Complete the clustering using the slack and break from the loop.
        for (ugraph_vertex_t i=0; i<n; i++) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].make_center(i);
            used_slack++;
            uncovered--;
          }
        }
        break;
      }
    }
    
    // Check if the clustering is valid
    if (uncovered == 0) {
      guesser.below();
      // this is a valid clustering, keep track of it
      best_clustering_iteration = iteration;
      for (ugraph_vertex_t i=0; i<n; i++) {
        valid_clustering[i] = vinfo[i];
      }
    } else {
      guesser.above();
    }
    
    LOG_INFO("Still " << uncovered << " nodes to cover");
    LOG_INFO("Cache hit rate: " << std::fixed << std::setprecision(2)
             << cccache.perc_hits() << "%");
    // update the probability
    p_curr = guesser.guess();
    iteration++;
  }

  EXPERIMENT_APPEND("algorithm-info", {{"used-slack", used_slack},
        {"p_curr", p_curr},
          {"best k-center iteration", best_clustering_iteration},});
  return valid_clustering;
  
}
