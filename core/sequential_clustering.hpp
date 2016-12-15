#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"

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

template<typename Sampler>
std::vector< ClusterVertex > sequential_cluster(const ugraph_t & graph,
                                                Sampler & sampler,
                                                const size_t k,
                                                const size_t slack,
                                                const double rate,
                                                const probability_t p_low,
                                                ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< probability_t > probabilities(n);
  probability_t p_curr = 1.0;
  size_t uncovered = n;
  size_t used_slack = 0;

  while (uncovered > 0 && p_curr > p_low) {
    LOG_INFO("Build clustering with p_curr=" << p_curr);
    std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());
    uncovered = n;
    used_slack = 0;
    sampler.min_probability(graph, p_curr);

    // Build the clustering. Start the count from 1 because the
    // stopping condition is _inside_ the cycle
    for (size_t center_cnt = 1; center_cnt < k; center_cnt++) {
      assert(uncovered == count_uncovered(vinfo));
      ugraph_vertex_t center = pick_vertex(graph, vinfo);
      vinfo[center].make_center(center);
      uncovered--;
      sampler.connection_probabilities(graph, center, probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i=0; i<n; i++) {
        if (probabilities[i] >= p_curr) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
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
    
    LOG_INFO("Still " << uncovered << " nodes to cover");
    // update the probability
    p_curr *= rate;
  }

  if (uncovered == 0) {
    experiment.append("algorithm-info", {{"used-slack", used_slack},
          {"p_curr", p_curr}});
    return vinfo;
  } else {
    throw std::logic_error("p_curr < p_low");
  }
}
