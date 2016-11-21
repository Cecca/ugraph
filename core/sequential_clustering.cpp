#include "sequential_clustering.hpp"

size_t count_uncovered(const std::vector< SequentialClusterVertex > & vinfo) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered()) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const std::vector< SequentialClusterVertex > & vinfo) {
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered()) {
      return i;
    }
  }
  throw std::logic_error("No uncovered node to select");
}

std::vector< SequentialClusterVertex > sequential_cluster(const ugraph_t & graph,
                                                          CCSampler & sampler,
                                                          const size_t k,
                                                          const size_t slack,
                                                          const double rate,
                                                          const probability_t p_low) {
  const size_t n = boost::num_vertices(graph);
  std::vector< SequentialClusterVertex > vinfo(n);
  std::vector< probability_t > probabilities(n);
  probability_t p_curr = 1.0;
  size_t uncovered = n;

  while (p_curr > p_low) {
    LOG_INFO("Build clustering with p_curr=" << p_curr);
    std::fill(vinfo.begin(), vinfo.end(), SequentialClusterVertex());
    uncovered = n;
    sampler.min_probability(graph, p_curr);

    // Build the clustering
    for (size_t center_cnt = 0; center_cnt < k; center_cnt++) {
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
        // Return the clustering
         for (ugraph_vertex_t i=0; i<n; i++) {
           if (!vinfo[i].is_covered()) {
             vinfo[i].make_center(i);
           }
         }
         return vinfo;
      }
    }

    LOG_INFO("Still " << uncovered << " nodes to cover");
    // update the probability
    p_curr *= rate;
  }

  throw std::logic_error("p_curr < p_low");
}
