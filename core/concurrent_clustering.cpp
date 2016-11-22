#include "concurrent_clustering.hpp"

size_t select_centers(std::vector< ConcurrentClusterVertex > & vinfo,
                      const probability_t prob,
                      Xorshift1024star & rnd) {
  const size_t n = vinfo.size();
  size_t cnt = 0;
  for (size_t v=0; v<n; v++) {
    if (!vinfo[v].is_covered() && rnd.next_double() <= prob) {
      cnt++;
      vinfo[v].make_center(v);
    }
  }
  return cnt;
}

std::vector< ConcurrentClusterVertex > concurrent_cluster(const ugraph_t & graph,
                                                          CCSampler & sampler,
                                                          const size_t batch,
                                                          const probability_t p_low,
                                                          Xorshift1024star & rnd,
                                                          ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ConcurrentClusterVertex > vinfo(n);
  probability_t p_curr = 1.0;
  size_t uncovered = n;
  // (probability, (center, node))
  std::vector< std::pair< probability_t, std::pair< ugraph_vertex_t, ugraph_vertex_t > > >
    probabilities_pq;
  
  while(uncovered > 0) {
    probability_t selection_prob = ((double) batch) / uncovered;
    select_centers(vinfo, selection_prob, rnd);
    
  }
  
  return vinfo;
}
