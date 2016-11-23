#include "concurrent_clustering.hpp"

size_t select_centers(std::vector< ClusterVertex > & vinfo,
                      std::vector< ugraph_vertex_t > & centers,
                      const probability_t prob,
                      Xorshift1024star & rnd) {
  const size_t n = vinfo.size();
  centers.clear();
  size_t cnt = 0;
  for (size_t v=0; v<n; v++) {
    if (!vinfo[v].is_covered() && rnd.next_double() <= prob) {
      cnt++;
      vinfo[v].make_center(v);
      centers.push_back(v);
    }
  }
  return cnt;
}

std::vector< ClusterVertex > concurrent_cluster(const ugraph_t & graph,
                                                          CCSampler & sampler,
                                                          const size_t batch,
                                                          const probability_t p_low,
                                                          Xorshift1024star & rnd,
                                                          ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n); // Vertex information
  std::vector< probability_t > probabilities(n);   // scratch vector for connection probabilities
  std::vector< ugraph_vertex_t > active_centers;   // scratch vector for IDs of active centers
  active_centers.reserve(2*batch);
  std::vector< bool > potential_cover_flags(n); // vector of flags to count the number of distinct nodes that can be covered
  std::vector< std::pair< probability_t, std::pair< ugraph_vertex_t, ugraph_vertex_t > > >
    probabilities_pq; // priority queue of center--node by decreasing probability: (probability, (center, node))
  
  probability_t p_curr = 1.0;   // current probability threshold
  size_t uncovered = n;         // Count of uncovered nodes
  
  while(uncovered > 0) {
    probability_t selection_prob = ((double) batch) / uncovered;
    LOG_INFO("Selecting centers with probability " << selection_prob);
    size_t num_selected = select_centers(vinfo, active_centers, selection_prob, rnd);
    uncovered -= num_selected;
    LOG_INFO("Still " << uncovered << "/" << n << " nodes after center selection (" << num_selected << " selected)");
    if (uncovered == 0) { break; } // early exit if all the uncovered nodes are turned into centers

    while (true) { // Termination condition at the end
      probability_t max_p = 0.0;
      sampler.min_probability(graph, p_curr);
      probabilities_pq.clear();
      std::fill(potential_cover_flags.begin(), potential_cover_flags.end(), false);
      for(ugraph_vertex_t c : active_centers) {
        sampler.connection_probabilities(graph, c, probabilities);
        for(ugraph_vertex_t v=0; v<n; v++) {
          if (!vinfo[v].is_covered()) {
            probability_t p = probabilities[v];
            max_p = (p>max_p)? p : max_p;
            if (p >= p_curr) {
              probabilities_pq.push_back(std::make_pair(p, std::make_pair(c, v)));
              potential_cover_flags[v] = true;
            }
          }
        }
      }
      size_t count_usable = 0;
      for(bool f : potential_cover_flags) { if(f) count_usable++; }
      if (count_usable >= uncovered / 2) {
        break;
      } else {
        LOG_INFO("Only " << count_usable << " 'usable nodes' (" << uncovered / 2
                         << " needed, p_curr=" << p_curr << ", p_max=" << max_p
                         << ")");
        p_curr = p_curr / 2;
      }
    }
    
    std::make_heap(probabilities_pq.begin(), probabilities_pq.end());
    size_t covered=0;
    while (covered<uncovered/2) {
      assert(!probabilities_pq.empty());
      std::pop_heap(probabilities_pq.begin(), probabilities_pq.end());
      auto to_cover = probabilities_pq.back();
      probabilities_pq.pop_back();
      if (!vinfo[to_cover.second.second].is_covered()) {
        vinfo[to_cover.second.second].cover(to_cover.second.first, to_cover.first);
        uncovered -= 1;
        covered++;
      }
    }
    LOG_INFO("After covering, still " << uncovered << "/" << n << " uncovered nodes");
  }
  
  return vinfo;
}
