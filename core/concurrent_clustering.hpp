#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include <boost/functional/hash.hpp>

void mark_reachable(const ugraph_t & graph,
                    const ugraph_vertex_t root,
                    std::vector< bool > & flags,
                    std::vector< ugraph_vertex_t > stack) {
  using namespace boost;
  stack.clear();
  stack.push_back(root);
  flags[root] = true;
  while(!stack.empty()) {
    ugraph_vertex_t v = stack.back();
    stack.pop_back();
    BGL_FORALL_OUTEDGES(v, e, graph, ugraph_t) {
      ugraph_vertex_t u = target(e, graph);
      if (!flags[u]) {
        flags[u] = true;
        // Enqueue u only if it is not already marked, otherwise we
        // already know the nodes reachable from there from a previous visit.
        stack.push_back(u);
      }
    }
  }
}

/**
 * Select centers according to the given probability. The selection
 * proceeds iteratively, until the selected batch of centers can reach
 * the target number of nodes.
 */
size_t select_centers(const ugraph_t & graph,
                      std::vector< ClusterVertex > & vinfo,
                      std::vector< ugraph_vertex_t > & centers,
                      const probability_t prob,
                      const size_t target,
                      Xorshift1024star & rnd,
                      std::vector< bool > & flags,
                      std::vector< ugraph_vertex_t > & stack) {
  const size_t n = vinfo.size();
  size_t tentative = 0;
  const size_t max_tentatives = 1024;
  size_t reachable = 0;

  centers.clear();

  while (reachable < target) {
    if (tentative > 0) {
      LOG_WARN("Could reach only " << reachable << "/" << target
                                   << " nodes from the " << centers.size()
                                   << " selected centers, doing tentative "
                                   << tentative << "/" << max_tentatives);
    }
    if (tentative >= max_tentatives) {
      throw std::logic_error("Exceeded maximum center selection tentatives");
    }
    centers.clear();
    reachable = 0;
    std::fill(flags.begin(), flags.end(), false);
    for (size_t v=0; v<n; v++) {
      if (!vinfo[v].is_covered() && rnd.next_double() <= prob) {
        assert(!vinfo[v].is_center());
        centers.push_back(v);
      }
    }
    for (auto c : centers) {
      mark_reachable(graph, c, flags, stack);
    }
    for(size_t i = 0; i<n; i++) {
      if(flags[i] && !vinfo[i].is_covered()) {
        reachable++;
      }
    }
    tentative++;
  }

  for (auto c : centers) {
    vinfo[c].make_center(c);
  }
  return centers.size();
}

template<typename Sampler>
std::vector< ClusterVertex > concurrent_cluster(const ugraph_t & graph,
                                                Sampler & sampler,
                                                const size_t batch,
                                                const probability_t p_low,
                                                Xorshift1024star & rnd,
                                                ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);

  // ----------------------------------
  // Data structures initialization
  
  std::vector< ClusterVertex > vinfo(n); // Vertex information
  std::vector< probability_t > probabilities(n);   // scratch vector for connection probabilities
  std::vector< ugraph_vertex_t > active_centers;   // scratch vector for IDs of active centers
  active_centers.reserve(2*batch);
  std::vector< bool > potential_cover_flags(n); // vector of flags to count the number of distinct nodes that can be covered
  std::vector< ugraph_vertex_t > stack(n);
  std::vector< std::pair< probability_t, std::pair< ugraph_vertex_t, ugraph_vertex_t > > >
    probabilities_pq; // priority queue of center--node by decreasing probability: (probability, (center, node))

  // ----------------------------------
  // Algorithm
  
  probability_t p_curr = 1.0;   // current probability threshold
  size_t uncovered = n;         // Count of uncovered nodes
  
  while(uncovered > 0) {
    probability_t selection_prob = ((double) batch) / uncovered;
    LOG_INFO("Selecting centers with probability " << selection_prob);
    size_t num_selected = select_centers(graph, vinfo, active_centers, selection_prob,
                                         std::max((size_t) 1, uncovered/2), rnd, potential_cover_flags, stack);
    assert(num_selected <= uncovered);
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
      size_t count_usable = num_selected;
      for(bool f : potential_cover_flags) { if(f) count_usable++; }
      if (count_usable >= uncovered / 2) {
        break;
      } else {
        LOG_INFO("Only " << count_usable << " 'usable nodes' (" << uncovered / 2
                         << " needed, p_curr=" << p_curr << ", p_max=" << max_p
                         << ")");
        p_curr = p_curr / 2;
        if (p_curr <= p_low) {
          throw std::logic_error("Could not find a clustering high enough connection probabilities");
        }
      }
    }
    
    std::make_heap(probabilities_pq.begin(), probabilities_pq.end());
    size_t covered=num_selected;
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

template<typename Sampler>
void shrink_clustering(const ugraph_t & graph,
                       Sampler & sampler,
                       const size_t target,
                       std::vector< ClusterVertex > & vinfo) {
  const size_t n = boost::num_vertices(graph);
  std::vector < ugraph_vertex_t > centers;
  std::unordered_map< ugraph_vertex_t, ugraph_vertex_t > center_mapping;
  std::unordered_map< ugraph_vertex_t, size_t > super_centers;
  std::vector< std::vector< probability_t > > probabilities(target); // Scratch vectors for probabilities
  for(size_t i=0; i<target; i++) {
    std::vector<probability_t> vec(n);
    probabilities[i] = vec;
  }
  
  probability_t p_alg = 1.0;
  size_t num_clusters = 0;
  for(ugraph_vertex_t v=0; v<n; ++v){
    p_alg = std::min(vinfo[v].probability(), p_alg);
    if (vinfo[v].is_center()) {
      num_clusters++;
      centers.push_back(v);
    }
  }

  LOG_DEBUG("Starting the guessing");
  probability_t guess = p_alg;
  while(center_mapping.size() < num_clusters) {
    LOG_INFO(">>>> guess " << guess);
    sampler.min_probability(graph, guess);
    center_mapping.clear();
    super_centers.clear();
    size_t super_center_idx = 0;
    for(ugraph_vertex_t c1 : centers) {
      if (center_mapping.count(c1) == 0) {
        super_centers[c1] = super_center_idx++;
        LOG_DEBUG("Getting probabilities " << c1 << " (" << super_centers[c1] << ")");
        auto & ps = probabilities.at(super_centers[c1]);
        assert(ps.size()==n);
        LOG_DEBUG("Got probabilities");
        sampler.connection_probabilities(graph, c1, ps);
        LOG_DEBUG("Computed probabilities");
        for (ugraph_vertex_t c2 : centers) {
          if (center_mapping.count(c2) == 0 && ps[c2] >= guess) {
            center_mapping[c2] = c1;
            LOG_DEBUG("Remapped " << c2 << " to " << c1);
          }
        }
      }
      if (super_center_idx >= target) {
        break;
      }
    }
    LOG_INFO("Center mapping of size " <<
             center_mapping.size() << "/" << centers.size());
    guess /= 2;
  }

  std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());

  for (const auto & center_pair : super_centers) {
    const ugraph_vertex_t c = center_pair.first;
    const auto & probs = probabilities[super_centers[c]];
    for (ugraph_vertex_t v=0; v<n; v++) {
      if (!vinfo[v].is_covered() || (vinfo[v].is_covered() && vinfo[v].probability() < probs[v])) {
        LOG_DEBUG("Cover " << v << " from " << c << " with p=" << probs[v]);
        vinfo[v].cover(c, probs[v]);
      }
    }
  }
}
