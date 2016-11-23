#include "sampler.hpp"
#include "logging.hpp"

namespace std {
  std::ostream & operator<<(std::ostream &os, CCSamplerThreadState & tstate) {
    os << "TState(" << tstate.rnd << ")";
    return os;
  }
}


void sample(const ugraph_t & g,
            CCSamplerThreadState & tstate)
{
  using namespace boost;
  
  BGL_FORALL_EDGES(e, g, ugraph_t) {
    const EdgeData & ed = g[e];
    tstate.edge_sample[ed.index] = tstate.rnd.next_double() <= ed.probability;
  }
}

void dfs(const ugraph_t & graph,
         const std::vector< bool > & smpl,
         std::vector< int > & components_map,
         std::vector< ugraph_vertex_t > & stack,
         ugraph_vertex_t root,
         int component_id) {
  using namespace boost;

  stack.clear();
  stack.push_back(root);
  components_map[root] = component_id;

  while(!stack.empty()) {
    ugraph_vertex_t u = stack.back();
    stack.pop_back();

    BGL_FORALL_OUTEDGES(u, e, graph, ugraph_t) {
      if (smpl[ graph[e].index ]) {
        ugraph_vertex_t v = target(e, graph);
        if (components_map[v] < 0) {
          components_map[v] = component_id;
          stack.push_back(v);
        }
      }
    }
  }
}

void connected_components(const ugraph_t & graph,
                          const std::vector< bool > & smpl,
                          std::vector< int > & components_map,
                          std::vector< ugraph_vertex_t > & stack) {
  using namespace boost;

  // initialize
  stack.clear();
  for (size_t i=0; i < components_map.size(); i++) {
    components_map[i] = -1;
  }
  int component_id = 0;
  BGL_FORALL_VERTICES(root, graph, ugraph_t) {
    if (components_map[root] < 0) {
      dfs(graph, smpl, components_map, stack, root, component_id++);
    }
  }
}

void CCSampler::min_probability(const ugraph_t & graph, probability_t prob) {
  const size_t n_samples = prob_to_samples(prob);
  sample_size(graph, n_samples);
  m_min_probability = (prob < m_min_probability)? prob : m_min_probability;
}

void CCSampler::sample_size(const ugraph_t & graph, size_t total_samples) {
  size_t new_samples = total_samples - m_samples.size();
  if (new_samples <= 0) {
    return;
  }
  LOG_INFO("Taking " << new_samples << " new samples");
  size_t start = m_samples.size();

  for (size_t i=0; i<new_samples; ++i) {
    // Build a new connected component vector in place
    m_samples.emplace_back(boost::num_vertices(graph), -1);
  }
  
#pragma omp parallel for default(none) shared(graph, start, new_samples)
  for (size_t i=start; i<start+new_samples; ++i) {
    auto tid = omp_get_thread_num();
    auto & tstate = m_thread_states[tid];
    sample(graph, tstate);
    auto & components = m_samples[i];
    connected_components(graph, tstate.edge_sample, components, tstate.stack);
  }
}


size_t CCSampler::connection_probabilities(const ugraph_t & graph,
                                           const ugraph_vertex_t from,
                                           std::vector< probability_t > & probabilities) {
  const size_t num_samples = m_samples.size();
  const size_t n = boost::num_vertices(graph);

  // Clear data structures
  std::fill(probabilities.begin(), probabilities.end(), 0.0);
  for (auto & tstate : m_thread_states) {
    std::fill(tstate.connection_counts.begin(), tstate.connection_counts.end(), 0);
  }
  
  // Accumulate, in parallel, the connection counts
#pragma omp parallel for default(none) shared(graph)
  for (size_t sample_idx=0; sample_idx < num_samples; sample_idx++) {
    auto tid = omp_get_thread_num();
    auto & tstate = m_thread_states[tid];
    
    const auto & smpl = m_samples[sample_idx];
    const int root_cc = smpl[from];
    for (size_t i=0; i< n; i++) {
      if (smpl[i] == root_cc) {
        tstate.connection_counts[i]++;
      }
    }
  }

  // Sum the partial counts together
  for (auto & tstate : m_thread_states) {
    for (size_t i=0; i< n; i++) {
      probabilities[i] += tstate.connection_counts[i];
    }
  }

  size_t cnt = 0;
  for (size_t i=0; i< n; i++) {
    probabilities[i] /= num_samples;
    if (probabilities[i] >= m_min_probability) {
      cnt++;
    }
  }

  return cnt;
}
