#include "cc_sampler.hpp"
#include "logging.hpp"
#include "connected_components.hpp"

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


void CCSampler::min_probability(const ugraph_t & graph, probability_t prob) {
  const size_t n_samples = prob_to_samples(prob);
  sample_size(graph, n_samples);
  m_min_probability = (prob < m_min_probability)? prob : m_min_probability;
}

void CCSampler::sample_size(const ugraph_t & graph, size_t total_samples) {
  m_used_samples = total_samples;
  if (total_samples <= m_samples.size()) {
    LOG_INFO("Using " << m_used_samples << " (no new samples)");
    return;
  }
  size_t new_samples = total_samples - m_samples.size();
  LOG_INFO("Using " << m_used_samples << " (taking " << new_samples << " new)");
  size_t start = m_samples.size();

  for (size_t i=0; i<new_samples; ++i) {
    // Build a new connected component vector in place
    m_samples.emplace_back(boost::num_vertices(graph), -1);
  }
  LOG_DEBUG("Allocated the new samples");
  
#pragma omp parallel for default(none) shared(graph, start, new_samples)
  for (size_t i=start; i<start+new_samples; ++i) {
    auto tid = omp_get_thread_num();
    auto & tstate = m_thread_states[tid];
    //sample(graph, tstate);
    auto & components = m_samples[i];
    union_find(graph, tstate.rnd, tstate.ranks, components);
    //connected_components(graph, tstate.edge_sample, components, tstate.stack);
  }
}


size_t CCSampler::connection_probabilities(const ugraph_t & graph,
                                           const ugraph_vertex_t from,
                                           std::vector< probability_t > & probabilities) {
  const size_t num_samples = m_used_samples;
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
    auto & connection_counts = m_thread_states[tid].connection_counts;
    
    const auto & smpl = m_samples[sample_idx];
    const int root_cc = smpl[from];
    for (size_t i=0; i< n; i++) {
      if (smpl[i] == root_cc) {
        connection_counts[i]++;
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

size_t CCSampler::connection_probabilities_cache(const ugraph_t & graph,
                                                 const ugraph_vertex_t from,
                                                 ConnectionCountsCache & cccache,
                                                 std::vector< probability_t > & probabilities) {
  const size_t num_samples = m_used_samples;
  const size_t n = boost::num_vertices(graph);

  // Clear data structures
  std::fill(probabilities.begin(), probabilities.end(), 0.0);
  for (auto & tstate : m_thread_states) {
    std::fill(tstate.connection_counts.begin(), tstate.connection_counts.end(), 0);
  }

  ConnectionCountsCacheElement & ccc_elem = cccache.get_or_new(from, n);

  const size_t starting_sample = ccc_elem.num_samples;
  LOG_DEBUG("Node " << from << ": starting from sample " <<
            starting_sample << " of " << num_samples << " used");
   
  // Accumulate, in parallel, the connection counts
#pragma omp parallel for default(none) shared(graph)
  for (size_t sample_idx=starting_sample; sample_idx < num_samples; sample_idx++) {
    auto tid = omp_get_thread_num();
    auto & connection_counts = m_thread_states[tid].connection_counts;
    
    const auto & smpl = m_samples[sample_idx];
    const int root_cc = smpl[from];
    for (size_t i=0; i< n; i++) {
      // Doing this way allows GCC (with -O3) to vectorize (probably)
      // the code. The if statement is slighlty slower.
      connection_counts[i] += ((smpl[i] == root_cc)? 1 : 0);
    }
  }

  // Sum the partial counts together
  for (auto & tstate : m_thread_states) {
    for (size_t i=0; i< n; i++) {
      ccc_elem.counts[i] += tstate.connection_counts[i];
    }
  }

  // Take the maximum because we may get back when doing binary
  // search, and we don't want to screw up estimates.
  ccc_elem.num_samples = std::max(num_samples, ccc_elem.num_samples);
  
  size_t cnt = 0;
  for (size_t i=0; i< n; i++) {
    probabilities[i] = ccc_elem.counts[i] / ((double) ccc_elem.num_samples);
    if (probabilities[i] >= m_min_probability) {
      cnt++;
    }
  }

  return cnt;
}


size_t CCSampler::connection_probabilities(const ugraph_t & graph,
                                           const ugraph_vertex_t from,
                                           const std::vector< ugraph_vertex_t > & targets,
                                           std::vector< probability_t > & probabilities) {
  const size_t num_samples = m_samples.size();

  // Clear data structures
  std::fill(probabilities.begin(), probabilities.end(), 0.0);
  for (auto & tstate : m_thread_states) {
    std::fill(tstate.connection_counts.begin(), tstate.connection_counts.end(), 0);
  }
  
  // Accumulate, in parallel, the connection counts
#pragma omp parallel for default(none) shared(graph, targets)
  for (size_t sample_idx=0; sample_idx < num_samples; sample_idx++) {
    auto tid = omp_get_thread_num();
    auto & tstate = m_thread_states[tid];
    
    const auto & smpl = m_samples[sample_idx];
    const int root_cc = smpl[from];
    for (const ugraph_vertex_t t : targets) {
      if (smpl[t] == root_cc) {
        tstate.connection_counts[t]++;
      }
    }
  }

  // Sum the partial counts together
  for (auto & tstate : m_thread_states) {
    for (const ugraph_vertex_t t : targets) {
      probabilities[t] += tstate.connection_counts[t];
    }
  }

  size_t cnt = 0;
  for (const ugraph_vertex_t t : targets) {
    probabilities[t] /= num_samples;
    if (probabilities[t] >= m_min_probability) {
      cnt++;
    }
  }

  return cnt;
}


probability_t CCSampler::connection_probability(const ugraph_t & graph,
                                                const std::vector< ugraph_vertex_t > & vertices) {
  const size_t num_samples = m_samples.size();

  const ugraph_vertex_t vertices_root = vertices[0];
  
  size_t count = 0;

  // Accumulate the counts in parallel
#pragma omp parallel for reduction(+:count) default(none) shared(vertices)
  for(size_t sample_idx=0; sample_idx<num_samples; ++sample_idx) {
    const auto & component = m_samples[sample_idx];
    const int component_id = component[vertices_root];
    bool connected = true;
    for(const ugraph_vertex_t v : vertices) {
      if (component[v] != component_id) {
        connected = false;
        break;
      }
    }
    if (connected) {
      count++;
    }
  }

  return count / ((double) num_samples);
}
