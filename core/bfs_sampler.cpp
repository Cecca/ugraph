#include "bfs_sampler.hpp"
#include "logging.hpp"

namespace std {
  std::ostream & operator<<(std::ostream &os, BfsSamplerThreadState & tstate) {
    os << "TState(" << tstate.rnd << ")";
    return os;
  }
}

const size_t g_infinite_distance = std::numeric_limits<size_t>::max();

void bfs(const ugraph_t & graph,
         const std::vector< bool > & smpl,
         std::vector< size_t > & distance_map,
         FixedCapacityQueue<ugraph_vertex_t> & queue,
         const ugraph_vertex_t root,
         const size_t max_dist) {
  queue.clear();
  queue.push(root);
  std::fill(distance_map.begin(), distance_map.end(), g_infinite_distance);
  distance_map[root] = 0;

  while(!queue.empty()) {
    ugraph_vertex_t v = queue.pop();
    size_t v_dist = distance_map[v];
    size_t new_dist = v_dist+1;
    if (new_dist <= max_dist) {
      BGL_FORALL_OUTEDGES(v, e, graph, ugraph_t) {
        if (smpl[ graph[e].index ]) {
          ugraph_vertex_t u = target(e, graph);
          size_t u_dist = distance_map[u];
          if (new_dist < u_dist) {
            distance_map[u] = v_dist+1;
            queue.push(u);
          }
        }
      }
    }
  }
  
}

void sample(const ugraph_t & g,
            BfsSamplerThreadState & tstate,
            std::vector<bool> & edge_sample)
{
  using namespace boost;
  
  BGL_FORALL_EDGES(e, g, ugraph_t) {
    const EdgeData & ed = g[e];
    edge_sample[ed.index] = tstate.rnd.next_double() <= ed.probability;
  }
}

void BfsSampler::min_probability(const ugraph_t & graph, probability_t prob) {
  const size_t n_samples = prob_to_samples(prob);
  sample_size(graph, n_samples);
  m_min_probability = (prob < m_min_probability)? prob : m_min_probability;
}


void BfsSampler::sample_size(const ugraph_t & graph, size_t total_samples) {
  size_t new_samples = total_samples - m_samples.size();
  if (new_samples <= 0) {
    return;
  }
  LOG_INFO("Taking " << new_samples << " new samples");
  size_t start = m_samples.size();

  for (size_t i=0; i<new_samples; ++i) {
    // Build a new edge sample vector in place
    m_samples.emplace_back(boost::num_edges(graph), false);
  }
  
#pragma omp parallel for default(none) shared(graph, start, new_samples)
  for (size_t i=start; i<start+new_samples; ++i) {
    auto tid = omp_get_thread_num();
    auto & tstate = m_thread_states[tid];
    auto & edge_sample = m_samples[i];
    sample(graph, tstate, edge_sample);
  }
}

size_t BfsSampler::connection_probabilities(const ugraph_t & graph,
                                            const ugraph_vertex_t from,
                                            const std::vector< ugraph_vertex_t > & targets,
                                            std::vector< probability_t > & probabilities) {
  return connection_probabilities(graph, from, probabilities);
}

size_t BfsSampler::connection_probabilities(const ugraph_t & graph,
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
    bfs(graph, smpl, tstate.distance_vector, tstate.queue, from, m_max_dist);
    for (size_t i=0; i< n; i++) {
      if (tstate.distance_vector[i] < g_infinite_distance) {
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
