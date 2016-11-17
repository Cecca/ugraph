#include "sampler.hpp"

void sample(ugraph_t & g,
            std::vector< bool > & edge_sample,
            Xorshift1024star & rnd)
{
  using namespace boost;
  
  BGL_FORALL_EDGES(e, g, ugraph_t) {
    EdgeData & ed = g[e];
    edge_sample[ed.index] = rnd.next_double() <= ed.probability;
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


bool next_undefined(const ugraph_t & graph,
                    const std::vector< int > & components_map,
                    ugraph_vertex_t & out) {
  using namespace boost;
  
  BGL_FORALL_VERTICES(v, graph, ugraph_t) {
    if (components_map[v] < 0) {
      out = v;
      return true;
    }
  }
  return false;
}

void connected_components(const ugraph_t & graph,
                          const std::vector< bool > & smpl,
                          std::vector< int > & components_map,
                          std::vector< ugraph_vertex_t > & stack) {
  using namespace boost;

  // initialize
  for (size_t i=0; i < components_map.size(); i++) {
    components_map[i] = -1;
  }
  int component_id = 0;
  ugraph_vertex_t root;
  while (next_undefined(graph, components_map, root)) {
    dfs(graph, smpl, components_map, stack, root, component_id++);
  }
}


void CCSampler::set_sample_size(ugraph_t & graph, size_t total_samples) {
  size_t new_samples = total_samples - m_samples.size();
  if (new_samples <= 0) {
    return;
  }
  size_t start = m_samples.size();

  for (size_t i=0; i<new_samples; ++i) {
    // Build a new connected component vector in place
    m_samples.emplace_back(boost::num_vertices(graph), -1);
  }

  auto tid = omp_get_thread_num();
  stack_t stack;
  component_vector_t components;
  edge_sample_t edge_sample;

#pragma omp parallel for default(none) shared(std::cout, graph, start, new_samples) private(stack, edge_sample, tid, components)
  for (size_t i=start; i<start+new_samples; ++i) {
    stack = s_stacks[tid];
    edge_sample = s_edge_samples[tid];
    sample(graph, edge_sample, m_rnds[tid]);
    std::cout << "HERE" << std::endl;
    components = m_samples[i];
    connected_components(graph, edge_sample, components, stack);
  }
}
