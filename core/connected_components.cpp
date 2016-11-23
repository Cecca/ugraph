#include "connected_components.hpp"

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
