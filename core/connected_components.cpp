#include "connected_components.hpp"
#include "logging.hpp"

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

void union_find(const ugraph_t & graph,
                Xorshift1024star & rnd,
                std::vector< size_t > & ranks,
                std::vector< int > & roots) {
  using namespace boost;
  
  std::fill(ranks.begin(), ranks.end(), 0);

  // Each node starts in it own connected component
  for (ugraph_vertex_t u=0; u < ranks.size(); ++u) {
    roots[u] = u;
  }

  BGL_FORALL_EDGES(e, graph, ugraph_t) {
    // sample the edge!
    if (rnd.next_double() <= graph[e].probability) {
      ugraph_vertex_t u = source(e, graph);
      ugraph_vertex_t v = target(e, graph);
      
      // find roots
      // TODO apply path compression
      while (u != roots[u]) { u = roots[u]; }
      while (v != roots[v]) { v = roots[v]; }

      if (u != v) {
        // the nodes are in different connected components
        size_t
          rank_u = ranks[u],
          rank_v = ranks[v];
        if (rank_u < rank_v) {
          roots[u] = v;
        } else if (rank_u > rank_v) {
          roots[v] = u;
        } else {
          roots[v] = u;
          ranks[u]++;
        }
      }
    }
  }

  // Compress the components
  for (ugraph_vertex_t u=0; u < ranks.size(); ++u) {
    ugraph_vertex_t r = u;
    while ( r != roots[r] ) { r = roots[r]; }
    roots[u] = r;
  }
}
