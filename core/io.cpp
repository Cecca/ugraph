#include "io.hpp"

inline
boost::graph_traits<ugraph_t>::vertex_descriptor
get_or_add(ugraph_t & g,
           std::unordered_map<std::string,
           boost::graph_traits<ugraph_t>::vertex_descriptor>& name_map,
           std::string& name) {
  using namespace boost;
  typedef graph_traits<ugraph_t>::vertex_descriptor vertex_t;

  vertex_t id;
  if (name_map.find(name) == name_map.end()) {
    id = add_vertex(VertexData(name), g);
    name_map[name] = id;
  } else {
    id = name_map[name];
  }
  return id;
}


void read_edge_list(ugraph_t & g, std::string path) {
  using namespace boost;

  // set to filter out duplicate edges
  std::set<std::pair<ugraph_vertex_t, ugraph_vertex_t>> dedup;

  // map from names to ids
  std::unordered_map<std::string, ugraph_vertex_t> names_map;
  
  EdgeDataFactory edgef;
  
  std::ifstream input(path);
  if (!input.good()) {
    throw std::runtime_error("File not found!");
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line[0] != '#' && line.size() > 0) { // filter comments and empty lines
      std::vector<std::string> tokens;
      split(tokens, line, is_any_of("\t "));

      ugraph_vertex_t src = get_or_add(g, names_map, tokens[0]);
      ugraph_vertex_t dst = get_or_add(g, names_map, tokens[1]);
      probability_t edge_probability =
        (tokens.size() > 2)? std::stof(tokens[2]) : 1.0;

      auto dedup_key =
        (src <= dst)? std::make_pair(src,dst) : std::make_pair(dst,src);
      
      if (dedup.find(dedup_key) == dedup.end()) {
        dedup.insert(dedup_key);
        add_edge(src, dst, edgef.build(edge_probability), g);
      }
    }
  }
}
