#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"

int main(int argc, char**argv) {
  std::string graph_path(argv[1]);
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  std::cout << "Num nodes " << boost::num_vertices(graph)
            << " num edges " << boost::num_edges(graph)
            << std::endl;
}
