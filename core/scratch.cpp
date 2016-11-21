#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "sampler.hpp"
#include "logging.hpp"


int main(int argc, char**argv) {
  std::string graph_path(argv[1]);
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_INFO("Num nodes " << boost::num_vertices(graph) <<
           " num edges " << boost::num_edges(graph));

  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");
  
  CCSampler sampler(graph, 123, omp_threads);
  sampler.log_states();
  sampler.sample_size(graph, 1024);
  LOG_INFO("==================");
  sampler.log_states();

  LOG_INFO("Computing connection probabilities");
  ugraph_vertex_t root = 1278;
  std::vector< probability_t > probabilities(boost::num_vertices(graph), 0.0);

  sampler.connection_probabilities(graph, root, probabilities);

  for (size_t i=0; i < probabilities.size(); i++) {
    std::cout << i << ":" << probabilities[i] << " ";
    if (i!= 0 && i % 10 == 0) {
      std::cout << std::endl;
    }
  }
}
