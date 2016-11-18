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
  sampler.set_sample_size(graph, 4);
  LOG_INFO("==================");
  sampler.log_states();
}
