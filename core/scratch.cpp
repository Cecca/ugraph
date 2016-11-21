#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "sampler.hpp"
#include "logging.hpp"
#include "sequential_clustering.hpp"


size_t prob_to_samples(probability_t prob, double epsilon, double delta) {
  return 1/(epsilon*epsilon*prob) * log(1/delta);
}

int main(int argc, char**argv) {
  std::string graph_path(argv[1]);
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_INFO("Num nodes " << boost::num_vertices(graph) <<
           " num edges " << boost::num_edges(graph));

  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");

  Xorshift1024star rnd(1234);
  CCSampler sampler(graph, 0.1, 0.01, prob_to_samples, 123, omp_threads);
  // sampler.log_states();
  // sampler.min_probability(graph, 0.5);
  // LOG_INFO("==================");
  // sampler.log_states();

  // LOG_INFO("Computing connection probabilities");
  // ugraph_vertex_t root = 1;
  // std::vector< probability_t > probabilities(boost::num_vertices(graph), 0.0);

  // size_t reliably_estimated = sampler.connection_probabilities(graph, root, probabilities);

  // for (size_t i=0; i < probabilities.size(); i++) {
  //   std::cout << i << ":" << probabilities[i] << " ";
  //   if (i!= 0 && i % 10 == 0) {
  //     std::cout << std::endl;
  //   }
  // }
  // std::cout << std::endl;

  // LOG_INFO("Of which " << reliably_estimated << " reliably estimated");

  sequential_cluster(graph, sampler, 128, 0, 0.9, 0.001, rnd);
}
