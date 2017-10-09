#include "types.hpp"
#include "io.hpp"
#include "cc_sampler.hpp"

boost::program_options::variables_map
parse_args(int argc, char** argv)
{
  using namespace boost;
  namespace po = boost::program_options;
  
  po::options_description desc("Clustering evaluation");
  desc.add_options()
    ("help", "produce help message")
    ("debug", "print debug output")
    ("trace", "print trace output")
    ("epsilon", po::value<double>()->default_value(0.1),
     "tolerated absolute error")
    ("delta", po::value<double>()->default_value(0.01),
     "error probability")
    ("graph", po::value<std::string>(),
     "file containing graph data");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    abort();
  }
  if (vm.count("trace")) {
    logging::set_level(logging::Level::Trace);
  } else if (vm.count("debug")) {
    logging::set_level(logging::Level::Debug);
  }
  return vm;
}


int main(int argc, char **argv) {
  auto vm = parse_args(argc, argv);

  std::random_device rd;
  uint64_t seed = rd();
  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");

  std::string path = vm["graph"].as<std::string>();
  
  ugraph_t graph;
  read_edge_list(graph, path);

  double
    epsilon = vm["epsilon"].as<double>(),
    delta = vm["delta"].as<double>();
  auto prob_to_samples = [epsilon, delta](double p) {return 1/(epsilon*epsilon*p) * log(1/delta);};
  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);
  sampler.min_probability(graph, 0.01);

  using namespace boost;

  std::vector<double> probabilities(num_vertices(graph));
  double probability_sum = 0.0;
  BGL_FORALL_VERTICES(v, graph, ugraph_t) {
    sampler.connection_probabilities(graph, v, probabilities);
    for (const double p : probabilities) {
      probability_sum += p;
    }
  }
  size_t n = num_vertices(graph);
  // We don't divide the denominator by 2 because the probability_sum
  // counts twice each pair
  double average_probability = probability_sum / (n*(n-1));

  LOG_INFO("Average connection probability is " << average_probability);
}
