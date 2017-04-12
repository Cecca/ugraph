#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "logging.hpp"
#include "git_info.hpp"
#include "average_probability_clustering.hpp"
#include "experiment_reporter.hpp"
#include "scores.hpp"


boost::program_options::variables_map
parse_args(int argc, char** argv)
{
  using namespace boost;
  namespace po = boost::program_options;

  po::options_description desc("Simple clustering");
  desc.add_options()
    ("help", "produce help message")
    ("revision", "get the git revision of the code")
    ("debug", "print debug output")
    ("trace", "print trace output")
    ("graph", po::value<std::string>(),
     "input graph")
    ("target,k", po::value<size_t>(),
     "desired number of clusters")
    ("rate,r", po::value<double>()->default_value(0.5),
     "decrease rate")
    ("depth", po::value<size_t>(), "BFS depth")
    ("epsilon", po::value<double>()->default_value(0.1),
     "tolerated absolute error")
    ("delta", po::value<double>()->default_value(0.01),
     "error probability")
    ("theory-samples-fraction", po::value<double>()->default_value(0.1),
     "Fraction of samples to be used with respect to the theory-defined formula")
    ("seed", po::value<uint64_t>(),
     "seed for random generator")
    ("fast-scores", "compute only fast-to-compute scores");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("revision")) {
    std::cout << g_GIT_SHA1 << std::endl;
    std::exit(0);
  }
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::exit(0);
  }
  if (vm.count("trace")) {
    logging::set_level(logging::Level::Trace);
  } else if (vm.count("debug")) {
    logging::set_level(logging::Level::Debug);
  }
  return vm;
}

void check_num_components(const ugraph_t & graph, const size_t target) {
  auto component_map = boost::make_vector_property_map<int>(boost::get(boost::vertex_index, graph));
  size_t num_components = boost::connected_components(graph, component_map);
  if (target < num_components) {
    LOG_ERROR("The target size ("
              << target
              << " + "
              << ") is smaller than the number of connected components ("
              << num_components << "): the algorithm can't terminate");
    throw std::logic_error("Target size too small");
  }
}

void add_clustering_info(const ugraph_t &graph,
                         const std::vector<ClusterVertex> &vinfo,
                         const std::string & table_name,
                         ExperimentReporter &exp) {
  size_t n = vinfo.size();
  for (ugraph_vertex_t v = 0; v < n; v++) {
    ugraph_vertex_t center = vinfo[v].center();
    exp.append(table_name.c_str(), {{"id", v},
                              {"center", center},
                              {"label", graph[v].label},
                              {"center label", graph[center].label},
                              {"probability", vinfo[v].probability()}});
  }
}

int main(int argc, char**argv) {
  auto args = parse_args(argc, argv);

  std::string graph_path(args["graph"].as<std::string>());

  uint64_t seed;
  if (args.count("seed")) {
    seed = args["seed"].as<uint64_t>();
  } else {
    std::random_device rd;
    seed = rd();
    LOG_INFO("Using random seed " << seed);
  }

  double
    epsilon = args["epsilon"].as<double>(),
    delta = args["delta"].as<double>(),
    rate = args["rate"].as<double>(),
    theory_samples_fraction = args["theory-samples-fraction"].as<double>(),
    p_low = 0.0001;

  size_t k = args["target"].as<size_t>();

  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");  

  ExperimentReporter exp;
  exp.tag("algorithm", std::string("average-probability"));
  exp.tag("input", graph_path);
  exp.tag("epsilon", epsilon);
  exp.tag("delta", delta);
  exp.tag("rate", rate);
  exp.tag("p_low", p_low);
  exp.tag("seed", seed);
  exp.tag("k", k);
  exp.tag("git-revision", std::string(g_GIT_SHA1));
  exp.tag("theory-samples-fraction", theory_samples_fraction);
  exp.tag("num-threads", omp_threads);
  
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_INFO("Loaded graph with " << boost::num_vertices(graph) <<
           " nodes and " << boost::num_edges(graph) << " edges");

  check_num_components(graph, k);
  
  auto prob_to_samples = [epsilon, delta, theory_samples_fraction](double p) {
    return theory_samples_fraction/(epsilon*epsilon*p) * log(1/delta);
  };

  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);
  
  auto start = std::chrono::steady_clock::now();

  auto clustering = average_probability_cluster(graph, sampler, k, rate, p_low, exp);

  auto end = std::chrono::steady_clock::now();
  double elapsed = std::chrono::duration_cast< std::chrono::milliseconds >(end - start).count();

  exp.append("performance", {{"time", elapsed},});
  
  add_clustering_info(graph, clustering, "clustering", exp);
  add_scores(graph, clustering, sampler, args.count("fast-scores"), exp);
  exp.save();
  LOG_INFO(elapsed << " ms elapsed.");
}
