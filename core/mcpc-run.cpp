#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "logging.hpp"
#include "git_info.hpp"
#include "mcpc.hpp"
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
    ("slack,s", po::value<size_t>()->default_value(0),
     "Allow a number `s` of nodes to be singletons clusters")
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
    ("with-acr", "also compute ACR measure")
    ("with-avpr", "also compute AVPR");

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

void check_num_components(const ugraph_t & graph, const size_t target, const size_t slack) {
  auto component_map = boost::make_vector_property_map<int>(boost::get(boost::vertex_index, graph));
  size_t num_components = boost::connected_components(graph, component_map);
  if (target + slack < num_components) {
    LOG_ERROR("The target size ("
              << target
              << " + "
              << slack
              << ") is smaller than the number of connected components ("
              << num_components << "): the algorithm can't terminate");
    throw std::logic_error("Target size too small");
  }
}

void add_clustering_info(const ugraph_t &graph,
                         const std::vector<ClusterVertex> &vinfo,
                         const std::string & table_name) {
  size_t n = vinfo.size();
  for (ugraph_vertex_t v = 0; v < n; v++) {
    ugraph_vertex_t center = vinfo[v].center();
    EXPERIMENT_APPEND(table_name.c_str(), {{"id", v},
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

  size_t
    k = args["target"].as<size_t>(),
    slack = args["slack"].as<size_t>();

  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");  

  
  EXPERIMENT_TAG("algorithm", std::string("k-center"));
  EXPERIMENT_TAG("input", graph_path);
  EXPERIMENT_TAG("epsilon", epsilon);
  EXPERIMENT_TAG("delta", delta);
  EXPERIMENT_TAG("rate", rate);
  EXPERIMENT_TAG("p_low", p_low);
  EXPERIMENT_TAG("seed", (size_t) seed);
  EXPERIMENT_TAG("k", k);
  EXPERIMENT_TAG("slack", slack);
  EXPERIMENT_TAG("git-revision", std::string(g_GIT_SHA1));
  EXPERIMENT_TAG("theory-samples-fraction", theory_samples_fraction);
  EXPERIMENT_TAG("num-threads", omp_threads);
  
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_INFO("Loaded graph with " << boost::num_vertices(graph) <<
           " nodes and " << boost::num_edges(graph) << " edges");

  check_num_components(graph, k, slack);
  
  auto prob_to_samples = [epsilon, delta, theory_samples_fraction](double p) {
    return theory_samples_fraction/(epsilon*epsilon*p) * log(1/delta);
  };

  Splitmix64 seeder(seed);
  Xorshift1024star rnd(seeder.next());
  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);
  
  std::vector<ClusterVertex> clustering;
  
  auto start = std::chrono::steady_clock::now();
  
  if (args.count("depth") > 0) {
    size_t depth = args["depth"].as<size_t>();
    EXPERIMENT_TAG("depth", depth);
    // Override the sampler, using the limited depth one
    BfsSampler sampler(graph, depth, prob_to_samples, seed, omp_threads);
    clustering = minimum_connection_probability_clustering(graph, sampler, k, slack, rate, p_low, rnd);
  } else {
    EXPERIMENT_TAG("depth", std::numeric_limits<double>::infinity());
    clustering = minimum_connection_probability_clustering(graph, sampler, k, slack, rate, p_low, rnd);
  }
  
  auto end = std::chrono::steady_clock::now();
  double elapsed = std::chrono::duration_cast< std::chrono::milliseconds >(end - start).count();

  EXPERIMENT_APPEND("performance", {{"time", elapsed},});
  
  add_clustering_info(graph, clustering, "clustering");
  add_scores(graph, clustering, sampler, args.count("with-acr") > 0, args.count("with-avpr") > 0);
  EXPERIMENT_SAVE();
  LOG_INFO(elapsed << " ms elapsed.");
}
