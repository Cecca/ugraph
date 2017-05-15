#include "cc_sampler.hpp"
#include "experiment_reporter.hpp"
#include "git_info.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "rand.hpp"
#include "scores.hpp"
#include "types.hpp"

#include <boost/graph/dijkstra_shortest_paths.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
    gmm_graph_t;

struct GmmVertexInfo {

  GmmVertexInfo()
      : center(0), distance(std::numeric_limits<double>::infinity()) {}

  GmmVertexInfo(ugraph_vertex_t c, double d) : center(c), distance(d) {}

  ugraph_vertex_t center;
  double distance;
};

void copy_graph(const ugraph_t &g_in, gmm_graph_t &g_out) {
  using namespace boost;
  size_t n = num_vertices(g_in);
  for (ugraph_vertex_t i; i < n; i++) {
    add_vertex(g_out);
  }

  BGL_FORALL_EDGES(e, g_in, ugraph_t) {
    double dist = std::log(1.0 / g_in[e].probability);
    add_edge(source(e, g_in), target(e, g_in), dist, g_out);
  }
}

void gmm(const gmm_graph_t &g, size_t k, Xorshift1024star &rnd,
         std::vector<GmmVertexInfo> &vertices) {
  size_t n = boost::num_vertices(g);
  vertices.resize(n);
  std::fill(vertices.begin(), vertices.end(), GmmVertexInfo());

  boost::vector_property_map<double> distances(n);

  ugraph_vertex_t center = std::floor(rnd.next_double() * n);
  size_t num_centers = 1;
  LOG_DEBUG("New center: " << center);
  boost::dijkstra_shortest_paths(g, center, boost::distance_map(distances));
  for (ugraph_vertex_t i = 0; i < n; i++) {
    if (vertices[i].distance > get(distances, i)) {
      vertices[i].distance = get(distances, i);
      vertices[i].center = center;
    }
  }

  while (num_centers < k) {
    size_t max_d = 0.0;
    size_t center = 0;
    for (ugraph_vertex_t i = 0; i < n; i++) {
      if (vertices[i].distance > max_d) {
        max_d = vertices[i].distance;
        center = i;
      }
    }
    num_centers++;
    LOG_DEBUG("New center: " << center);
    boost::dijkstra_shortest_paths(g, center, boost::distance_map(distances));
    for (ugraph_vertex_t i = 0; i < n; i++) {
      if (vertices[i].distance > get(distances, i)) {
        vertices[i].distance = get(distances, i);
        vertices[i].center = center;
      }
    }
  }
}

void add_clustering_info(const ugraph_t &graph,
                         const std::vector<ClusterVertex> &vinfo,
                         const std::string &table_name) {
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

boost::program_options::variables_map parse_args(int argc, char **argv) {
  using namespace boost;
  namespace po = boost::program_options;

  po::options_description desc("Simple clustering");
  desc.add_options()
    ("help", "produce help message")
    ("revision", "get the git revision of the code")
    ("debug", "print debug output")
    ("trace", "print trace output")
    ("graph", po::value<std::string>(), "input graph")
    ("target,k", po::value<size_t>(), "desired number of clusters")
    ("epsilon", po::value<double>()->default_value(0.1),
     "tolerated absolute error")
    ("delta", po::value<double>()->default_value(0.01), "error probability")
    ("theory-samples-fraction", po::value<double>()->default_value(0.1),
     "Fraction of samples to be used with respect to the theory-defined "
     "formula")
    ("seed", po::value<uint64_t>(), "seed for random generator")
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

int main(int argc, char **argv) {
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

  double epsilon = args["epsilon"].as<double>(),
         delta = args["delta"].as<double>(),
         theory_samples_fraction = args["theory-samples-fraction"].as<double>(),
         p_low = 0.1;

  size_t k = args["target"].as<size_t>();

  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");

  
  EXPERIMENT_TAG("algorithm", std::string("gmm"));
  EXPERIMENT_TAG("input", graph_path);
  EXPERIMENT_TAG("p_low", p_low);
  EXPERIMENT_TAG("epsilon", epsilon);
  EXPERIMENT_TAG("delta", delta);
  EXPERIMENT_TAG("seed", seed);
  EXPERIMENT_TAG("k", k);
  EXPERIMENT_TAG("git-revision", std::string(g_GIT_SHA1));
  EXPERIMENT_TAG("theory-samples-fraction", theory_samples_fraction);
  EXPERIMENT_TAG("num-threads", omp_threads);

  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_INFO("Loaded graph with " << boost::num_vertices(graph) << " nodes and "
                                << boost::num_edges(graph) << " edges");

  gmm_graph_t gmm_graph;
  copy_graph(graph, gmm_graph);

  auto prob_to_samples = [epsilon, delta, theory_samples_fraction](double p) {
    return theory_samples_fraction / (epsilon * epsilon * p) * log(1 / delta);
  };
  Splitmix64 seeder(seed);
  Xorshift1024star rnd(seeder.next());
  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);

  std::vector<GmmVertexInfo> gmm_clustering;
  auto start = std::chrono::steady_clock::now();
  gmm(gmm_graph, k, rnd, gmm_clustering);
  auto end = std::chrono::steady_clock::now();

  double elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  EXPERIMENT_APPEND("performance", {{"time", elapsed},});

  size_t n = boost::num_vertices(graph);
  std::vector<probability_t> probabilities(n);
  sampler.min_probability(graph, p_low);
  // Now compute the probabilities
  std::vector<ClusterVertex> clustering(n);
  for (ugraph_vertex_t i = 0; i < n; i++) {
    if (gmm_clustering[i].center == i) {
      clustering[i].force_make_center(i);
      sampler.connection_probabilities(graph, i, probabilities);
      for (ugraph_vertex_t j = 0; j < n; j++) {
        auto & v = clustering[j];
        if (gmm_clustering[j].center == i &&
            (!v.is_covered() || v.probability() <= probabilities[j])) {
          v.cover(i, probabilities[j]);
        }
      }
    }
  }

  add_clustering_info(graph, clustering, "clustering");
  add_scores(graph, clustering, sampler, args.count("fast-scores"));
  EXPERIMENT_SAVE();
  LOG_INFO(elapsed << " ms elapsed.");
}
