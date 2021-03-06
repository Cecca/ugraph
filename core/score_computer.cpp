#include "scores.hpp"
#include "io.hpp"
#include "json.hpp"
#include <boost/filesystem.hpp>

ugraph_vertex_t _find_id(const ugraph_t & graph, std::string label) {
  using namespace boost;
  BGL_FORALL_VERTICES(v, graph, ugraph_t) {
    if (graph[v].label == label) {
      return v;
    }
  }
  throw std::logic_error("Label not found");
}


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
    ("with-acr", "also compute ACR measure")
    ("with-avpr", "also compute AVPR")
    ("graph", po::value<std::string>(),
     "file containing graph data")
    ("odir", po::value<std::string>(), "output directory")
    ("clustering", po::value<std::string>(),
     "json file containing the clustering");
  
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



std::vector< ClusterVertex > load_clustering(const ugraph_t & graph,
                                             const std::string & path,
                                             CCSampler & sampler) {
  std::unordered_map<std::string, ugraph_vertex_t> id_map;
  using namespace boost;
  BGL_FORALL_VERTICES(v, graph, ugraph_t) {
    id_map[graph[v].label] = v;
  }
  
  LOG_INFO("Loading json data");
  using json = nlohmann::json;
  std::ifstream input(path);
  json data;
  input >> data;
  input.close();
  LOG_INFO("Loaded json data");
  auto clustering_table = data["tables"]["clustering"];
  std::unordered_map<ugraph_vertex_t, std::vector<ugraph_vertex_t> > clusters_map;
  LOG_INFO("Building cluster map");
  const size_t num_nodes = clustering_table.size();
  const size_t one_cent = num_nodes / 100;
  size_t i=0;
  for (const auto & row : clustering_table) {
    if (i % (10*one_cent) == 0) {
      LOG_INFO("Progress: " << (i / one_cent) << "% (" << i << "/" << num_nodes << ")"); 
    }
    i++;
    ugraph_vertex_t
      center = id_map[row["center label"]],// _find_id(graph, row["center label"]),
      vertex = id_map[row["label"]];// _find_id(graph, row["label"]);
    LOG_DEBUG("Loading vertex " << vertex << " assigned to " << center);
    if (clusters_map.count(center) == 0) {
      clusters_map[center] = std::vector<ugraph_vertex_t>();
      clusters_map.reserve(1024);
    }
    clusters_map[center].push_back(vertex);
  }
  LOG_INFO("Built cluster map");

  return build_cluster_vertices(graph, clusters_map, sampler);
}

bool should_run(boost::program_options::variables_map args) {
  if (args.count("overwrite") > 0 || args.count("odir") > 0) {
    return true;
  }
  using json = nlohmann::json;
  std::string path(args["clustering"].as<std::string>());
  
  std::ifstream input(path);
  json data;
  input >> data;
  return !(data["tables"]["scores"][0].count("p_min") > 0 &&
           data["tables"]["scores"][0].count("acr") > 0 &&
           data["tables"]["scores"][0].count("avpr") > 0 );  
}


int main(int argc, char *argv[]) {
  
  auto args = parse_args(argc, argv);

  if (!should_run(args)) {
    LOG_INFO("Scores are already present, not running");
    return 0;
  }
  
  std::random_device rd;
  uint64_t seed = rd();
  auto omp_threads = omp_get_max_threads();
  LOG_INFO("Running with " << omp_threads << " threads");

  std::string clustering_path = args["clustering"].as<std::string>();

  LOG_INFO("Loading json data");
  using json = nlohmann::json;
  std::ifstream input(clustering_path);
  json data;
  input >> data;
  input.close();
  LOG_INFO("Loaded json data");
  
  std::string graph_path;
  if (args.count("graph")) {
    graph_path = args["graph"].as<std::string>();
  } else if (data["tags"].count("input")) {
    graph_path = data["tags"]["input"];
  } else if (data["tags"].count("graph")) {
    graph_path = data["tags"]["graph"];
  } else {
    LOG_ERROR("Could not retrieve input graph from clustering file");
    return -1;
  }
  LOG_DEBUG("Loading graph from " << graph_path);
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  LOG_DEBUG("Loaded graph");

  double
    epsilon = args["epsilon"].as<double>(),
    delta = args["delta"].as<double>();
  auto prob_to_samples = [epsilon, delta](double p) {return 1/(epsilon*epsilon*p) * log(1/delta);};
  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);
  sampler.min_probability(graph, 0.01);
  
  auto vinfo = load_clustering(graph, clustering_path, sampler);
  auto clusters = build_clusters(vinfo);
  LOG_DEBUG("Built clusters");

  if (data["tables"]["scores"][0].count("p_min") == 0) {
    LOG_INFO("Computing minimum probability");
    probability_t min_p = min_probability(vinfo);
    data["tables"]["scores"][0]["p_min"] = min_p;
  }
  if (data["tables"]["scores"][0].count("average probability") == 0) {
    probability_t sum_p = sum_probability(vinfo);
    data["tables"]["scores"][0]["average probability"] = sum_p / boost::num_vertices(graph);
  }
  if (data["tables"]["scores"][0].count("num clusters") == 0) {
    size_t num_clusters = num_centers(vinfo);
    data["tables"]["scores"][0]["num clusters"] = num_clusters;
  }
  
  if (data["tables"]["scores"][0].count("acr") == 0 && args.count("with-acr") > 0) {
    LOG_INFO("Computing ACR");
    double acr = average_cluster_reliability(graph, clusters, sampler);
    data["tables"]["scores"][0]["acr"] = acr;
  }
  if (data["tables"]["scores"][0].count("inner-avpr") == 0 && args.count("with-avpr") > 0) {
    LOG_INFO("Computing AVPR");
    AVPR avpr = average_vertex_pairwise_reliability(graph, vinfo, sampler);

    data["tables"]["scores"][0]["inner-avpr"] = avpr.inner;
    data["tables"]["scores"][0]["outer-avpr"] = avpr.outer;
  }

  if (args.count("odir")) {
    boost::filesystem::path dir(args["odir"].as<std::string>());
    if (!boost::filesystem::is_directory(dir)) {
      LOG_INFO("Creating directory " << dir);
      boost::filesystem::create_directory(dir);
    }
    boost::filesystem::path clustering_fs_path (clustering_path);
    boost::filesystem::path out = dir / clustering_fs_path.filename();
    LOG_INFO("Writing the result in file " << out);
    std::ofstream output(out.string());
    json clustering_table = json::array();
    for (ugraph_vertex_t v = 0; v < boost::num_vertices(graph); v++) {
      ugraph_vertex_t center = vinfo[v].center();
      clustering_table.push_back({{"id", v},
                              {"center", center},
                              {"label", graph[v].label},
                              {"center label", graph[center].label},
                              {"probability", vinfo[v].probability()}});
    }
    data["tables"]["clustering"] = clustering_table;
    output << data << std::endl;
    output.close();
  } else {
    LOG_INFO("Writing the result in place in file " << clustering_path);
    std::ofstream output(clustering_path);
    output << data << std::endl;
    output.close();
  }
  return 0;
}

