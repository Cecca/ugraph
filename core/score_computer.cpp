#include "scores.hpp"
#include "io.hpp"
#include "json.hpp"

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
    ("epsilon", po::value<double>()->default_value(0.1),
     "tolerated absolute error")
    ("delta", po::value<double>()->default_value(0.01),
     "error probability")
    ("graph", po::value<std::string>(),
     "file containing graph data")
    ("clustering", po::value<std::string>(),
     "json file containing the clustering");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    abort();
  }
  return vm;
}



std::vector< ClusterVertex > load_clustering(const ugraph_t & graph,
                                             const std::string & path,
                                             CCSampler & sampler) {
  using json = nlohmann::json;
  std::ifstream input(path);
  json data;
  input >> data;
  input.close();
  auto clustering_table = data["tables"]["clustering"];
  std::unordered_map<ugraph_vertex_t, std::vector<ugraph_vertex_t> > clusters_map;
  
  for (const auto & row : clustering_table) {
    ugraph_vertex_t
      center = _find_id(graph, row["center label"]),
      vertex = _find_id(graph, row["label"]);
    if (clusters_map.count(center) == 0) {
      clusters_map[center] = std::vector<ugraph_vertex_t>();
    }
    clusters_map[center].push_back(vertex);
  }

  return build_cluster_vertices(graph, clusters_map, sampler);
}

bool should_run(boost::program_options::variables_map args) {
  if (args.count("overwrite") > 0) {
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

  std::string graph_path = args["graph"].as<std::string>();
  ugraph_t graph;
  read_edge_list(graph, graph_path);

  double
    epsilon = args["epsilon"].as<double>(),
    delta = args["delta"].as<double>();
  auto prob_to_samples = [epsilon, delta](double p) {return 1/(epsilon*epsilon*p) * log(1/delta);};
  CCSampler sampler(graph, prob_to_samples, seed, omp_threads);
  sampler.min_probability(graph, 0.01);

  std::string clustering_path = args["clustering"].as<std::string>();
  auto vinfo = load_clustering(graph, clustering_path, sampler);
  auto clusters = build_clusters(vinfo);

  LOG_INFO("Computing minimum probability");
  probability_t min_p = min_probability(vinfo);
  LOG_INFO("Computing ACR");
  double acr = average_cluster_reliability(graph, clusters, sampler);
  LOG_INFO("Computing AVPR");
  double avpr = average_vertex_pairwise_reliability(graph, clusters, sampler);
  
  LOG_INFO("Clustering with:" <<
           "\n\tp_min = " << min_p <<
           "\n\tavpr  = " << avpr <<
           "\n\tacr   = " << acr);

  using json = nlohmann::json;
  // add the scores to the json file
  std::ifstream input(clustering_path);
  json data;
  input >> data;
  input.close();
  data["tables"]["scores"][0]["p_min"] = min_p;
  data["tables"]["scores"][0]["acr"] = acr;
  data["tables"]["scores"][0]["avpr"] = avpr;
  
  LOG_INFO("Writing the result in place in file " << clustering_path);
  std::ofstream output(clustering_path);
  output << data << std::endl;
  output.close();
  return 0;
}

