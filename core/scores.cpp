#include "scores.hpp"

std::vector< ClusterVertex > build_cluster_vertices(const ugraph_t & graph,
                                                    const std::unordered_map< ugraph_vertex_t, std::vector<ugraph_vertex_t> > & clusters,
                                                    CCSampler & sampler) {
  std::vector< probability_t > probabilities(boost::num_vertices(graph), 0.0);
  std::vector< ClusterVertex > vinfo(boost::num_vertices(graph));

  LOG_INFO("Building ClusterVertex");
  
  for (const auto & entry : clusters) {
    const ugraph_vertex_t center = entry.first;
    auto cluster = entry.second;
    std::sort(cluster.begin(), cluster.end());
    LOG_DEBUG("Cluster " << center << " ------------");
    std::stringstream sstr;
    for(auto v : cluster) {
      sstr << v << " ";
    }
    LOG_DEBUG("" << sstr.str());
    sampler.connection_probabilities(graph, center, probabilities);
    for(const ugraph_vertex_t v : cluster) {
      if (v == center) {
        vinfo[v].make_center(v);
      } else {
        LOG_DEBUG("Cover node " << v << " from " << center << " with p=" << probabilities[v]);
        vinfo[v].cover(center, probabilities[v]);
      }
    }
  }

  for (size_t i =0; i<vinfo.size(); i++) {
    //LOG_INFO("Node:" << i << " prob: " << vinfo[i].probability() << " center:" << vinfo[i].center());
    if (!vinfo[i].is_covered()) {
      LOG_INFO("Node " << i << " is uncovered (" << graph[i].label << ")");
    }
  }

  LOG_INFO("Built ClusterVertex");
  
  return vinfo;
}

probability_t min_probability(const std::vector< ClusterVertex > & vinfo) {
  probability_t min_p = 1.0;
  for (const auto & v : vinfo) {
    probability_t p = v.probability();
    min_p = std::min(min_p, p);
  }
  return min_p;
}

double average_cluster_reliability(const ugraph_t & graph,
                                   const std::vector< std::vector<ugraph_vertex_t> > & clusters,
                                   CCSampler & sampler) {

  double sum = 0.0;
  // The denominator is the sum of the sizes of all the clusters. In
  // case of overlapping clusterings, it may be greater than the number
  // of nodes in the graph.
  size_t denominator = 0;
  for(const auto & cluster : clusters) {
    probability_t reliability = sampler.connection_probability(graph, cluster);
    sum += cluster.size() * reliability;
    denominator += cluster.size();
  }
  
  double acr = sum / denominator;
  REQUIRE(acr <= 1.0, "ACR greater than one!");
  return acr;
}

double average_vertex_pairwise_reliability(const ugraph_t & graph,
                                           std::vector<std::vector<ugraph_vertex_t> > & clusters,
                                           CCSampler & sampler) {
  std::vector< probability_t > probabilities(boost::num_vertices(graph), 0.0);

  double numerator = 0.0;
  double denominator = 0.0;

  for (const auto & cluster : clusters) {
    if (cluster.size() >= 2) {
      for (ugraph_vertex_t u : cluster) {
        sampler.connection_probabilities(graph, u, cluster, probabilities);
        for (ugraph_vertex_t v : cluster) {
          if (u < v) {
            numerator += probabilities[v];
          }
        }
      }
    }
    const size_t size = cluster.size();
    denominator += size*(size-1);
  }
  
  double avpr = (2*numerator) / denominator;
  REQUIRE(avpr <= 1.0, "AVPR score greater than one!");
  return avpr;
}

std::vector< std::vector< ugraph_vertex_t > >
build_clusters(const std::vector< ClusterVertex > & vinfo) {
  const size_t n = vinfo.size();
  std::map< ugraph_vertex_t, std::vector< ugraph_vertex_t > > clusters_map;
  for (ugraph_vertex_t v=0; v<n; v++) {
    clusters_map[vinfo[v].center()].push_back(v);
  }
  std::vector< std::vector< ugraph_vertex_t > > clusters;
  for (const auto & entry : clusters_map) {
    clusters.push_back(entry.second);
  }
  return clusters;
}

void add_scores(const ugraph_t & graph,
                const std::vector< ClusterVertex > & vinfo,
                CCSampler & sampler,
                const bool only_p_min,
                ExperimentReporter & experiment) {
  std::vector< std::vector< ugraph_vertex_t > > clusters = build_clusters(vinfo);
  
  LOG_INFO("Computing minimum probability");
  probability_t min_p = min_probability(vinfo);
  sampler.min_probability(graph, min_p);
  if (only_p_min) {
    LOG_INFO("Clustering with:" << "\n\tp_min = " << min_p);
    experiment.append("scores", {{"p_min", min_p},});
    return;
  }
  LOG_INFO("Computing ACR");
  double acr = average_cluster_reliability(graph, clusters, sampler);
  LOG_INFO("Computing AVPR");
  double avpr = average_vertex_pairwise_reliability(graph, clusters, sampler);

  LOG_INFO("Clustering with:" <<
           "\n\tp_min = " << min_p <<
           "\n\tavpr  = " << avpr <<
           "\n\tacr   = " << acr);
  experiment.append("scores", {{"acr", acr}, {"p_min", min_p}, {"avpr", avpr}});
}
