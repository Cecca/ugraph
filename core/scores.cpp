#include "sampler.hpp"
#include "cluster_vertex.hpp"

probability_t min_probability(const ugraph_t & graph,
                              const std::vector< std::vector<ugraph_vertex_t> > & clusters,
                              CCSampler & sampler) {
  
}

double average_cluster_reliability(const ugraph_t & graph,
                                   const std::vector< std::vector<ugraph_vertex_t> > & clusters,
                                   CCSampler & sampler) {

  double sum = 0.0;
  // The denominator is the sum of the sizes of all the clusters. In
  // cas of overlapping clusterings, it may be greater than the number
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
        sampler.connection_probabilities(graph, u, probabilities);
        for (ugraph_vertex_t v : cluster) {
          numerator += probabilities[v];
        }
      }
    }
    const size_t size = cluster.size();
    denominator += size*(size-1);
  }
  
  double avpr = numerator / denominator;
  REQUIRE(avpr <= 1.0, "AVPR score greater than one!");
  return avpr;
}
