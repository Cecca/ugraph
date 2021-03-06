#pragma once

#include "cc_sampler.hpp"
#include "cluster_vertex.hpp"
#include "experiment_reporter.hpp"

std::vector< ClusterVertex > build_cluster_vertices(const ugraph_t & graph,
                                                    const std::unordered_map< ugraph_vertex_t, std::vector<ugraph_vertex_t> > & clusters,
                                                    CCSampler & sampler);

probability_t min_probability(const std::vector< ClusterVertex > & vinfo);

probability_t sum_probability(const std::vector< ClusterVertex > & vinfo);

size_t num_centers(const std::vector< ClusterVertex > & vinfo);

double average_cluster_reliability(const ugraph_t & graph,
                                   const std::vector< std::vector<ugraph_vertex_t> > & clusters,
                                   CCSampler & sampler);

// Deprecated
double average_vertex_pairwise_reliability_old(const ugraph_t & graph,
                                           std::vector<std::vector<ugraph_vertex_t> > & clusters,
                                           CCSampler & sampler);

struct AVPR {
  double inner;
  double outer;

  AVPR(): inner(-1), outer(-1) {}
  AVPR(double inner, double outer): inner(inner), outer(outer) {
    REQUIRE(inner <= 1.0, "Inner AVPR should be less than 1.0");
    REQUIRE(outer <= 1.0, "Outer AVPR should be less than 1.0");
  }
};

/// Computes the Average Vertex Pairwise Reliability
AVPR average_vertex_pairwise_reliability(const ugraph_t & graph,
            const std::vector<ClusterVertex> & vinfo,
            CCSampler & sampler);

std::vector< std::vector< ugraph_vertex_t > >
build_clusters(const std::vector< ClusterVertex > & vinfo);

void add_scores(const ugraph_t & graph,
                const std::vector< ClusterVertex > & vinfo,
                CCSampler & sampler,
                const bool with_acr,
                const bool with_avpr);
