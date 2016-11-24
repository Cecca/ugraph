#pragma once

#include "sampler.hpp"
#include "cluster_vertex.hpp"
#include "experiment_reporter.hpp"

std::vector< ClusterVertex > build_cluster_vertices(const ugraph_t & graph,
                                                    const std::unordered_map< ugraph_vertex_t, std::vector<ugraph_vertex_t> > & clusters,
                                                    CCSampler & sampler);
probability_t min_probability(const std::vector< ClusterVertex > & vinfo);

double average_cluster_reliability(const ugraph_t & graph,
                                   const std::vector< std::vector<ugraph_vertex_t> > & clusters,
                                   CCSampler & sampler);

double average_vertex_pairwise_reliability(const ugraph_t & graph,
                                           std::vector<std::vector<ugraph_vertex_t> > & clusters,
                                           CCSampler & sampler);

std::vector< std::vector< ugraph_vertex_t > >
build_clusters(const std::vector< ClusterVertex > & vinfo);

void add_scores(const ugraph_t & graph,
                const std::vector< ClusterVertex > & vinfo,
                CCSampler & sampler,
                ExperimentReporter & experiment);
