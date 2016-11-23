#pragma once

#include "types.hpp"
#include "require.hpp"
#include "sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"

std::vector< ClusterVertex > concurrent_cluster(const ugraph_t & graph,
                                                CCSampler & sampler,
                                                const size_t batch,
                                                const probability_t p_low,
                                                Xorshift1024star & rnd,
                                                ExperimentReporter & experiment);
