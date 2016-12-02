#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"

std::vector< ClusterVertex > sequential_cluster(const ugraph_t & graph,
                                                CCSampler & sampler,
                                                const size_t k,
                                                const size_t slack,
                                                const double rate,
                                                const probability_t p_low,
                                                ExperimentReporter & experiment);
