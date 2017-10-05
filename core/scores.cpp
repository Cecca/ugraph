#include "scores.hpp"

std::vector< ClusterVertex > build_cluster_vertices(const ugraph_t & graph,
                                                    const std::unordered_map< ugraph_vertex_t, std::vector<ugraph_vertex_t> > & clusters,
                                                    CCSampler & sampler) {
  std::vector< probability_t > probabilities(boost::num_vertices(graph), 0.0);
  std::vector< ClusterVertex > vinfo(boost::num_vertices(graph));

  LOG_INFO("Building ClusterVertex");

  const size_t n_clusters = clusters.size();
  const size_t n_clusters_one_cent = n_clusters / 100;
  size_t cluster_idx=0;

  LOG_DEBUG("Starting loop over clusters");
  for (const auto & entry : clusters) {
    if (n_clusters_one_cent != 0 && cluster_idx % n_clusters_one_cent == 0) {
      LOG_INFO("  " << (cluster_idx / n_clusters_one_cent) << "%");
    }
    cluster_idx++;
    const ugraph_vertex_t center = entry.first;
    auto cluster = entry.second;
    LOG_DEBUG("Looking at probabilities");
    sampler.connection_probabilities(graph, center, cluster, probabilities);
    LOG_DEBUG("Looked at probabilities");
    for(const ugraph_vertex_t v : cluster) {
      if (v == center) {
        LOG_DEBUG("Turning " << v << " into a center");
        vinfo[v].make_center(v);
      } else {
        LOG_DEBUG("Cover node " << v << " from " << center << " with p=" << probabilities[v]);
        vinfo[v].cover(center, probabilities[v]);
      }
    }
  }

  for (size_t i =0; i<vinfo.size(); i++) {
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

probability_t sum_probability(const std::vector< ClusterVertex > & vinfo) {
  probability_t sum = 0.0;
  for (const auto & v : vinfo) {
    probability_t p = v.probability();
    sum += p;
  }
  return sum;
}

size_t num_centers(const std::vector< ClusterVertex > & vinfo) {
  size_t sum = 0;
  for (const auto & v : vinfo) {
    if (v.is_center()) {
      sum++;
    }
  }
  return sum;
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

AVPR average_vertex_pairwise_reliability(const ugraph_t & graph,
            const std::vector<ClusterVertex> & vinfo,
            CCSampler & sampler) {
  // Instead of looking at the connection probabilities of single
  // pairs, we consider the connected components of each sample and,
  // for each cluster, we add (x \choose 2) to the counter of that
  // cluster, where x is the size of the intersection of the cluster
  // with a connected component of the sample.

  const size_t n = vinfo.size();
  
  // First, map cluster centers to contiguous identifiers
  std::unordered_map<ugraph_vertex_t, size_t> cluster_ids;
  size_t cur_id = 0;
  for (const auto & v : vinfo) {
    if (v.is_center()) {
      cluster_ids[v.center()] = cur_id++;
    }
  }
  const size_t n_clusters = cluster_ids.size();
  std::vector<size_t> cluster_sizes(n_clusters, 0);
  for(const auto & v: vinfo) {
    ugraph_vertex_t center = v.center();
    cluster_sizes[cluster_ids[center]]++;
  }
  
  const size_t n_threads = omp_get_max_threads();
  
  // The vectors that will hold the counts for the clusters, one for
  // each thread
  std::vector<std::vector<size_t>> t_cluster_inner_counts(
      n_threads, std::vector<size_t>(n_clusters, 0));
  std::vector<std::vector<size_t>> t_cluster_outer_counts(
      n_threads, std::vector<size_t>(n_clusters, 0));
  // The samples
  const std::vector<CCSampler::component_vector_t> &samples = sampler.get_samples();
  const size_t n_samples = samples.size();
  size_t processed_samples = 0;
  
  // For each sample in parallel, accumulate counts
#pragma omp parallel for
  for(size_t sample_idx=0; sample_idx < n_samples; sample_idx++){
    const auto & sample = samples[sample_idx];
    REQUIRE(sample.size() == n, "Samples are of the wrong size!");
    std::unordered_map<size_t, size_t> connected_components_ids;
    size_t cur_comp_id = 0;
    for(const auto cc_id : sample) {
      // Set the ID if not already done
      if (connected_components_ids.count(cc_id) == 0) {
        connected_components_ids[cc_id] = cur_comp_id++;
      }
    }
    const size_t num_connected_components = connected_components_ids.size();
    std::vector<size_t> connected_components_sizes(num_connected_components, 0);
    for(const auto cc_id : sample)  {
      connected_components_sizes[connected_components_ids[cc_id]]++;
    }

    const auto tid = omp_get_thread_num();
    auto& cluster_inner_counts = t_cluster_inner_counts[tid];
    auto& cluster_outer_counts = t_cluster_outer_counts[tid];

    // A (sparse) matrix of `num_clusters` x `num_components` elements that
    // contains in element (i,j) the number of elements of cluster i
    // belonging to the connected component j.
    std::vector<std::unordered_map<size_t, size_t>> intersection_sizes(n_clusters);
    // for (size_t i=0; i<n_clusters; i++) {
    //   std::vector<size_t> vec;
    //   for(size_t j=0; j<num_connected_components; j++) {
    //     vec.push_back(0);
    //   }
    //   intersection_sizes.push_back(vec);
    // }
    // Populate the intersection counts
    for (size_t i=0; i<n; i++) {
      const size_t cluster_id = cluster_ids[vinfo[i].center()];
      const size_t component_id = connected_components_ids[sample[i]];
      if (intersection_sizes.at(cluster_id).count(component_id) == 0) {
        intersection_sizes[cluster_id][component_id] = 1;
      } else {
        intersection_sizes[cluster_id][component_id] += 1;
      }
    }
    
    for (size_t cluster_idx=0; cluster_idx<n_clusters; cluster_idx++) {
      size_t inner_cnt=0;
      size_t outer_cnt=0;
      auto & i_sizes = intersection_sizes[cluster_idx];
      // TODO: Apply SIMD reduction
      for (size_t component_idx=0; component_idx<num_connected_components; component_idx++){
        size_t intersection = 0;
        if (i_sizes.count(component_idx) > 0) {
          intersection = i_sizes[component_idx];
        }
        size_t difference = connected_components_sizes[component_idx] - intersection;
        inner_cnt += intersection*(intersection-1)/2;
        size_t outer_pairs = intersection*difference;
        size_t max_pairs = cluster_sizes[cluster_idx]*connected_components_sizes[component_idx];
        outer_cnt += outer_pairs;
      }
      cluster_inner_counts[cluster_idx] += inner_cnt;
      cluster_outer_counts[cluster_idx] += outer_cnt;
    }
  }

  double inner_denominator = 0.0;
  for(const auto & size : cluster_sizes) {
    inner_denominator += (size*(size-1))/2.0;
  }
  double outer_denominator = 0.0;
  for(const auto in_cluster_size : cluster_sizes) {
    size_t out_cluster_size = n - in_cluster_size;
    outer_denominator += in_cluster_size*out_cluster_size;
  }
  
  double inner_numerator = 0.0;
  double outer_numerator = 0.0;
  for (const auto & cnts : t_cluster_inner_counts) {
    for (const size_t cnt : cnts) {
      inner_numerator += (((double) cnt) / n_samples);
    }
  }
  for (const auto & cnts : t_cluster_outer_counts) {
    for (const size_t cnt : cnts) {
      outer_numerator += (((double) cnt) / n_samples);
    }
  }
  AVPR avpr;
  avpr.inner = inner_numerator / inner_denominator;
  avpr.outer = outer_numerator / outer_denominator;
  REQUIRE(avpr.inner <= 1.0, "Inner AVPR must be smaller than 1");
  REQUIRE(avpr.outer <= 1.0, "Outer AVPR must be smaller than 1");
  return avpr;
}

double average_vertex_pairwise_reliability_old(const ugraph_t & graph,
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
                const bool with_acr,
                const bool with_avpr) {
  std::vector< std::vector< ugraph_vertex_t > > clusters = build_clusters(vinfo);

  size_t num_clusters = 0;
  for(const auto & v : vinfo) {
    if (v.is_center()) num_clusters++;
  }
  
  LOG_INFO("Computing minimum probability");
  probability_t min_p = min_probability(vinfo);
  probability_t sum_p = sum_probability(vinfo);
  probability_t avg_p = sum_p / boost::num_vertices(graph);
  LOG_INFO("Sum_p " << sum_p << " avg_p " << avg_p);
    
  double acr=-1;
  AVPR avpr;
  if (with_acr) {
    sampler.min_probability(graph, min_p);
    LOG_INFO("Computing ACR");
    acr = average_cluster_reliability(graph, clusters, sampler);
  }
  if (with_avpr) {
    sampler.min_probability(graph, min_p);
    LOG_INFO("Computing AVPR");
    avpr = average_vertex_pairwise_reliability(graph, vinfo, sampler);
  }
    
  LOG_INFO("Clustering with:" <<
           "\n\t# clusters = " << num_clusters << 
           "\n\tp_min = " << min_p <<
           "\n\taverage p = " << avg_p <<
           "\n\tinner-avpr      = " << avpr.inner <<
           "\n\touter-avpr      = " << avpr.outer <<
           "\n\tacr   = " << acr);
  EXPERIMENT_APPEND("scores", {{"acr", acr},
                               {"p_min", min_p},
                               {"inner-avpr", avpr.inner},
                               {"outer-avpr", avpr.outer},
                               {"probabilities sum", sum_p},
                               {"average probability", avg_p},
                               {"num clusters", num_clusters}});
}
