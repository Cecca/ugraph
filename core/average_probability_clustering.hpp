#pragma once

#include "types.hpp"
#include "require.hpp"
#include "cc_sampler.hpp"
#include "bfs_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include "guesser.hpp"
#include "counts_cache.hpp"
#include "termcolor.hpp"

class APCExponentialGuesser {

public:
  APCExponentialGuesser(const probability_t gamma, const probability_t p_low)
    : m_p_low(p_low), m_gamma(gamma), m_upper(1.0), m_lower(1.0), m_avg_p(0.0),
      binary_search(false), m_i(0) {}

  void update(probability_t avg_p) {
    if (binary_search) {
      if (avg_p > m_avg_p) {
        m_avg_p = avg_p;
        m_lower = (m_upper + m_lower) / 2;
        LOG_DEBUG("[APCExponentialGuesser] Below in binary search (new m_lower " << m_lower << ")");
      } else {
        m_upper = (m_upper + m_lower) / 2;
        LOG_DEBUG("[APCExponentialGuesser] Above in binary search (new m_upper " << m_upper << ")");
      }
    } else {
      // We decrease the guess until the score starts to decrease
      if (avg_p > m_avg_p) {
        m_avg_p = avg_p;
        m_upper = (m_i >= 2)? (1.0 - m_gamma * (1 << (m_i - 2))) : 1.0;
        m_lower = 1.0 - m_gamma * (1 << m_i);
        m_i++;
        if (m_lower <= m_p_low) {
          m_lower = m_p_low;
          binary_search = true;
        }
        LOG_DEBUG("[APCExponentialGuesser] Above in exponential search. New m_lower " << m_lower << " new m_upper " << m_upper);
      } else {
        // As soon as avg_p decreases, we start the binary search
        LOG_DEBUG("[APCExponentialGuesser] Below in exponential search: start binary search");
        binary_search = true; // Start binary search on the next call of `guess`
      }
    }
  }
  
  probability_t guess() {
    REQUIRE(m_lower <= m_upper, "Upper and lower bounds inverted!");
    if (binary_search) {
      // The last guess is the lower bound, so to have always a
      // probability less than the necessary one
      // if (stop()) { return m_lower; }
      return (m_upper + m_lower) / 2;
    } else {
      return m_lower;
    }
  }

  bool stop() const {
    if (binary_search) {
      LOG_DEBUG("[APCExponentialGuesser] Stopping condition: " << (1.0-m_lower/m_upper) << "<=" << m_gamma);
    }
    return binary_search && (1.0-m_lower/m_upper) <= m_gamma;
  }

private:
  const probability_t m_p_low;
  probability_t m_gamma;
  probability_t m_upper;
  probability_t m_lower;
  probability_t m_avg_p;

  bool binary_search;
  size_t m_i;
};


size_t count_uncovered(const std::vector< ClusterVertex > & vinfo) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered()) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            const std::vector< ClusterVertex > & vinfo) {
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered() || vinfo[i].probability() < p_curr) {
      return i;
    }
  }
  throw std::logic_error("No uncovered node to select");
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            ConnectionCountsCache & cccache,
                            const std::vector< ClusterVertex > & vinfo) {
  size_t n = vinfo.size();
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (cccache.contains(i)) {
      if (!vinfo[i].is_covered() || vinfo[i].probability() < p_curr) {
        return i;
      } else if (vinfo[i].is_covered() && !vinfo[i].is_center()) {
        // reset counter to mark for eviction for nodes that are
        // covered but are not centers.
        cccache.set_accessed(i, 0);
      }
    }
  }
  
  return pick_vertex(graph, p_curr, vinfo);
}

double sum_center_connection_probabilities(const std::vector< ClusterVertex > & vinfo) {
  double sum = 0.0;
  for (const ClusterVertex & v : vinfo) {
    if (v.is_covered()) {
      sum += v.probability();
    }
  }
  return sum;
}

template<typename Sampler>
std::vector< ClusterVertex >
average_probability_cluster(const ugraph_t & graph,
                            Sampler & sampler,
                            const size_t k,
                            const double rate,
                            const probability_t p_low,
                            ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > probabilities(n);
  ConnectionCountsCache cccache(std::min(k, 500ul));
  size_t iteration = 0;
  probability_t p_curr = 1.0;
  probability_t reliable_estimate_lower_bound = 1.0;
  APCExponentialGuesser guesser(rate, p_low);
  size_t uncovered = n;
  double max_sum = 0.0;

  while (!guesser.stop()) {
    LOG_INFO(">>> Build clustering with p_curr=" << p_curr);
    cccache.cleanup();
    LOG_DEBUG(cccache.str());
    std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());
    uncovered = n;
    sampler.min_probability(graph, p_curr);

    // Build the clustering. Start the count from 1 because the
    // stopping condition is _inside_ the cycle
    for (size_t center_cnt = 1; center_cnt < k; center_cnt++) {
      assert(uncovered == count_uncovered(vinfo));
      ugraph_vertex_t center = pick_vertex(graph, p_curr, cccache, vinfo);
      if (!vinfo[center].is_covered()) {
        uncovered--;
      }
      vinfo[center].force_make_center(center);
      sampler.connection_probabilities_cache(graph, center, cccache, probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i=0; i<n; i++) {
        if (probabilities[i] >= reliable_estimate_lower_bound) {
          if (!vinfo[i].is_covered()) {
            vinfo[i].cover(center, probabilities[i]);
            uncovered--;
          } else if (vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
          }
        }
      }

      if (uncovered == 0) {
        break;
      }

    }
    double prob_sum = sum_center_connection_probabilities(vinfo);

    if (prob_sum >= max_sum) {
      for (ugraph_vertex_t i=0; i<n; i++) {
        valid_clustering[i] = vinfo[i];
      }
      max_sum = prob_sum;
      LOG_INFO("Average connection probability " <<
               termcolor::green << prob_sum / n << termcolor::reset);
    } else {
      LOG_INFO("Average connection probability " << prob_sum / n);
    }
    guesser.update(max_sum / n);
    
    LOG_INFO("Cache hit rate: " << std::fixed << std::setprecision(2)
             << cccache.perc_hits() << "%");
    // update the probability
    p_curr = guesser.guess();
    reliable_estimate_lower_bound =
      (p_curr < reliable_estimate_lower_bound)? p_curr : reliable_estimate_lower_bound;
    iteration++;
  }

  // Fix uncovered nodes
  // get the first center, so to assign to it uncovered nodes
  ugraph_vertex_t first_center_idx=0;
  for (; first_center_idx<n || valid_clustering[first_center_idx].is_center(); first_center_idx++){}
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!valid_clustering[i].is_covered()) {
      valid_clustering[i].cover(first_center_idx, 0.0);
    }
  }

  return valid_clustering;
  
}
