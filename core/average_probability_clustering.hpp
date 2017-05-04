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

class UniformGuesser {
public:
  UniformGuesser(const probability_t gamma, const probability_t p_low)
    : m_p_low(p_low), m_gamma(gamma), m_current(1.0), m_below(false) {}

  void update(probability_t avg_p) {
    m_current -= m_gamma;
  }
  
  probability_t guess() {
    return m_current;
  }

  bool stop() const {
    return m_below || m_current <= m_p_low;
  }

private:
  const probability_t m_p_low;
  const probability_t m_gamma;
  probability_t m_current;
  bool m_below;
};

class APCExponentialGuesser {

public:
  APCExponentialGuesser(const probability_t gamma, const probability_t p_low)
    : m_p_low(p_low), m_gamma(gamma), m_upper(1.0), m_lower(1.0), m_avg_p(0.0),
      binary_search(false), m_i(0) {}

  void update(probability_t avg_p) {
    LOG_DEBUG("[APCExponentialGuesser] Updating with " << avg_p);
    if (binary_search) {
      if (avg_p < m_avg_p) {
        m_upper = (m_upper + m_lower) / 2;
        LOG_DEBUG("[APCExponentialGuesser] Going down in binary search (new m_lower " <<
                 m_lower << " m_avg_p " << m_avg_p << ")");
      } else {
        m_avg_p = avg_p;
        m_lower = (m_upper + m_lower) / 2;
        LOG_DEBUG("[APCExponentialGuesser] Going up in binary search (new m_upper " <<
                 m_upper << " m_avg_p " << m_avg_p << ")");
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
        LOG_DEBUG("[APCExponentialGuesser] Going down in exponential search. New m_lower " <<
                 m_lower << " new m_upper " << m_upper  << " m_avg_p " << m_avg_p);
      } else {
        // As soon as avg_p decreases, we start the binary search
        LOG_DEBUG("[APCExponentialGuesser] Found bottom in exponential search: start binary search" <<
                 " m_avg_p " << m_avg_p);
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

class DirectionalGuesser {

public:
  DirectionalGuesser(const probability_t gamma, const probability_t p_low)
    : m_p_low(p_low), m_gamma(gamma), m_upper(1.0), m_lower(1.0),
      m_max_avg_p(0.0), m_last_avg_p(0.0), m_second_to_last_avg_p(0.0),
      binary_search(false), direction_down(true), m_i(0) {}

  void update(probability_t avg_p) {
    if (binary_search) {
      if (avg_p > m_max_avg_p) {
        // continue in the same direction
        m_max_avg_p = avg_p;
        LOG_INFO("[Directional Guesser] continue same direction");
      } else if (avg_p > m_last_avg_p) {
        LOG_INFO("[Directional Guesser] continue same direction (without updating the max)");
      } else {
        // reverse the direction
        direction_down = !direction_down;
        LOG_INFO("[Directional Guesser] reverse direction");
      }
      if (direction_down) {
        m_upper = (m_upper + m_lower) / 2;
      } else {
        m_lower = (m_upper + m_lower) / 2;
      }
    } else {
      if (avg_p >= m_max_avg_p) {
        m_max_avg_p = avg_p;
        m_upper = (m_i >= 2)? (1.0 - m_gamma * (1 << (m_i - 2))) : 1.0;
        m_lower = 1.0 - m_gamma * (1 << m_i);
        m_i++;
        if (m_lower <= m_p_low) {
          m_lower = m_p_low;
          binary_search = true;
          direction_down = true;
          LOG_WARN("[Directional Guesser] lower bound below p_low, start search upward");
        }
      } else {
        binary_search = true;
        direction_down = false;
        LOG_INFO("[Directional Guesser] start search upward");
      }
    }
    m_second_to_last_avg_p = m_last_avg_p;
    m_last_avg_p = avg_p;
  }
  
  probability_t guess() {
    REQUIRE(m_lower <= m_upper, "Upper and lower bounds inverted!");
    if (binary_search) {
      return (m_upper + m_lower) / 2;
    } else {
      return m_lower;
    }
  }
  
  bool stop() const {
    if (binary_search) {
      LOG_DEBUG("[APCExponentialGuesser] Stopping condition: " << (1.0-m_lower/m_upper) << "<=" << m_gamma);
    }
    return (binary_search && (1.0-m_lower/m_upper) <= m_gamma) || (binary_search && m_last_avg_p == m_second_to_last_avg_p);
    //return (binary_search && (1.0-m_lower/m_upper) <= m_gamma);
  }
  
private:
  const probability_t m_p_low;
  probability_t m_gamma;
  probability_t m_upper;
  probability_t m_lower;
  probability_t m_max_avg_p;
  probability_t m_last_avg_p;
  probability_t m_second_to_last_avg_p;

  bool binary_search;
  bool direction_down;
  size_t m_i;
  
};


size_t count_uncovered(const std::vector< ClusterVertex > & vinfo, double p_curr) {
  size_t cnt = 0;
  for (const auto & v : vinfo) {
    if (!v.is_covered() || v.probability() < p_curr) {
      cnt++;
    }
  }
  return cnt;
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            const std::vector< ClusterVertex > & vinfo) {
  auto n = boost::num_vertices(graph);
  // Give precedence to really uncovered nodes
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered()) {
      return i;
    }
  }
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (vinfo[i].probability() < p_curr) {
      return i;
    }
  }
  throw std::logic_error("No uncovered node to select");
}

ugraph_vertex_t pick_vertex_rnd(const ugraph_t & graph,  
                                Xorshift1024star & rnd,
                                ConnectionCountsCache & cccache,
                                probability_t p_curr,
                                std::vector<ugraph_vertex_t> & uncovered_scratch,
                                const std::vector< ClusterVertex > & vinfo) {
  uncovered_scratch.clear();
  auto n = boost::num_vertices(graph);
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!vinfo[i].is_covered() || vinfo[i].probability() < p_curr) {
      uncovered_scratch.push_back(i);
    } else if (vinfo[i].is_covered() && cccache.contains(i)) {
      // mark the node for eviction from the cache
      cccache.set_accessed(i, 0);
    }
  }
  if (uncovered_scratch.size() == 0) {
    throw std::logic_error("No uncovered node to select");
  }
  size_t i = (size_t) std::floor(rnd.next_double()*uncovered_scratch.size());
  return uncovered_scratch[i];
}

ugraph_vertex_t pick_vertex(const ugraph_t & graph,
                            const probability_t p_curr,
                            ConnectionCountsCache & cccache,
                            const std::vector< ClusterVertex > & vinfo) {
  size_t n = vinfo.size();
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (cccache.contains(i)) {
      if (!vinfo[i].is_covered()) {
        return i;
      } else if (vinfo[i].is_covered() && !vinfo[i].is_center()) {
        // reset counter to mark for eviction for nodes that are
        // covered but are not centers.
        cccache.set_accessed(i, 0);
      }
    }
  }
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (cccache.contains(i)) {
      if (vinfo[i].probability() < p_curr) {
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
                            Xorshift1024star & rnd,
                            const size_t k,
                            const double rate,
                            const probability_t p_low,
                            ExperimentReporter & experiment) {
  const size_t n = boost::num_vertices(graph);
  std::vector< ClusterVertex > vinfo(n);
  std::vector< ClusterVertex > valid_clustering(n);
  std::vector< probability_t > probabilities(n);
  std::vector< ugraph_vertex_t > uncovered_scratch(n);
  ConnectionCountsCache cccache(std::min(k, 1000ul));
  size_t iteration = 0;
  probability_t p_curr = 1.0;
  probability_t reliable_estimate_lower_bound = 1.0;
  //DirectionalGuesser guesser(rate, p_low);
  UniformGuesser guesser(rate, p_low);
  // FIXME: Keep track of uncovered nodes
  size_t uncovered = n;
  double max_sum = 0.0;

  while (!guesser.stop()) {
    LOG_INFO(">>> Build clustering with p_curr=" << p_curr);
    cccache.cleanup();
    LOG_DEBUG(cccache.str());
    std::fill(vinfo.begin(), vinfo.end(), ClusterVertex());
    sampler.min_probability(graph, p_curr);

    // Build the clustering. Start the count from 1 because the
    // stopping condition is _inside_ the cycle
    for (size_t center_cnt = 1; center_cnt < k; center_cnt++) {
      ugraph_vertex_t center = pick_vertex_rnd(graph, rnd, cccache, p_curr, uncovered_scratch, vinfo);
      vinfo[center].force_make_center(center);
      sampler.connection_probabilities_cache(graph, center, cccache, probabilities);
      // Cover the nodes
      for (ugraph_vertex_t i=0; i<n; i++) {
        if (probabilities[i] >= reliable_estimate_lower_bound) {
          if (!vinfo[i].is_covered() || vinfo[i].probability() < probabilities[i]) {
            vinfo[i].cover(center, probabilities[i]);
          }
        }
      }

      if (count_uncovered(vinfo, p_curr) == 0) {
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
    guesser.update(prob_sum / n);
    
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
  for (; first_center_idx<n; first_center_idx++){
    if (valid_clustering[first_center_idx].is_center()) {
      break;
    }
  }
  for (ugraph_vertex_t i=0; i<n; i++) {
    if (!valid_clustering[i].is_covered()) {
      valid_clustering[i].cover(first_center_idx, 0.0);
    }
  }

  return valid_clustering;
  
}
