#pragma once

#include "types.hpp"
#include "require.hpp"
#include "sampler.hpp"
#include "experiment_reporter.hpp"

// FIXME Unify with the class SequentialClusterVertex: we don'e need
// the distinction, since there is no need for the m_is_active field.
class ConcurrentClusterVertex {

public:

  ConcurrentClusterVertex(): m_center(0), m_is_center(false), m_probability(-1.0) {};

  bool is_covered() const {
    return m_probability > 0.0;
  }

  ugraph_vertex_t center() const {
    REQUIRE(is_covered(), "The node is uncovered");
    return m_center;
  }

  bool is_center() const {
    return m_is_center;
  }

  probability_t probability() const {
    REQUIRE(is_covered(), "The node is uncovered");
    return m_probability;
  }

  void make_center(ugraph_vertex_t id) {
    assert(!is_covered());
    m_center = id;
    m_is_center = true;
    m_probability = 1.0;
  }

  void cover(ugraph_vertex_t center, probability_t prob) {
    REQUIRE(!is_covered(), "The node is covered");
    m_center = center;
    m_probability = prob;
  }
  
private:
  
  ugraph_vertex_t m_center;

  bool m_is_center;

  probability_t m_probability;

};

std::vector< ConcurrentClusterVertex > concurrent_cluster(const ugraph_t & graph,
                                                          CCSampler & sampler,
                                                          const size_t batch,
                                                          const probability_t p_low,
                                                          Xorshift1024star & rnd,
                                                          ExperimentReporter & experiment);
