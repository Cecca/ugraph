#pragma once

#include "types.hpp"
#include "require.hpp"

class ClusterVertex {

public:

  ClusterVertex(): m_center(0), m_is_center(false), m_probability(-100.0) {};

  bool is_covered() const {
    return m_probability >= 0.0;
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

  probability_t unreliable_probability() const {
    REQUIRE(!is_covered(), "The node is covered, its probability is reliable");
    if (m_probability < -1.0) {
      return 0.0;
    } else {
      return -m_probability;
    }
  }
  
  void make_center(ugraph_vertex_t id) {
    REQUIRE(!is_covered(), "The node is already covered");
    m_center = id;
    m_is_center = true;
    m_probability = 1.0;
  }

  void force_make_center(ugraph_vertex_t id) {
    m_center = id;
    m_is_center = true;
    m_probability = 1.0;
  }
  
  void cover(ugraph_vertex_t center, probability_t prob) {
    m_center = center;
    m_probability = prob;
  }

  void unreliable_cover(ugraph_vertex_t center, probability_t prob) {
    m_center = center;
    m_probability = -prob;
  }
  
private:
  
  ugraph_vertex_t m_center;

  bool m_is_center;

  probability_t m_probability;

};
