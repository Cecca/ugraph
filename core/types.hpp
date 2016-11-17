#pragma once

#include "prelude.hpp"

typedef double probability_t;

struct EdgeData {

  probability_t probability;
  size_t index;

  EdgeData() : probability(1.0), index(0) {}

  EdgeData(probability_t p, size_t idx)
    : probability(p), index(idx) {}
  
};

class EdgeDataFactory {
private:
  size_t m_idx;

public:
  EdgeDataFactory(): m_idx(0) {};

  EdgeData build(probability_t p) {
    return EdgeData(p, m_idx++);
  }
};


struct VertexData {

  std::string label;

  VertexData(std::string label) : label(label) {}

  VertexData(): label("no label") {}

};

/**
 * An uncertain graph has weighted indexed edges, with probabilities attached.
 */
typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  VertexData,
  EdgeData
  > ugraph_t;

/**
 * Uncertain graph vertex type
 */
typedef boost::graph_traits<ugraph_t>::vertex_descriptor
ugraph_vertex_t;

/**
 * Uncertain graph edge type
 */
typedef boost::graph_traits<ugraph_t>::edge_descriptor
ugraph_edge_t;
