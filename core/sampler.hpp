#pragma once

// Code in this header is responsible of handling progressive
// sampling, so to answer queries on connection probabilities.
//
// The code is parametric in the representation of samples:
//  - connected components
//  - edge flags

#include "prelude.hpp"
#include "types.hpp"
#include "rand.hpp"

/**
 * A sampler representing each sample as a connected component
 */
class CCSampler {

  typedef std::vector< int > component_vector_t;
  typedef std::vector< ugraph_vertex_t > stack_t;
  typedef std::vector< bool > edge_sample_t;

public:

  CCSampler(ugraph_t & graph,
            uint64_t seed,
            size_t num_threads)
    : m_samples(std::vector< component_vector_t >()),
      s_stacks(std::vector< stack_t >(num_threads)),
      s_edge_samples(std::vector< edge_sample_t >(num_threads)),
      m_rnds(std::vector< Xorshift1024star >())
  {
    for (size_t i=0; i<num_threads; ++i) {
      s_stacks[i].reserve(boost::num_vertices(graph));
      s_edge_samples[i].reserve(boost::num_edges(graph));
    }
    Xorshift1024star rnd(seed);
    for (size_t i=0; i<num_threads; ++i) {
      rnd.jump();
      m_rnds.push_back(rnd);
    }
  }

  // Add samples, if needed
  void set_sample_size(ugraph_t & graph, size_t total_samples);

 
private:
  std::vector< component_vector_t > m_samples;

  // Stacks for dfs, one for each thread
  std::vector< stack_t > s_stacks;

  // Edge samples, one for each thread
  std::vector< edge_sample_t > s_edge_samples;

public:
  // Random generators, one for each thread
  std::vector< Xorshift1024star > m_rnds;
  
};
