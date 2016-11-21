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
#include "logging.hpp"

struct CCSamplerThreadState {
  typedef std::vector< int > component_vector_t;
  typedef std::vector< ugraph_vertex_t > stack_t;
  typedef std::vector< bool > edge_sample_t;

  CCSamplerThreadState(const ugraph_t & graph,
                       Xorshift1024star randgen) : stack(stack_t(boost::num_vertices(graph))),
                                                   edge_sample(edge_sample_t(boost::num_edges(graph))),
                                                   connection_counts(std::vector< size_t >(boost::num_vertices(graph))),
                                                   rnd(randgen)
  {};
  
  
  // Stacks for dfs, one for each thread
  stack_t stack;

  // Edge samples, one for each thread
  edge_sample_t edge_sample;

  std::vector< size_t > connection_counts;

  // Random generators, one for each thread
  Xorshift1024star rnd;
  
};

namespace std {
  std::ostream & operator<<(std::ostream &os, CCSamplerThreadState & tstate);
}

/**
 * A sampler representing each sample as a connected component
 */
class CCSampler {

public:

  typedef std::vector< int > component_vector_t;
  
  CCSampler(const ugraph_t & graph,
            double epsilon,
            double delta,
            std::function<size_t(double, double, double)> prob_to_samples,
            uint64_t seed,
            size_t num_threads)
    : m_epsilon(epsilon),
      m_delta(delta),
      prob_to_samples(prob_to_samples),
      m_samples(std::vector< component_vector_t >()),
      m_thread_states(std::vector< CCSamplerThreadState >()) {
    Xorshift1024star rnd(seed);
    for (size_t i=0; i<num_threads; ++i) {
      rnd.jump();
      m_thread_states.emplace_back(graph, rnd);
    }
  }

  void min_probability(const ugraph_t & graph, probability_t prob);

  // Add samples, if needed
  void sample_size(const ugraph_t & graph, size_t total_samples);

  void log_states() {
    for (auto & tstate : m_thread_states) {
      LOG_INFO(tstate);
    }
  }
  
  size_t connection_probabilities(const ugraph_t & graph,
                                  const ugraph_vertex_t from,
                                  std::vector< probability_t > & probabilities);
 
private:

  const double m_epsilon;

  const double m_delta;

  std::function<size_t(double, double, double)> prob_to_samples;

  std::vector< component_vector_t > m_samples;

  std::vector< CCSamplerThreadState > m_thread_states;

  // The minimum connection probability that is estimate reliably
  probability_t m_min_probability = 1.0;
};
