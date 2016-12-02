#pragma once

// Code in this header is responsible of handling progressive
// sampling, so to answer queries on connection probabilities.
// In particular, this implements limited-depth sampling

#include "logging.hpp"
#include "prelude.hpp"
#include "rand.hpp"
#include "types.hpp"

template <typename T> class FixedCapacityQueue {

public:
  FixedCapacityQueue(size_t capacity)
      : m_capacity(capacity), m_begin(0), m_end(0),
        m_storage(std::vector<T>(capacity)) {}

  bool empty() const { return m_begin == m_end; }

  void clear() {
    m_begin = 0;
    m_end = 0;
  }
  
  void push(T elem) {
    size_t new_end = (m_end + 1) % m_capacity;
    if (new_end == m_begin) {
      throw std::logic_error("Queue capacity exceeded");
    }
    m_storage[m_end] = elem;
    m_end = new_end;
  }

  T pop() {
    assert(!empty());
    T elem = m_storage[m_begin];
    m_begin = (m_begin + 1) % m_capacity;
    return elem;
  }

private:
  size_t m_capacity;
  size_t m_begin;
  size_t m_end;

  std::vector<T> m_storage;
};

struct BfsSamplerThreadState {
  typedef std::vector<size_t> distance_vector_t;
  typedef std::vector<bool> edge_sample_t;

  BfsSamplerThreadState(const ugraph_t &graph, Xorshift1024star randgen)
    : queue(FixedCapacityQueue<ugraph_vertex_t>(boost::num_vertices(graph))),
      connection_counts(std::vector<size_t>(boost::num_vertices(graph))),
      distance_vector(std::vector<size_t>(boost::num_vertices(graph))),
      rnd(randgen){};

  FixedCapacityQueue<ugraph_vertex_t> queue;

  std::vector<size_t> connection_counts;

  std::vector<size_t> distance_vector;

  // Random generators, one for each thread
  Xorshift1024star rnd;
};

namespace std {
std::ostream &operator<<(std::ostream &os, BfsSamplerThreadState &tstate);
}

class BfsSampler {

public:
  typedef std::vector<bool> edge_sample_t;

  BfsSampler(const ugraph_t &graph,
             const size_t max_depth,
             std::function<size_t(double)> prob_to_samples, uint64_t seed,
             size_t num_threads)
      : prob_to_samples(prob_to_samples),
        m_max_dist(max_depth),
        m_samples(std::vector<edge_sample_t>()),
        m_thread_states(std::vector<BfsSamplerThreadState>()) {
    Xorshift1024star rnd(seed);
    for (size_t i = 0; i < num_threads; ++i) {
      rnd.jump();
      m_thread_states.emplace_back(graph, rnd);
    }
  }

  void min_probability(const ugraph_t &graph, probability_t prob);

  // Add samples, if needed
  void sample_size(const ugraph_t &graph, size_t total_samples);

  void log_states() {
    for (auto &tstate : m_thread_states) {
      LOG_INFO(tstate);
    }
  }

  size_t connection_probabilities(const ugraph_t &graph,
                                  const ugraph_vertex_t from,
                                  std::vector<probability_t> &probabilities);

  /** The probability that a given set of nodes is connected */
  probability_t
  connection_probability(const ugraph_t &graph,
                         const std::vector<ugraph_vertex_t> &vertices) {
    throw std::logic_error("Not implemented");
  }

private:
  std::function<size_t(double)> prob_to_samples;

  size_t m_max_dist;
  
  std::vector<edge_sample_t> m_samples;

  std::vector<BfsSamplerThreadState> m_thread_states;

  // The minimum connection probability that is estimate reliably
  probability_t m_min_probability = 1.0;
};
