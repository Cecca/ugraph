#pragma once

#include "prelude.hpp"
#include "types.hpp"
#include "require.hpp"
#include "logging.hpp"
#include "cluster_vertex.hpp"

// A cache storing connection counts from a set of centers to all the
// others. Has a maximum size, and centers are evicted on an LRU basis

struct ConnectionCountsCacheElement {

  ConnectionCountsCacheElement(): times_accessed(0), num_samples(0),
                                  counts(std::vector<size_t>()) {
    // It's here only to satisfy the compiler
    throw std::logic_error("Never use this constructor");
  }

  ConnectionCountsCacheElement(size_t n): times_accessed(0), num_samples(0),
                                          counts(std::vector<size_t>(n)) {}
  
  size_t times_accessed;
  size_t num_samples;
  std::vector<size_t> counts;
};

class ConnectionCountsCache {

public:
  ConnectionCountsCache(size_t max_size)
    : m_max_size(max_size),
      m_cache(std::unordered_map<ugraph_vertex_t, ConnectionCountsCacheElement>()) {}

  bool contains(ugraph_vertex_t v) const {
    return m_cache.count(v) > 0;
  }

  ConnectionCountsCacheElement& get(ugraph_vertex_t v) {
    REQUIRE(contains(v), "Cache does not contain the requested element");
    m_cache[v].times_accessed++;
    return m_cache[v];
  }

  void add_new(ugraph_vertex_t v, size_t n) {
    REQUIRE(!contains(v), "Cache already contains the requested element");
    LOG_DEBUG("Adding " << v << " to the cache");
    m_cache.emplace(v, ConnectionCountsCacheElement(n));
  }

  ConnectionCountsCacheElement& get_or_new(ugraph_vertex_t v, size_t n) {
    if (contains(v)) {
      return get(v);
    } else {
      add_new(v, n);
      return get(v);
    }
  }

  size_t size() const {
    return m_cache.size();
  }

  std::string str() const {
    std::stringstream sstr;
    sstr << "Cache size: " << m_cache.size() << " ::: ";
    for (auto it=m_cache.begin(); it != m_cache.end(); it++) {
      sstr << it->first << "(" << it->second.times_accessed << ") ";
    }
    return sstr.str();
  }

  int uncovered_node(const std::vector< ClusterVertex > & vinfo) {
    size_t n = vinfo.size();
    for (int i=0; i<n; i++) {
      if (m_cache.count(i) > 0) {
        if (!vinfo[i].is_covered()) {
          return i;
        } else if (vinfo[i].is_covered() && !vinfo[i].is_center()) {
          // reset counter to mark for eviction for nodes that are
          // covered but are not centers.
          m_cache[i].times_accessed = 0;
        }
      }
    }
    return -1;
  }
  
  void cleanup() {
    while (m_cache.size() > m_max_size) {
      // remove the cache entry with the lowest times_accessed count
      size_t min_accessed = std::numeric_limits<size_t>::max();
      ugraph_vertex_t min_v = -1;
      for (auto it=m_cache.begin(); it != m_cache.end(); it++) {
        if (it->second.times_accessed < min_accessed) {
          min_accessed = it->second.times_accessed;
          min_v = it->first;
        }
      }
      LOG_DEBUG("Removing from cache vertex " << min_v <<
                " which was accessed " << min_accessed << " times");
      m_cache.erase(min_v);
    }
  }
  
private:
  size_t m_max_size;
  
  std::unordered_map<ugraph_vertex_t, ConnectionCountsCacheElement> m_cache;
  
};
