#pragma once

#include "types.hpp"
#include "rand.hpp"

void dfs(const ugraph_t & graph,
         const std::vector< bool > & smpl,
         std::vector< int > & components_map,
         std::vector< ugraph_vertex_t > & stack,
         ugraph_vertex_t root,
         int component_id);

void union_find(const ugraph_t & graph,
                Xorshift1024star & rnd,
                std::vector< size_t > & ranks,
                std::vector< int > & roots);

void connected_components(const ugraph_t & graph,
                          const std::vector< bool > & smpl,
                          std::vector< int > & components_map,
                          std::vector< ugraph_vertex_t > & stack);
