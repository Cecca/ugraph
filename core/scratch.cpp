#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "sampler.hpp"

namespace std {
  std::ostream &
  operator<<(std::ostream &os,
             Xorshift1024star & rnd)
  {
    std::stringstream ss;
    ss << "Xorshift1024*(0x" << std::uppercase << std::setfill('0') << std::setw(4) << std::hex;
    auto state = rnd.state();
    for(size_t i=0; i<2; i++) {
      ss << state[i];
    }
    ss << ")";
    os << ss.str();
    return os;
  }
}


int main(int argc, char**argv) {
  std::string graph_path(argv[1]);
  ugraph_t graph;
  read_edge_list(graph, graph_path);
  std::cout << "Num nodes " << boost::num_vertices(graph)
            << " num edges " << boost::num_edges(graph)
            << std::endl;

  auto omp_threads = omp_get_max_threads();
  std::cout << "Running with " << omp_threads << " threads" << std::endl;
  
  CCSampler sampler(graph, 123, omp_threads);
  for(Xorshift1024star & rnd : sampler.m_rnds) {
    std::cout << rnd << std::endl;
  }
  sampler.set_sample_size(graph, 128);
  std::cout << "==================" << std::endl;
  for(Xorshift1024star & rnd : sampler.m_rnds) {
    std::cout << rnd << std::endl;
  }
}
