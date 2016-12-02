#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "bfs_sampler.hpp"
#include "logging.hpp"
#include "sequential_clustering.hpp"


int main(int argc, char**argv) {
  FixedCapacityQueue<int> queue(12);
  int elems[] = {1,2,3,4,5,6,7,8,9,8};
  for (int e : elems) {
    queue.push(e);
  }
  std::cout << queue.pop() << std::endl;
  std::cout << queue.pop() << std::endl;
  queue.push(1001);
  queue.push(1002);
  while(!queue.empty()) {
    std::cout << queue.pop() << std::endl;
  }
}
