#include "prelude.hpp"
#include "rand.hpp"

int main(int argc, char**argv) {
  Xorshift1024star rnd(11091872049813);

  std::cout << "First sequence" << std::endl;
  for(size_t i=0; i<3; ++i) {
    std::cout << rnd.next_double() << std::endl;
  }

  Xorshift1024star rnd2 = rnd;
  rnd2.jump();
  std::cout << "Second sequence" << std::endl;
  for(size_t i=0; i<3; ++i) {
    std::cout << rnd2.next_double() << std::endl;
  }

  std::cout << "First sequence?" << std::endl;
  for(size_t i=0; i<3; ++i) {
    std::cout << rnd.next_double() << std::endl;
  }

  std::cout << "Reference" << std::endl;
  Xorshift1024star rnd_ref(1);
  for(size_t i=0; i<6; ++i) {
    std::cout << rnd_ref.next_double() << std::endl;
  }
  
  return 0;
}
