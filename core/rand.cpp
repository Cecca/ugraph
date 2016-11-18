#include "rand.hpp"

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

