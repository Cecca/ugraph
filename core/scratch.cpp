#include "prelude.hpp"
#include "guesser.hpp"


int main(int argc, char**argv) {
  if (argc < 2) std::exit(1);
  probability_t target = std::stod(std::string(argv[1]));

  LOG_INFO("Target is " << target);
  
  ExponentialGuesser g(0.01, 0.001);

  probability_t guess;
  while(!g.stop()) {
    guess = g.guess();
    LOG_INFO("" << guess);
    if (guess <= target) {
      g.below();
    } else {
      g.above();
    }
  }

  LOG_INFO(">>> Final guess " << guess);
}
