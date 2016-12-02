#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "cc_sampler.hpp"
#include "logging.hpp"
#include "sequential_clustering.hpp"


size_t prob_to_samples(probability_t prob, double epsilon, double delta) {
  return 1/(epsilon*epsilon*prob) * log(1/delta);
}

int main(int argc, char**argv) {
  
}
