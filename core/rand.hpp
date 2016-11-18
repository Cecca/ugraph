#pragma once

#include <stdint.h>
#include "prelude.hpp"

/**
 * Taken from http://xoroshiro.di.unimi.it/splitmix64.c
 */
class Splitmix64 {
private:
  uint64_t x;

public:
  Splitmix64(uint64_t seed): x(seed) {};

  uint64_t next() {
    uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
  }

};


/**
 * Taken from http://xoroshiro.di.unimi.it/xorshift1024star.c
 */
class Xorshift1024star {
private:
  uint64_t s[16]; 
  int p;

public:

  Xorshift1024star(uint64_t seed): p(0) {
    Splitmix64 seeder(seed);
    for (size_t i = 0; i < 16; ++i) {
      s[i] = seeder.next();
    }
  }

  std::array<uint64_t, 16> state() const {
    std::array<uint64_t, 16> to_ret;
    for (size_t i=0; i<16; ++i) {
      to_ret[i] = s[i];
    }
    return to_ret;
  }
  
  uint64_t next() {
    const uint64_t s0 = s[p];
    uint64_t s1 = s[p = (p + 1) & 15];
    s1 ^= s1 << 31; // a
    s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
    return s[p] * UINT64_C(1181783497276652981);
  }

  /**
   * From http://xoroshiro.di.unimi.it/:
   *
   * With the exception of generators designed to provide directly
   * double-precision floating-point numbers, the fastest way to
   * convert in C99 a 64-bit unsigned integer x to a 64-bit double is
   * presented in this function.
   *
   * The code above cooks up by bit manipulation a real number in the
   * interval [1..2), and then subtracts one to obtain a real number
   * in the interval [0..1). If x is chosen uniformly among 64-bit
   * integers, d is chosen uniformly among dyadic rationals of the
   * form k / 2^{-52}.
   */
  double next_double() {
    uint64_t x = next();
    const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
  }

  void jump() {
    static const uint64_t JUMP[] = { 0x84242f96eca9c41d,
                                     0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
                                     0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
                                     0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
                                     0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
                                     0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
    };

    uint64_t t[16] = { 0 };
    for(size_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
      for(size_t b = 0; b < 64; b++) {
        if (JUMP[i] & 1ULL << b)
          for(size_t j = 0; j < 16; j++)
            t[j] ^= s[(j + p) & 15];
        next();
      }
    
    for(int j = 0; j < 16; j++)
      s[(j + p) & 15] = t[j];
  }

};


namespace std {
  std::ostream & operator<<(std::ostream &os, Xorshift1024star & rnd);
}
