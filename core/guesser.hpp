#pragma once

#include "logging.hpp"
#include "types.hpp"
#include "require.hpp"

class Guesser {

public:
  Guesser(const probability_t gamma, const probability_t p_low)
      : m_p_low(p_low), m_gamma(gamma), m_upper(1.0), m_lower(1.0),
        binary_search(false), m_i(0) {}

  /** Set the state of the guesser when the guess is below the solution */
  void below() {
    if (binary_search) {
      m_lower = (m_upper + m_lower) / 2;
      LOG_DEBUG("[Guesser] Below in binary search (new m_lower" << m_lower << ")");
    } else {
      LOG_DEBUG("[Guesser] Below in exponential search: start binary search");
      binary_search = true; // Start binary search on the next call of `guess`
    }
  }

  /** Set the state of the guesser when the guess is above the solution */
  void above() {
    if (binary_search) {
      m_upper = (m_upper + m_lower) / 2;
      LOG_DEBUG("[Guesser] Above in binary search (new m_upper" << m_lower << ")");
    } else {
      m_upper = m_lower;
      m_lower = 1.0 - m_gamma * (1 << m_i);
      m_i++;
      if (m_lower <= m_p_low) {
        m_lower = m_p_low;
        binary_search = true;
      }
      LOG_DEBUG("[Guesser] Above in exponential search. New m_lower " << m_lower);
    }
  }

  probability_t guess() {
    REQUIRE(m_lower <= m_upper, "Upper and lower bounds inverted!");
    if (binary_search) {
      // The last guess is the lower bound, so to have always a
      // probability less than the necessary one
      // if (stop()) { return m_lower; }
      return (m_upper + m_lower) / 2;
    } else {
      return m_lower;
    }
  }

  bool stop() const {
    if (binary_search) {
      LOG_DEBUG("[Guesser] Stopping condition: " << (1.0-m_lower/m_upper) << "<=" << m_gamma);
    }
    return binary_search && (1.0-m_lower/m_upper) <= m_gamma;
  }

private:
  const probability_t m_p_low;
  probability_t m_gamma;
  probability_t m_upper;
  probability_t m_lower;

  bool binary_search;
  size_t m_i;
};


class GeometricGuesser {
  public:
  GeometricGuesser(const probability_t gamma, const probability_t p_low)
    : m_p_low(p_low), m_gamma(gamma), m_guess(1.0), m_below(false) {}

  /** Set the state of the guesser when the guess is below the solution */
  void below() {
    m_below = true;
  }

  /** Set the state of the guesser when the guess is above the solution */
  void above() {
    m_guess *= m_gamma;
  }

  probability_t guess() {
    return m_guess;
  }

  bool stop() const {
    return m_below || m_guess < m_p_low;
  }

private:
  const probability_t m_p_low;
  probability_t m_gamma;
  probability_t m_guess;
  bool m_below;

};
