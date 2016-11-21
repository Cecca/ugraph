#pragma once

#include <exception>

#define REQUIRE(condition, message)             \
  if(!(condition)){                              \
    throw std::logic_error(message);             \
  }
