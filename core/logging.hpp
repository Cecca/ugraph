#pragma once

#include "prelude.hpp"

namespace logging {

  enum Level : int {
    Trace = 0,
    Debug = 1,
    Info  = 2,
    Warn  = 3,
    Error = 4
  };
  
  void set_level(logging::Level level);
  
  void error(std::string & msg);
  void warn(std::string & msg);
  void info(std::string & msg);
  void debug(std::string & msg);
  void trace(std::string & msg);

  std::string vm_size();
  
}


#define LOG_TRACE(x)\
  do { std::stringstream s; \
    s << x; \
    std::string msg = s.str(); \
    logging::trace(msg);\
  } while (false);


#define LOG_DEBUG(x)\
  do { std::stringstream s; \
    s << x; \
    std::string msg = s.str(); \
    logging::debug(msg);\
  } while (false);


#define LOG_INFO(x)\
  do { std::stringstream s; \
    s << x; \
    std::string msg = s.str(); \
    logging::info(msg);\
  } while (false);


#define LOG_WARN(x)\
  do { std::stringstream s; \
    s << x; \
    std::string msg = s.str(); \
    logging::warn(msg);\
  } while (false);


#define LOG_ERROR(x)\
  do { std::stringstream s; \
    s << x; \
    std::string msg = s.str(); \
    logging::error(msg);\
  } while (false);

#define PRINT_CLUSTERING(clustering)                                           \
  printf("[");                                                                 \
  for (ugraph_vertex_t ___j = 0; ___j < n; ___j++) {                           \
    if (clustering[___j].is_covered()) {                                       \
      printf("+");                                                             \
    } else {                                                                   \
      printf("-");                                                             \
    }                                                                          \
  }                                                                            \
  printf("]\n");                                                               \



// End of logging.hpp
