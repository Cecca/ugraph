#include "logging.hpp"
#include "termcolor.hpp"

logging::Level log_level = logging::Level::Info;

void logging::set_level(logging::Level level) {
  log_level = level;
}

void logging::trace(std::string & msg) {
  if (log_level <= logging::Level::Trace){
    std::cerr << termcolor::grey << "- " << termcolor::reset
              << msg << std::endl;
  }
}

void logging::debug(std::string & msg) {
  if (log_level <= logging::Level::Debug){
    std::cerr << termcolor::cyan << "= " << termcolor::reset
              << msg << std::endl;
  }
}

void logging::info(std::string & msg) {
  if (log_level <= logging::Level::Info){
    std::cerr << termcolor::green << "* " << termcolor::reset
              << msg << std::endl;
  }
}

void logging::warn(std::string & msg) {
  if (log_level <= logging::Level::Warn){
    std::cerr << termcolor::yellow << "! " << termcolor::reset
              << msg << std::endl;
  }
}

void logging::error(std::string & msg) {
  if (log_level <= logging::Level::Error){
    std::cerr << termcolor::red << "X " << termcolor::reset
              << msg << std::endl;
  }
}
