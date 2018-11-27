#include "logging.hpp"
#include "termcolor.hpp"

std::string logging::vm_size() {
  std::ifstream info("/proc/self/status");
  std::string line;
  while (std::getline(info, line)) {
    if (line.find("VmSize") != std::string::npos){
      break;
    }
  }
  info.close();
  std::vector<std::string> tokens;
  boost::split(tokens, line, boost::is_space());
  size_t ntok = tokens.size();
  if (ntok < 2) {
    return "xxx Kb";
  }
  try {
    int num_kb = stoi(tokens[ntok - 2]);
    std::stringstream sstr;
    if (num_kb > 1048576) {
      sstr << std::setprecision(3) << (num_kb / 1048576.0) << " Gb";
    } else if (num_kb > 1024) {
      sstr << std::setprecision(3) << (num_kb / 1024.0) << " Mb";
    } else {
      sstr << num_kb << " Kb";
    }
    return sstr.str();
  } catch (const std::invalid_argument& e) {
    return "xxx Kb";
  }
}  

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
              /* << "[" << vm_size() << "] " */
              << msg
              << std::endl;
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
