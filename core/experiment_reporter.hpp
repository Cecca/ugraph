#pragma once

#include "json.hpp"
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/variant.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <map>
#include <exception>

struct ToJsonVisitor : public boost::static_visitor<nlohmann::json> {

  template <typename T> nlohmann::json operator()(T &i) const {
    return nlohmann::json(i);
  }
};

class ExperimentReporter {

  typedef boost::variant<size_t, int, double, std::string, bool>
      element_value_t;
  typedef std::map<std::string, element_value_t> table_row_t;
  typedef std::vector<table_row_t> table_t;

public:
  ExperimentReporter() : date(boost::posix_time::second_clock::local_time()) {}

public:
  static ExperimentReporter& get_instance() {
    static ExperimentReporter instance;
    return instance;
  }
  
  void tag(const std::string &key, const element_value_t &val) {
    tags[key] = val;
  }

  void append(const std::string &table, table_row_t row) {
    if (tables.find(table) == tables.end()) {
      tables[table] = table_t();
    } else {
      table_row_t first_row = tables[table].front();
      std::set<std::string> keys;
      for (auto it = first_row.begin(); it != first_row.end(); ++it) {
        keys.insert(it->first);
      }
      if (row.size() != keys.size()) {
        throw std::runtime_error(
            "New row and existing ones have a different number of columns");
      }
      for (auto it = row.begin(); it != row.end(); ++it) {
        if (keys.find(it->first) == keys.end()) {
          throw std::runtime_error(
              "New row has different headings from the existing ones");
        }
      }
    }
    tables[table].push_back(row);
  }

  void save(std::ostream &out) {
    nlohmann::json root;
    
    root["date"] = boost::posix_time::to_iso_extended_string(date);

    ToJsonVisitor jvis;

    for (auto it = tags.begin(); it != tags.end(); ++it) {
      root["tags"][it->first] = boost::apply_visitor(jvis, it->second);
    }

    for (auto tables_it = tables.begin(); tables_it != tables.end();
         ++tables_it) {
      std::string tname = tables_it->first;
      table_t table = tables_it->second;
      nlohmann::json jtable;
      for (auto row_it = table.begin(); row_it != table.end(); ++row_it) {
        nlohmann::json jrow;
        for (auto it = row_it->begin(); it != row_it->end(); ++it) {
          jrow[it->first] = boost::apply_visitor(jvis, it->second);
        }
        jtable.push_back(jrow);
      }
      root["tables"][tname] = jtable;
    }

    out << root.dump(2) << std::endl;
  }

  void save(const std::string &path) {
    std::ofstream out_stream;
    out_stream.open(path + ".bz2");
    boost::iostreams::filtering_ostream out;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(out_stream);
    save(out);
  }

  void save() {
    std::string path = boost::posix_time::to_iso_string(date);
    path += ".json";
    save(path);
  }

private:
  boost::posix_time::ptime date;
  std::map<std::string, element_value_t> tags;
  std::map<std::string, table_t> tables;
};


#define EXPERIMENT_APPEND(table, ...)                           \
  ExperimentReporter::get_instance().append(table, __VA_ARGS__)

#define EXPERIMENT_TAG(key, val)                        \
  ExperimentReporter::get_instance().tag(key, val)

#define EXPERIMENT_SAVE()                       \
  ExperimentReporter::get_instance().save()
