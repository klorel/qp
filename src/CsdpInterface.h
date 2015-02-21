#ifndef CSDP_INTERFACE_H_
#define CSDP_INTERFACE_H_

#include <map>
#include <vector>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Problem.h"

class CsdpInterface {
 public:
  CsdpInterface();
  ~CsdpInterface();
 public:
  void clear();
  void add(size_t, double);
  void add(size_t, size_t, double);
  void bestReformulation(std::map<size_t, double> & result) const;
 private:
  size_t id(size_t);
  void print(std::ostream & stream = std::cout) const;
  std::string const & name(size_t id) const;
 private:
};
#endif
