/*
 * ConvexRelaxation.cpp
 *
 *  Created on: 21 févr. 2015
 *      Author: mruiz
 */

#include "ConvexRelaxation.h"
#include "Problem.h"

ConvexRelaxation::ConvexRelaxation(Problem const & problem)
    : _problem(problem) {
  // TODO Auto-generated constructor stub

}

ConvexRelaxation::~ConvexRelaxation() {
  // TODO Auto-generated destructor stub
}

void ConvexRelaxation::identifyBlocks() {
  QuadraticTerms qTerms;
  _problem.getQuadraticIndexes(qTerms);
// get the partition induced by the quadratic terms
  std::vector<std::list<Variable>> links(_problem.variables().size());
  for (auto const & v : _problem.variables()) {
    links[v.id()].push_back(v);
  }
  for (auto const & term : qTerms) {
    auto & list1(links[term.first.first.id()]);
    auto & list2(links[term.first.second.id()]);
    list1.splice(list1.end(), list2);
  }
  std::vector<size_t> blockId(links.size(), -1);
  size_t nBlock(0);
  size_t maxSize(0);
  for (auto const & list : links) {
    if (!list.empty())
      std::cout << "block ";
    for (auto const & v : list) {
      std::cout << v.name() << " ";
      blockId[v.id()] = nBlock;
    }
    if (!list.empty()) {
      std::cout << std::endl;
      ++nBlock;
      maxSize = std::max(maxSize, list.size());
    }
  }
  std::cout << "nBlock  : " << nBlock << std::endl;
  std::cout << "maxSize : " << maxSize << std::endl;
}
