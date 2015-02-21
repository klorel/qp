/*
 * ConvexRelaxation.h
 *
 *  Created on: 21 févr. 2015
 *      Author: mruiz
 */

#ifndef SRC_CONVEXRELAXATION_H_
#define SRC_CONVEXRELAXATION_H_

class Problem;

class ConvexRelaxation {
 public:
  ConvexRelaxation(Problem const &);
  virtual ~ConvexRelaxation();
 public:
  void identifyBlocks();
 private:
  Problem const & _problem;
};

#endif /* SRC_CONVEXRELAXATION_H_ */
