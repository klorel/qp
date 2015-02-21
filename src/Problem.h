#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <cassert>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <utility>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <memory>

double const MyAmplInf = 1e20;
double const MyEps = 1e-10;
double const MyInfinity = 1e50;

enum SENSE {
  LEQ,
  EQ,
  GEQ
};

class QuadraticIndex;
class Variable;
class Function;
class Constraint;
class Problem;

typedef std::map<QuadraticIndex, double> QuadraticTerms;
typedef std::map<Variable, double> LinearTerms;

typedef std::shared_ptr<Variable> VariablePtr;
typedef std::shared_ptr<QuadraticTerms> QuadraticTermsPtr;
typedef std::shared_ptr<LinearTerms> LinearTermsPtr;
typedef std::shared_ptr<double> DoublePtr;
typedef std::shared_ptr<SENSE> SensePtr;
typedef std::shared_ptr<std::string> StringPtr;
typedef std::shared_ptr<size_t> IndexPtr;

typedef std::vector<Function> Objectifs;
typedef std::vector<Variable> Variables;
typedef std::vector<Constraint> Constraints;

Function operator-(Variable const &);
Function operator-(Function const &);

Function operator+(Variable const &);
Function operator+(Function const &);

Function operator+(double, Variable const &);
Function operator+(Variable const &, double);
Function operator+(Variable const &, Variable const &);

Function operator-(double, Variable const &);
Function operator-(Variable const &, double);
Function operator-(Variable const &, Variable const &);

Function operator*(double, Variable const &);
Function operator*(Variable const &, double);
Function operator*(Variable const &, Variable const &);

Function operator+(double, Function const &);
Function operator+(Function const &, double);
Function operator+(Function const &, Function const &);

Function operator-(double, Function const &);
Function operator-(Function const &, double);
Function operator-(Function const &, Function const &);

Function operator*(double, Function const &);
Function operator*(Function const &, double);
Function operator*(Function const &, Function const &);

Function operator+(Variable const &, Function const &);
Function operator+(Function const &, Variable const &);

Function operator-(Variable const &, Function const &);
Function operator-(Function const &, Variable const &);

Function operator*(Variable const &, Function const &);
Function operator*(Function const &, Variable const &);

Constraint operator<=(Function const &, double);
Constraint operator==(Function const &, double);
Constraint operator>=(Function const &, double);

Constraint operator<=(double, Function const &);
Constraint operator==(double, Function const &);
Constraint operator>=(double, Function const &);

std::ostream & operator<<(std::ostream &, SENSE);
std::ostream & operator<<(std::ostream &, Variable const &);
std::ostream & operator<<(std::ostream &, QuadraticIndex const &);
std::ostream & operator<<(std::ostream &, QuadraticTerms const &);
std::ostream & operator<<(std::ostream &, LinearTerms const &);
std::ostream & operator<<(std::ostream &, Function const &);
std::ostream & operator<<(std::ostream &, Constraint const &);
std::ostream & operator<<(std::ostream &, Problem const &);

void operator<<(Problem &, Constraint const &);
void operator>>(Constraint const &, Problem&);

class Variable {
 public:
  std::string const & name() const;
  std::string & name();
  double lb() const;
  double &lb();
  double ub() const;
  double &ub();

  size_t id() const;
  explicit Variable(size_t);
  Variable(Variable const &);
  bool operator<(Variable const & rhs) const {
    return id() < rhs.id();
  }
 private:
  IndexPtr _id;
  StringPtr _name;
  DoublePtr _lb;
  DoublePtr _ub;
};
class QuadraticIndex : public std::pair<Variable, Variable> {
 public:
  QuadraticIndex(Variable const & lhs, Variable const & rhs)
      : std::pair<Variable, Variable>(lhs.id() < rhs.id() ? lhs : rhs,
                                      lhs.id() < rhs.id() ? rhs : lhs) {
  }
  void print(std::ostream &) const;
};

class Function {
 public:
  void add(Variable const &, Variable const &, double);
  void add(QuadraticIndex, double);
  void add(Variable const &, double);
  void add(double);
  void clear();

  void operator+=(Variable const &);
  void operator+=(Function const &);
  void operator+=(double);

  void operator-=(Variable const &);
  void operator-=(Function const &);
  void operator-=(double);

  void operator*=(Variable const &);
  void operator*=(Function const &);
  void operator*=(double);

  bool isCst() const;
  bool isLinear() const;

  void print(std::ostream & stream) const;
 public:
  Function();
  Function(Function const & rhs);
  virtual ~Function();
 public:
  QuadraticTerms const & qTerms() const;
  LinearTerms const & lTerms() const;
  double const & cst() const;
 private:
  QuadraticTerms & qTerms();
  LinearTerms & lTerms();
  double & cst();
 private:
  QuadraticTermsPtr _qTerms;
  LinearTermsPtr _lTerms;
  DoublePtr _cst;
};

class Constraint : public Function {
 public:
  void setSense(SENSE);
  void print(std::ostream & stream) const;

  Constraint();
  Constraint(Constraint const & rhs);
  Constraint(Function const &, SENSE, double);
  Constraint(Function const &, double, double);
  virtual ~Constraint();
  SENSE sense() const;
 protected:
  SENSE &sense();
 private:
  SensePtr _sense;
};
class Problem {
 public:
  Variable const & newVariable(std::string const &);
  Variable const & newVariable(double, double);
  void print(std::ostream &) const;
  void newConstraint(Constraint const &);

  size_t addObj(Function const &);

  Variable & variable(size_t idx);

  Objectifs & objectifs();
  Objectifs const & objectifs() const;

  Variables & variables();
  Variables const & variables() const;

  Constraints & constraints();
  Constraints const & constraints() const;

  void clear();

  void getQuadraticIndexes(QuadraticTerms &) const;
 protected:
  Objectifs _objectifs;
  Variables _variables;
  Constraints _constraints;
};
#endif
