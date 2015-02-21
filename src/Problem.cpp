#include "Problem.h"
#include "Formatter.h"
#include <cassert>
/**
 *
 */
Variable::Variable(Variable const & rhs)
    : _id(rhs._id),
      _name(rhs._name),
      _lb(rhs._lb),
      _ub(rhs._ub) {

}
Variable::Variable(size_t id)
    : _id(new size_t(id)),
      _name(new std::string),
      _lb(new double(-MyAmplInf)),
      _ub(new double(MyAmplInf)) {

}
std::string const & Variable::name() const {
  return *_name;
}
std::string & Variable::name() {
  return *_name;
}

size_t Variable::id() const {
  return *_id;
}
double Variable::lb() const {
  return *_lb;
}
double &Variable::lb() {
  return *_lb;
}
double Variable::ub() const {
  return *_ub;
}
double &Variable::ub() {
  return *_ub;
}
/**
 *
 */
Function::Function()
    : _qTerms(new QuadraticTerms),
      _lTerms(new LinearTerms),
      _cst(new double(0)) {

}

Function::Function(Function const & rhs)
    : _qTerms(rhs._qTerms),
      _lTerms(rhs._lTerms),
      _cst(rhs._cst) {
}
Function::~Function() {

}

QuadraticTerms const & Function::qTerms() const {
  return *_qTerms;
}
LinearTerms const & Function::lTerms() const {
  return *_lTerms;
}
double const & Function::cst() const {
  return *_cst;
}
QuadraticTerms & Function::qTerms() {
  return *_qTerms;
}
LinearTerms & Function::lTerms() {
  return *_lTerms;
}
double & Function::cst() {
  return *_cst;
}
void Function::add(Variable const & i, Variable const & j, double v) {
  if (v != 0)
    (*_qTerms)[QuadraticIndex(i, j)] += v;
}
void Function::add(QuadraticIndex ij, double v) {
  if (v != 0)
    (*_qTerms)[ij] += v;
}
void Function::add(Variable const & i, double v) {
  if (v != 0)
    (*_lTerms)[i] += v;
}
void Function::add(double v) {
  if (v != 0)
    (*_cst) += v;
}

void Function::operator+=(Variable const & rhs) {
  add(rhs, 1);
}
void Function::operator+=(Function const & rhs) {
  for (auto const & term : rhs.qTerms()) {
    add(term.first, term.second);
  }
  for (auto const & term : rhs.lTerms()) {
    add(term.first, term.second);
  }
  add(rhs.cst());
}
void Function::operator+=(double rhs) {
  add(rhs);
}

void Function::operator-=(Variable const & rhs) {
  add(rhs, -1);
}
void Function::operator-=(Function const & rhs) {
  for (auto const & term : rhs.lTerms()) {
    add(term.first, -term.second);
  }
  for (auto const & term : rhs.qTerms()) {
    add(term.first, -term.second);
  }
  add(-rhs.cst());
}
void Function::operator-=(double rhs) {
  add(-rhs);
}

void Function::operator*=(Variable const & rhs) {
  assert(isLinear());
  for (auto const & term : lTerms()) {
    add(term.first, rhs, term.second);
  }
  lTerms().clear();
  add(rhs, cst());
  cst() = 0;
}
void Function::operator*=(Function const & rhs) {
  assert(isLinear() && rhs.isLinear());
  Function result;
  for (auto const & term1 : lTerms()) {
    for (auto const & term2 : rhs.lTerms()) {
      result.add(term1.first, term2.first, term1.second * term2.second);
    }
  }
  for (auto const & term : lTerms()) {
    result.add(term.first, term.second * rhs.cst());
  }
  for (auto const & term : rhs.lTerms()) {
    result.add(term.first, term.second * cst());
  }
  result.cst() *= cst() * rhs.cst();

  cst() = result.cst();
  lTerms() = result.lTerms();
  qTerms() = result.qTerms();
}
void Function::operator*=(double rhs) {
  cst() *= rhs;
  for (auto & term : qTerms()) {
    term.second *= rhs;
  }
  for (auto & term : lTerms()) {
    term.second *= rhs;
  }
}
bool Function::isCst() const {
  return qTerms().empty() && lTerms().empty();
}
bool Function::isLinear() const {
  return qTerms().empty();
}

void Function::print(std::ostream & stream) const {
  if (cst() != 0)
    stream << cst();
  stream << qTerms();
  stream << lTerms();
}
void Function::clear() {
  *_cst = 0;
  lTerms().clear();
  qTerms().clear();
}

Constraint::Constraint()
    : Function(),
      _sense(new SENSE(EQ)) {

}
Constraint::Constraint(Constraint const & rhs)
    : Function(rhs),
      _sense(rhs._sense) {

}
Constraint::Constraint(Function const & f, SENSE sense, double rhs)
    : Function(f),
      _sense(new SENSE(sense)) {
  add(-rhs);
}
Constraint::Constraint(Function const & function, double lb, double ub)
    : Function(function),
      _sense(new SENSE(EQ)) {
//  std::cout << "lb : " << lb << std::endl;
//  std::cout << "ub : " << ub << std::endl;
  assert(lb <= -MyAmplInf || ub >= MyAmplInf || lb == ub);

  if (lb == ub) {
    add(-lb);
    sense() = EQ;
  } else if (lb <= -MyAmplInf) {
    add(-ub);
    sense() = LEQ;
  } else if (ub >= +MyAmplInf) {
    sense() = GEQ;
    add(-lb);
  }
}
Constraint::~Constraint() {

}
SENSE Constraint::sense() const {
  return *_sense;
}
SENSE & Constraint::sense() {
  return *_sense;
}
void Constraint::setSense(SENSE sense) {
  *_sense = sense;
}

void Constraint::print(std::ostream & stream) const {
  stream << "0" << sense();
  Function::print(stream);
}
/**
 *
 */

Function operator-(Variable const & rhs) {
  Function result;
  result -= rhs;
  return result;
}
Function operator-(Function const & rhs) {
  Function result;
  result -= rhs;
  return result;
}

Function operator+(Variable const & rhs) {
  Function result;
  result += rhs;
  return result;
}
Function operator+(Function const & rhs) {
  return rhs;
}
/**
 *
 */
Function operator+(Variable const & lhs, double rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}
Function operator+(Variable const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}
Function operator+(double lhs, Variable const &rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}

Function operator-(Variable const & lhs, double rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator-(Variable const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator-(double lhs, Variable const &rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator*(Variable const & lhs, double rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
Function operator*(Variable const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
Function operator*(double lhs, Variable const &rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
/**
 *
 */

Function operator+(Function const & lhs, double rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}
Function operator+(Function const & lhs, Function const & rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}
Function operator+(double lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}

Function operator-(Function const & lhs, double rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator-(Function const & lhs, Function const & rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator-(double lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator*(Function const & lhs, double rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
Function operator*(Function const & lhs, Function const & rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
Function operator*(double lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
/**
 *
 */

Function operator+(Function const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}
Function operator+(Variable const & lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result += rhs;
  return result;
}

Function operator-(Function const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator-(Variable const & lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result -= rhs;
  return result;
}
Function operator*(Function const & lhs, Variable const & rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
Function operator*(Variable const & lhs, Function const &rhs) {
  Function result;
  result += lhs;
  result *= rhs;
  return result;
}
/**
 *
 */
Constraint operator<=(Function const & function, double rhs) {
  std::cout << "function " << function << std::endl;
  return Constraint(function, LEQ, rhs);
}
Constraint operator==(Function const & function, double rhs) {
  return Constraint(function, EQ, rhs);
}
Constraint operator>=(Function const & function, double rhs) {
  return Constraint(function, GEQ, rhs);
}

Constraint operator<=(double lhs, Function const & function) {
  return Constraint(function, GEQ, lhs);
}
Constraint operator==(double lhs, Function const & function) {
  return Constraint(function, EQ, lhs);
}
Constraint operator>=(double lhs, Function const & function) {
  return Constraint(function, LEQ, lhs);
}
/**
 *
 */
std::ostream & operator<<(std::ostream & stream, SENSE sense) {
  switch (sense) {
    case GEQ:
      stream << ">=";
      break;
    case EQ:
      stream << "==";
      break;
    case LEQ:
      stream << "<=";
      break;
    default:
      assert(false && "UNKOWN ENUM CASE");
      break;
  }
  return stream;
}
std::ostream & operator<<(std::ostream & stream, Variable const &rhs) {
  stream << rhs.name();
  return stream;
}
std::ostream & operator<<(std::ostream & stream, QuadraticIndex const &rhs) {
  rhs.print(stream);
  return stream;
}
std::ostream & operator<<(std::ostream & stream, QuadraticTerms const &rhs) {

  for (auto & term : rhs) {
    stream << Formatter::Print(term.second);
    stream << term.first;
  }
  return stream;
}
std::ostream & operator<<(std::ostream & stream, LinearTerms const &rhs) {
  for (auto & term : rhs) {
    stream << Formatter::Print(term.second);
    stream << term.first;
  }
  return stream;
}
std::ostream & operator<<(std::ostream & stream, Function const &rhs) {
  rhs.print(stream);
  return stream;
}
std::ostream & operator<<(std::ostream & stream, Constraint const & rhs) {
  rhs.print(stream);
  return stream;
}
std::ostream & operator<<(std::ostream & stream, Problem const & rhs) {
  rhs.print(stream);
  return stream;
}

void operator<<(Problem & problem, Constraint const & constraint) {
  problem.newConstraint(constraint);

}
void operator>>(Constraint const & constraint, Problem& problem) {
  problem.newConstraint(constraint);

}
/**
 *
 */

Variable const & Problem::newVariable(std::string const & name) {
  _variables.push_back(Variable(_variables.size()));
  _variables.back().name() = name;
  return _variables.back();
}
Variable const & Problem::newVariable(double lb, double ub) {
  _variables.push_back(Variable(_variables.size()));
  _variables.back().name() = "x_" + Formatter::ToString(_variables.size());
  _variables.back().lb() = lb;
  _variables.back().ub() = ub;
  return _variables.back();
}

void QuadraticIndex::print(std::ostream &stream) const {
  stream << first.name();
  stream << ".";
  stream << second.name();
}

void Problem::print(std::ostream & stream) const {
  if (!_objectifs.empty())
    stream << "minimize " << _objectifs.front() << std::endl;
  else
    stream << "minimize 0" << std::endl;
  for (auto const& constraint : _constraints) {
    stream << constraint << std::endl;
  }

}
void Problem::newConstraint(Constraint const & constraint) {
  _constraints.push_back(constraint);
}
size_t Problem::addObj(Function const & function) {
  _objectifs.push_back(function);
  return _objectifs.size() - 1;
}

Variable & Problem::variable(size_t idx) {
  return _variables[idx];
}
Objectifs & Problem::objectifs() {
  return _objectifs;
}
Objectifs
const & Problem::objectifs() const {
  return _objectifs;
}
Variables & Problem::variables() {
  return _variables;
}
Variables
const & Problem::variables() const {
  return _variables;
}
void Problem::clear() {
  _constraints.clear();
  _objectifs.clear();
  _variables.clear();
}

void Problem::getQuadraticIndexes(QuadraticTerms & result) const {
  result.clear();
  for (auto const & f : _objectifs) {
    for (auto const & term : f.qTerms()) {
      result[term.first] += 1;
    }
  }
  for (auto const & f : _constraints) {
    for (auto const & term : f.qTerms()) {
      result[term.first] += 1;
    }
  }
  std::cout << "result : " << result << std::endl;
}
Constraints & Problem::constraints() {
  return _constraints;
}
Constraints
const & Problem::constraints() const {
  return _constraints;
}
