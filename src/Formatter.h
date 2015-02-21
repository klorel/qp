/*
 * Formatter.h
 *
 *  Created on: 8 févr. 2015
 *      Author: mruiz
 */

#ifndef SRC_FORMATTER_H_
#define SRC_FORMATTER_H_
#include <string>
#include <sstream>
class Formatter {
 public:
  static std::string Print(double );
  template<class T> static std::string ToString(T value);
};

template<class T> inline std::string Formatter::ToString(T value){
  std::stringstream buffer;
  buffer << value;
  return buffer.str();
}
#endif /* SRC_FORMATTER_H_ */
