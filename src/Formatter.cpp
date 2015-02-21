/*
 * Formatter.cpp
 *
 *  Created on: 8 févr. 2015
 *      Author: mruiz
 */

#include "Formatter.h"

std::string Formatter::Print(double value) {
  if (value == 1)
    return "+";
  else if (value == -1)
    return "-";
  else if (value > 0)
    return "+" + ToString(value);
  else
    return ToString(value);
}
