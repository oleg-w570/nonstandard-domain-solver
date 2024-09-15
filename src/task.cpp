#include "task.hpp"

#include <cmath>

Task::Task() {
  a = -1;
  b = 1;
  c = -1;
  d = 1;
}

double Task::U(double x, double y) const { return std::exp(1 - x * x - y * y); }

double Task::F(double x, double y) const {
  return -4 * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
}