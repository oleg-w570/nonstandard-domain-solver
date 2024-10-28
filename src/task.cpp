// Copyright Zorin Oleg
#include "task.hpp"

#include <cmath>

Task::Task() {
  a = -1;
  b = 1;
  c = -1;
  d = 1;
}

double Task::U(const double x, const double y) {
  return std::exp(1 - x * x - y * y);
}

double Task::F(const double x, const double y) {
  return -4 * (1 - x * x - y * y) * std::exp(1 - x * x - y * y);
}

std::function<bool(std::size_t, std::size_t)> Task::GetBorderPredicate(
    std::size_t n, std::size_t m) {
  const std::size_t bottom = m / 4;
  const std::size_t top = 3 * m / 4;
  const std::size_t right = 7 * n / 8;
  const std::size_t in_left = n / 4;
  const std::size_t in_right = n / 2;

  return [=](const std::size_t i, const std::size_t j) {
    return (i == 0 || j == 0 || j == m ||
            (j >= bottom && j <= top &&
             (i == in_right || i == in_left || i == right)) ||
            (i == n && (j <= bottom || j >= top)) ||
            ((j == bottom || j == top) &&
             (i >= in_left && i <= in_right || i >= right)));
  };
}

std::function<bool(std::size_t, std::size_t)> Task::GetSkipNodePredicate(
    const std::size_t n, const std::size_t m) {
  const std::size_t bottom = m / 4;
  const std::size_t top = 3 * m / 4;
  const std::size_t right = 7 * n / 8;
  const std::size_t in_left = n / 4;
  const std::size_t in_right = n / 2;

  return [=](const std::size_t i, const std::size_t j) {
    return (((i >= in_left && i <= in_right) || i >= right) && j >= bottom &&
            j <= top);
  };
}