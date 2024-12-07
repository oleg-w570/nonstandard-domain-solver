// Copyright Zorin Oleg
#pragma once
#include <functional>

class Task {
 public:
  double a, b, c, d;

  Task();
  [[nodiscard]] static double U(double x, double y);
  [[nodiscard]] static double F(double x, double y);
  [[nodiscard]] static std::function<bool(std::size_t, std::size_t)> GetBorderPredicate(
      std::size_t n, std::size_t m);
  [[nodiscard]] static std::function<bool(std::size_t, std::size_t)> GetSkipNodePredicate(
      std::size_t n, std::size_t m);
};