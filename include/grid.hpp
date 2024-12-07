// Copyright Zorin Oleg
#pragma once
#include <vector>

struct Grid {
  std::size_t n, m;
  double h, k;
  double h2, k2, A;
  std::vector<double> x, y;

  Grid(std::size_t n, std::size_t m, double a, double b, double c, double d);
};
