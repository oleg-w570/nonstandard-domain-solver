// Copyright Zorin Oleg
#include "grid.hpp"

Grid::Grid(const std::size_t n, const std::size_t m, const double a, const double b, const double c,
           const double d)
    : n(n), m(m), x(n + 1), y(m + 1) {
  h = (b - a) / static_cast<double>(n);
  k = (d - c) / static_cast<double>(m);

  h2 = 1 / (h * h);
  k2 = 1 / (k * k);
  A = -2 * (h2 + k2);

  double xi = a;
  for (std::size_t i = 0; i < n + 1; ++i) {
    x[i] = xi;
    xi += h;
  }

  double yj = c;
  for (std::size_t j = 0; j < m + 1; ++j) {
    y[j] = yj;
    yj += k;
  }
}
