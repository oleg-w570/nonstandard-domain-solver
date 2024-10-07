// Copyright Zorin Oleg
#include "chebyshev_method.hpp"

#include <cmath>
#include <numbers>

ChebyshevMethod::ChebyshevMethod(const Task &task, const Grid &grid, const double eps,
                                 const unsigned max_iter, const unsigned K)
    : task(task),
      grid(grid),
      eps(eps),
      max_iter(max_iter),
      K(K),
      tau(K),
      accuracy(eps + 1.0),
      n_iter(0) {
  InitializeChebyshevParameters();
}

void ChebyshevMethod::run(std::vector<std::vector<double>> &v) {
  const auto &skip_node = task.GetNodeValidPredicate(grid.n, grid.m);
  v_local_prev = v;
  v_global_prev = v;

  while (n_iter < max_iter && accuracy > eps) {
    v_global_prev = v;

    for (const auto &t : tau) {
      std::swap(v_local_prev, v);

      for (std::size_t i = 1; i < grid.n; ++i) {
        for (std::size_t j = 1; j < grid.m; ++j) {
          if (skip_node(i, j)) continue;
          v[i][j] = v_local_prev[i][j] * (1 + t * grid.A) +
                    t * (grid.k2 * (v_local_prev[i][j - 1] + v_local_prev[i][j + 1]) +
                         grid.h2 * (v_local_prev[i - 1][j] + v_local_prev[i + 1][j]) -
                         task.F(grid.x[i], grid.y[j]));
        }
      }
    }

    accuracy = MaxDifference(v, v_global_prev);
    n_iter += K;
  }
}

void ChebyshevMethod::InitializeChebyshevParameters() {
  constexpr double pi = std::numbers::pi;

  const double M_min = 4 * grid.h2 * std::sin(pi / static_cast<double>(2 * grid.n)) *
                           std::sin(pi / static_cast<double>(2 * grid.n)) +
                       4 * grid.k2 * std::sin(pi / static_cast<double>(2 * grid.m)) *
                           std::sin(pi / static_cast<double>(2 * grid.m));
  const double M_max =
      4 * grid.h2 *
          std::sin(pi * static_cast<double>(grid.n - 1) / static_cast<double>(2 * grid.n)) *
          std::sin(pi * static_cast<double>(grid.n - 1) / static_cast<double>(2 * grid.n)) +
      4 * grid.k2 *
          std::sin(pi * static_cast<double>(grid.m - 1) / static_cast<double>(2 * grid.m)) *
          std::sin(pi * static_cast<double>(grid.m - 1) / static_cast<double>(2 * grid.m));

  for (int s = 0; s < K; ++s) {
    tau[s] =
        1 / (0.5 * (M_max + M_min) + 0.5 * (M_max - M_min) * std::cos(pi * (1 + 2 * s) / (2 * K)));
  }
}

double ChebyshevMethod::MaxDifference(const std::vector<std::vector<double>> &v1,
                                      const std::vector<std::vector<double>> &v2) const {
  double max_diff = 0.0;

  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      const double curr_diff = std::abs(v1[i][j] - v2[i][j]);
      if (curr_diff > max_diff) max_diff = curr_diff;
    }
  }

  return max_diff;
}
