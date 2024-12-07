// Copyright Zorin Oleg
#include "chebyshev_method.hpp"

#include <cmath>

ChebyshevMethod::ChebyshevMethod(const Grid &grid_, const double eps, const unsigned max_iter, const unsigned K)
    : Method(grid_, eps, max_iter), K_(K), tau_(K) {
  InitializeChebyshevParameters();
}

void ChebyshevMethod::run(std::vector<std::vector<double>> &v, const std::vector<bool> &node_mask,
                          const std::vector<double> &f_values) {
  v_local_prev_ = v;
  v_global_prev_ = v;

#pragma omp parallel
  {
    while (n_iter_ < max_iter_ && accuracy_ > eps_) {
#pragma omp single
      { v_global_prev_ = v; }

      for (const auto t : tau_) {
#pragma omp single
        { std::swap(v_local_prev_, v); }

        const auto expr = 1.0 + t * grid_.A;
#pragma omp for collapse(2) schedule(auto)
        for (std::size_t i = 1; i < grid_.n; ++i) {
          for (std::size_t j = 1; j < grid_.m; ++j) {
            if (node_mask[i * (grid_.m + 1) + j]) continue;
            v[i][j] = v_local_prev_[i][j] * expr + t * (grid_.k2 * (v_local_prev_[i][j - 1] + v_local_prev_[i][j + 1]) +
                                                        grid_.h2 * (v_local_prev_[i - 1][j] + v_local_prev_[i + 1][j]) -
                                                        f_values[i * (grid_.m + 1) + j]);
          }
        }
      }

#pragma omp single
      {
        accuracy_ = 0.0;
        n_iter_ += K_;
      }

#pragma omp for reduction(max : accuracy_) collapse(2)
      for (std::size_t i = 1; i < grid_.n; ++i) {
        for (std::size_t j = 1; j < grid_.m; ++j) {
          const auto curr_diff = std::abs(v[i][j] - v_global_prev_[i][j]);
          if (curr_diff > accuracy_) accuracy_ = curr_diff;
        }
      }
    }
  }
}

void ChebyshevMethod::InitializeChebyshevParameters() {
  constexpr auto pi = 3.14159265358979323846;
  const auto M_min =
      4 * grid_.h2 * std::sin(pi / static_cast<double>(2 * grid_.n)) * std::sin(pi / static_cast<double>(2 * grid_.n)) +
      4 * grid_.k2 * std::sin(pi / static_cast<double>(2 * grid_.m)) * std::sin(pi / static_cast<double>(2 * grid_.m));
  const auto M_max = 4 * grid_.h2 * std::sin(pi * static_cast<double>(grid_.n - 1) / static_cast<double>(2 * grid_.n)) *
                         std::sin(pi * static_cast<double>(grid_.n - 1) / static_cast<double>(2 * grid_.n)) +
                     4 * grid_.k2 * std::sin(pi * static_cast<double>(grid_.m - 1) / static_cast<double>(2 * grid_.m)) *
                         std::sin(pi * static_cast<double>(grid_.m - 1) / static_cast<double>(2 * grid_.m));

  for (unsigned s = 0; s < K_; ++s) {
    tau_[s] = 1 / (0.5 * (M_max + M_min) + 0.5 * (M_max - M_min) * std::cos(pi * (1 + 2 * s) / (2 * K_)));
  }
}
