// Copyright Zorin Oleg
#include "solver.hpp"

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Solver::Solver(int n, int m, double eps, int max_iter, int K)
    : n(n),
      m(m),
      eps(eps),
      max_iter(max_iter),
      K(K),
      u(n + 1, m + 1),
      v(n + 1, m + 1),
      diff(n + 1, m + 1) {
  SetUpGrid();
  SetUpChebishevParameters();
}

void Solver::SetUpGrid() {
  h = (task.b - task.a) / n;
  k = (task.d - task.c) / m;

  h2 = 1 / (h * h);
  k2 = 1 / (k * k);
  A = -2 * (h2 + k2);

  bottom = m / 4;
  top = 3 * m / 4;
  right = 7 * n / 8;
  in_left = n / 4;
  in_right = n / 2;

  x = Vector(n + 1);
  y = Vector(m + 1);

  double xi = task.a;
  for (std::size_t i = 0; i < n + 1; ++i) {
    x[i] = xi;
    xi += h;
  }

  double yj = task.c;
  for (std::size_t j = 0; j < m + 1; ++j) {
    y[j] = yj;
    yj += k;
  }
}

void Solver::SetUpChebishevParameters() {
  tau = Vector(K);

  const double Mmin = 4 * h2 * std::sin(M_PI / (2 * n)) * std::sin(M_PI / (2 * n)) +
                      4 * k2 * std::sin(M_PI / (2 * m)) * std::sin(M_PI / (2 * m));
  const double Mmax =
      4 * h2 * std::sin(M_PI * (n - 1) / (2 * n)) * std::sin(M_PI * (n - 1) / (2 * n)) +
      4 * k2 * std::sin(M_PI * (m - 1) / (2 * m)) * std::sin(M_PI * (m - 1) / (2 * m));

  for (int s = 0; s < K; ++s) {
    tau[s] =
        1 / (0.5 * (Mmax + Mmin) + 0.5 * (Mmax - Mmin) * std::cos(M_PI * (1 + 2 * s) / (2 * K)));
  }
}

void Solver::CalculateBorder(Matrix &z) {
  std::size_t i, j;

  for (j = 0; j <= bottom; ++j) {
    z(n, j) = task.U(task.b, y[j]);
  }
  for (j = bottom; j <= top; ++j) {
    z(right, j) = task.U(x[right], y[j]);
  }
  for (j = top; j < m + 1; ++j) {
    z(n, j) = task.U(task.b, y[j]);
  }
  for (i = right; i < n + 1; ++i) {
    z(i, bottom) = task.U(x[i], y[bottom]);
    z(i, top) = task.U(x[i], y[top]);
  }

  for (j = bottom; j <= top; ++j) {
    z(in_left, j) = task.U(x[in_left], y[j]);
    z(in_right, j) = task.U(x[in_right], y[j]);
  }
  for (i = in_left; i <= in_right; ++i) {
    z(i, bottom) = task.U(x[i], y[bottom]);
    z(i, top) = task.U(x[i], y[top]);
  }

  for (i = 0; i < n + 1; ++i) {
    z(i, 0) = task.U(x[i], task.c);
    z(i, m) = task.U(x[i], task.d);
  }
  for (j = 0; j < m + 1; ++j) {
    z(0, j) = task.U(task.a, y[j]);
  }
}

void Solver::CalculateTrueSolution() {
  CalculateBorder(u);
  for (std::size_t idx = 0; idx < u.size(); ++idx) {
    std::size_t i = u.row(idx);
    std::size_t j = u.col(idx);
    if (i == 0 || j == 0 || i == n || j == m ||
        (i >= in_left && i <= in_right || i >= right) && j >= bottom && j <= top)
      continue;
    u[idx] = task.U(x[i], y[j]);
  }
}

inline double Solver::ComputeNextValue(std::size_t i, std::size_t j, double t) const {
  return (v_local_prev(i, j) * (1.0 + t * A) +
          t * (k2 * (v_local_prev(i, j - 1) + v_local_prev(i, j + 1)) +
               h2 * (v_local_prev(i - 1, j) + v_local_prev(i + 1, j)) - task.F(x[i], y[j])));
}

void Solver::ChebishevLocalIteration(double t) {
  #pragma omp parallel for
  for (int idx = 0; idx < v.size(); ++idx) {
    std::size_t i = v.row(idx);
    std::size_t j = v.col(idx);
    if (i == 0 || j == 0 || i == n || j == m ||
        (i >= in_left && i <= in_right || i >= right) && j >= bottom && j <= top)
      continue;
    v[idx] = ComputeNextValue(i, j, t);
  }
}

void Solver::ChebishevGlobalIteration() {
  v_global_prev = v;

  for (const auto &t : tau) {
    std::swap(v, v_local_prev);
    ChebishevLocalIteration(t);
  }
}

void Solver::ChebishevMethod() {
  v_global_prev = Matrix(n + 1, m + 1);
  v_local_prev = Matrix(n + 1, m + 1);
  CalculateBorder(v_local_prev);

  n_iter = 0;
  accuracy = eps + 1;

  while (n_iter < max_iter && accuracy > eps) {
    ChebishevGlobalIteration();
    accuracy = MaxDifference(v, v_global_prev);
    n_iter += K;
  }

  v_local_prev.clear();
  v_global_prev.clear();
}

void Solver::CalculateDiffSolutions() {
  max_diff = 0;

  for (std::size_t idx = 0; idx < u.size(); ++idx) {
    const double curr_diff = std::abs(u[idx] - v[idx]);
    diff[idx] = curr_diff;
    if (curr_diff > max_diff) {
      max_diff = curr_diff;
      max_diff_x = x[u.row(idx)];
      max_diff_y = y[u.col(idx)];
    }
  }
}

inline double Solver::ComputeDiscrepancyValue(std::size_t i, std::size_t j) const {
  return std::abs(A * v(i, j) + h2 * (v(i - 1, j) + v(i + 1, j)) +
                  k2 * (v(i, j - 1) + v(i, j + 1)) - task.F(x[i], y[j]));
}

double Solver::CalculateDiscrepancy() {
  double R_max = 0;

  for (std::size_t idx = 0; idx < v.size(); ++idx) {
    std::size_t i = v.row(idx);
    std::size_t j = v.col(idx);
    if (i == 0 || j == 0 || i == n || j == m ||
        (i >= in_left && i <= in_right || i >= right) && j >= bottom && j <= top)
      continue;
    const double R = ComputeDiscrepancyValue(i, j);
    if (R > R_max) R_max = R;
  }

  return R_max;
}

void Solver::Solve() {
  CalculateTrueSolution();

  CalculateBorder(v);
  R_null = CalculateDiscrepancy();

  ChebishevMethod();
  R_res = CalculateDiscrepancy();

  CalculateDiffSolutions();
}
