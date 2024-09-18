#include "solver.hpp"

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Solver::Solver(int n, int m, double eps, int max_iter, int K)
    : n(n), m(m), eps(eps), max_iter(max_iter), K(K) {
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
    z[n][j] = task.U(task.b, y[j]);
  }
  for (j = bottom; j <= top; ++j) {
    z[right][j] = task.U(x[right], y[j]);
  }
  for (j = top; j < m + 1; ++j) {
    z[n][j] = task.U(task.b, y[j]);
  }
  for (i = right; i < n + 1; ++i) {
    z[i][bottom] = task.U(x[i], y[bottom]);
    z[i][top] = task.U(x[i], y[top]);
  }

  for (j = bottom; j <= top; ++j) {
    z[in_left][j] = task.U(x[in_left], y[j]);
    z[in_right][j] = task.U(x[in_right], y[j]);
  }
  for (i = in_left; i <= in_right; ++i) {
    z[i][bottom] = task.U(x[i], y[bottom]);
    z[i][top] = task.U(x[i], y[top]);
  }

  for (i = 0; i < n + 1; ++i) {
    z[i][0] = task.U(x[i], task.c);
    z[i][m] = task.U(x[i], task.d);
  }
  for (j = 0; j < m + 1; ++j) {
    z[0][j] = task.U(task.a, y[j]);
  }
}

void Solver::CalculateTrueSolution() {
  CalculateBorder(u);
  std::size_t i, j;

  for (j = 1; j < bottom; ++j) {
    for (i = 1; i < n; ++i) {
      u[i][j] = task.U(x[i], y[j]);
    }
  }
  for (j = bottom; j < top + 1; ++j) {
    for (i = 1; i < in_left; ++i) {
      u[i][j] = task.U(x[i], y[j]);
    }
    for (i = in_right + 1; i < right; ++i) {
      u[i][j] = task.U(x[i], y[j]);
    }
  }
  for (j = top + 1; j < m; ++j) {
    for (i = 1; i < n; ++i) {
      u[i][j] = task.U(x[i], y[j]);
    }
  }
}

double Solver::VectorDiffNorm(const Matrix &v1, const Matrix &v2) const {
  double max_diff = 0;
  for (std::size_t i = 0; i < n + 1; ++i) {
    for (std::size_t j = 0; j < m + 1; ++j) {
      const double curr_diff = std::abs(v1[i][j] - v2[i][j]);
      if (curr_diff > max_diff) max_diff = curr_diff;
    }
  }
  return max_diff;
}

inline double Solver::ComputeNextValue(std::size_t i, std::size_t j, double t) const {
  return (v_local_prev[i][j] * (1 + t * A) +
          t * (k2 * (v_local_prev[i][j - 1] + v_local_prev[i][j + 1]) +
               h2 * (v_local_prev[i - 1][j] + v_local_prev[i + 1][j]) - task.F(x[i], y[j])));
}

void Solver::ChebishevLocalIteration(double t) {
  #pragma omp parallel
  {
    #pragma omp sections
    {
      #pragma omp section
      {
        // Параллелизация внутри секции
        #pragma omp parallel for schedule(dynamic)
        for (int i = 1; i < in_left; ++i) {
          for (std::size_t j = 1; j < m; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
        }
      }

      #pragma omp section
      {
        // Параллелизация внутри секции
        #pragma omp parallel for schedule(dynamic)
        for (int i = in_left; i <= in_right; ++i) {
          for (std::size_t j = 1; j < bottom; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
          for (std::size_t j = top + 1; j < m; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
        }
      }

      #pragma omp section
      {
        // Параллелизация внутри секции
        #pragma omp parallel for schedule(dynamic)
        for (int i = in_right + 1; i < right; ++i) {
          for (std::size_t j = 1; j < m; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
        }
      }

      #pragma omp section
      {
        // Параллелизация внутри секции
        #pragma omp parallel for schedule(dynamic)
        for (int i = right; i < n; ++i) {
          for (std::size_t j = 1; j < bottom; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
          for (std::size_t j = top + 1; j < m; ++j) {
            v[i][j] = ComputeNextValue(i, j, t);
          }
        }
      }
    }  // end of sections
  }  // end of parallel
}


inline void Solver::CopyMatrix(const Matrix &from, Matrix &to) {
  for (std::size_t i = 0; i < to.size(); ++i) {
    std::copy(from[i].begin(), from[i].end(), to[i].begin());
  }
}

void Solver::ChebishevGlobalIteration() {
  CopyMatrix(v, v_global_prev);

  for (const auto &t : tau) {
    std::swap(v, v_local_prev);
    ChebishevLocalIteration(t);
  }
}

void Solver::ChebishevMethod() {
  v_global_prev = Matrix(n + 1, Vector(m + 1));
  v_local_prev = Matrix(n + 1, Vector(m + 1));
  CalculateBorder(v_local_prev);

  n_iter = 0;
  accuracy = eps + 1;

  while (n_iter < max_iter && accuracy > eps) {
    ChebishevGlobalIteration();
    accuracy = VectorDiffNorm(v, v_global_prev);
    n_iter += K;
  }

  v_local_prev.clear();
  v_global_prev.clear();
}

void Solver::CalculateDiffSolutions() {
  max_diff = 0;

  for (std::size_t i = 0; i < n + 1; ++i) {
    for (std::size_t j = 0; j < m + 1; ++j) {
      const double curr_diff = std::abs(u[i][j] - v[i][j]);
      diff[i][j] = curr_diff;
      if (curr_diff > max_diff) {
        max_diff = curr_diff;
        max_diff_x = x[i];
        max_diff_y = y[j];
      }
    }
  }
}

inline double Solver::ComputeDiscrepancyValue(std::size_t i, std::size_t j) const {
  return std::abs(A * v[i][j] + h2 * (v[i - 1][j] + v[i + 1][j]) +
                  k2 * (v[i][j - 1] + v[i][j + 1]) - task.F(x[i], y[j]));
}

double Solver::CalculateDiscrepancy() {
  std::size_t i, j;
  double R_max = 0;

  for (j = 1; j < bottom; ++j) {
    for (i = 1; i < n; ++i) {
      const double R = ComputeDiscrepancyValue(i, j);
      if (R > R_max) R_max = R;
    }
  }
  for (j = bottom; j < top + 1; ++j) {
    for (i = 1; i < in_left; ++i) {
      const double R = ComputeDiscrepancyValue(i, j);
      if (R > R_max) R_max = R;
    }
    for (i = in_right + 1; i < right; ++i) {
      const double R = ComputeDiscrepancyValue(i, j);
      if (R > R_max) R_max = R;
    }
  }
  for (j = top + 1; j < m; ++j) {
    for (i = 1; i < n; ++i) {
      const double R = ComputeDiscrepancyValue(i, j);
      if (R > R_max) R_max = R;
    }
  }
  return R_max;
}

void Solver::Solve() {
  u = Matrix(n + 1, Vector(m + 1));
  CalculateTrueSolution();

  v = Matrix(n + 1, Vector(m + 1));
  CalculateBorder(v);
  R_null = CalculateDiscrepancy();

  ChebishevMethod();
  R_res = CalculateDiscrepancy();

  diff = Matrix(n + 1, Vector(m + 1));
  CalculateDiffSolutions();
}
