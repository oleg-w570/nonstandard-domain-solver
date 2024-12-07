// Copyright Zorin Oleg
#include "solver.hpp"

#include <memory>

#include "chebyshev_method.hpp"

Solver::Solver(std::unique_ptr<Task> task, MethodType method_type, const std::size_t n, const std::size_t m, double eps,
               unsigned max_iter)
    : task(std::move(task)),
      grid(n, m, task->a(), task->b(), task->c(), task->d()),
      exact_solution(n + 1, std::vector<double>(m + 1)),
      numerical_solution(n + 1, std::vector<double>(m + 1)),
      diff(n + 1, std::vector<double>(m + 1)),
      f_values((n + 1) * (m + 1)),
      skip_node_mask((n + 1) * (m + 1)) {
  const auto K = 16u;
  switch (method_type) {
    case MethodType::Chebyshev:
      method = std::make_unique<ChebyshevMethod>(grid, eps, max_iter, K);
      break;
    default:
      throw std::invalid_argument("Unknown method type");
      break;
  }
}

void Solver::InitializeNodeMask() {
  const auto skip_node = task->SkipNodePredicate(grid.n, grid.m);

  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      skip_node_mask[i * (grid.m + 1) + j] = skip_node(i, j);
    }
  }
}

void Solver::InitializeFValues() {
  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      f_values[i * (grid.m + 1) + j] = task->F(grid.x[i], grid.y[j]);
    }
  }
}

void Solver::ComputeExactSolution() {
  const auto is_border = task->BorderPredicate(grid.n, grid.m);
  const auto skip_node = task->SkipNodePredicate(grid.n, grid.m);

  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      if (is_border(i, j) || !skip_node_mask[i * (grid.m + 1) + j]) {
        exact_solution[i][j] = task->U(grid.x[i], grid.y[j]);
      }
    }
  }
}

void Solver::InitializeNumericalSolution() {
  const auto is_border = task->BorderPredicate(grid.n, grid.m);

  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      if (is_border(i, j)) {
        numerical_solution[i][j] = task->U(grid.x[i], grid.y[j]);
      }
    }
  }
}

double Solver::ComputeDiscrepancy() {
  auto max_discrepancy = 0.0;

  for (std::size_t i = 1; i < grid.n; ++i) {
    for (std::size_t j = 1; j < grid.m; ++j) {
      if (skip_node_mask[i * (grid.m + 1) + j]) continue;
      const auto discrepancy = std::abs(
          grid.A * numerical_solution[i][j] + grid.h2 * (numerical_solution[i - 1][j] + numerical_solution[i + 1][j]) +
          grid.k2 * (numerical_solution[i][j - 1] + numerical_solution[i][j + 1]) - f_values[i * (grid.m + 1) + j]);
      if (discrepancy > max_discrepancy) {
        max_discrepancy = discrepancy;
      }
    }
  }

  return max_discrepancy;
}

void Solver::ComputeDiffSolutions() {
  for (std::size_t i = 0; i < grid.n + 1; ++i) {
    for (std::size_t j = 0; j < grid.m + 1; ++j) {
      const auto curr_diff = std::abs(exact_solution[i][j] - numerical_solution[i][j]);
      diff[i][j] = curr_diff;
      if (curr_diff > max_diff) {
        max_diff = curr_diff;
        max_diff_x = grid.x[i];
        max_diff_y = grid.y[j];
      }
    }
  }
}

void Solver::Solve() {
  InitializeNodeMask();
  InitializeFValues();
  ComputeExactSolution();
  InitializeNumericalSolution();

  initial_discrepancy = ComputeDiscrepancy();

  method->run(numerical_solution, skip_node_mask, f_values);

  result_discrepancy = ComputeDiscrepancy();

  ComputeDiffSolutions();
}

const std::vector<std::vector<double>> &Solver::GetExactSolution() const { return exact_solution; }

const std::vector<std::vector<double>> &Solver::GetNumericalSolution() const { return numerical_solution; }

const std::vector<std::vector<double>> &Solver::GetDiff() const { return diff; }

double Solver::GetMaxDiff() const { return max_diff; }

double Solver::GetMaxDiffX() const { return max_diff_x; }

double Solver::GetMaxDiffY() const { return max_diff_y; }

double Solver::GetInitialDiscrepancy() const { return initial_discrepancy; }

double Solver::GetResultDiscrepancy() const { return result_discrepancy; }

double Solver::GetAccuracy() const { return method->accuracy(); }

double Solver::GetIterationCount() const { return method->iteration_count(); }
