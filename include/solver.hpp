// Copyright Zorin Oleg
#pragma once
#include <memory>
#include <vector>

#include "grid.hpp"
#include "method.hpp"
#include "task.hpp"

enum class MethodType {
  Chebyshev,
};

class Solver {
  std::unique_ptr<Task> task;
  Grid grid;
  std::unique_ptr<Method> method;

  std::vector<std::vector<double>> exact_solution;
  std::vector<std::vector<double>> numerical_solution;
  std::vector<std::vector<double>> diff;
  std::vector<double> f_values;
  std::vector<bool> skip_node_mask;
  double max_diff, max_diff_x, max_diff_y;
  double initial_discrepancy;
  double result_discrepancy;

  void ComputeExactSolution();
  void InitializeNumericalSolution();
  void InitializeFValues();
  void InitializeNodeMask();
  double ComputeDiscrepancy();
  void ComputeDiffSolutions();

 public:
  Solver(std::unique_ptr<Task> task, MethodType method_type, std::size_t n, std::size_t m, double eps, unsigned max_iter);
  void Solve();
  [[nodiscard]] const std::vector<std::vector<double>> &GetExactSolution() const;
  [[nodiscard]] const std::vector<std::vector<double>> &GetNumericalSolution() const;
  [[nodiscard]] const std::vector<std::vector<double>> &GetDiff() const;
  [[nodiscard]] double GetMaxDiff() const;
  [[nodiscard]] double GetMaxDiffX() const;
  [[nodiscard]] double GetMaxDiffY() const;
  [[nodiscard]] double GetInitialDiscrepancy() const;
  [[nodiscard]] double GetResultDiscrepancy() const;
  [[nodiscard]] double GetAccuracy() const;
  [[nodiscard]] double GetIterationCount() const;
};
