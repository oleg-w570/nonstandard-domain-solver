// Copyright Zorin Oleg
#pragma once
#include <vector>

#include "chebyshev_method.hpp"
#include "grid.hpp"
#include "task.hpp"

class Solver {
  Task task;
  Grid grid;
  ChebyshevMethod method;

  std::vector<std::vector<double>> exact_solution;
  std::vector<std::vector<double>> numerical_solution;
  std::vector<std::vector<double>> diff;
  std::vector<double> f_values;
  double max_diff, max_diff_x, max_diff_y;
  double initial_discrepancy;
  double result_discrepancy;

  void CalculateExactSolution();
  void InitializeNumericalSolution();
  void InitializeFValues();
  double CalculateDiscrepancy();
  void CalculateDiffSolutions();

 public:
  Solver(std::size_t n, std::size_t m, double eps, unsigned max_iter,
         unsigned K);
  void Solve();
  [[nodiscard]] const std::vector<std::vector<double>>& GetExactSolution() const;
  [[nodiscard]] const std::vector<std::vector<double>>& GetNumericalSolution() const;
  [[nodiscard]] const std::vector<std::vector<double>>& GetDiff() const;
  [[nodiscard]] double GetMaxDiff() const;
  [[nodiscard]] double GetMaxDiffX() const;
  [[nodiscard]] double GetMaxDiffY() const;
  [[nodiscard]] double GetInitialDiscrepancy() const;
  [[nodiscard]] double GetResultDiscrepancy() const;
  [[nodiscard]] double GetAccuracy() const;
  [[nodiscard]] unsigned GetIterationCount() const;
};
