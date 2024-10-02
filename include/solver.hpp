// Copyright Zorin Oleg
#pragma once
#include <vector>

#include "matrix.hpp"
#include "task.hpp"

using Vector = std::vector<double>;

class Solver {
 private:
  Task task;
  double eps;
  int max_iter;

  std::size_t n, m;
  std::size_t bottom, top, right, in_left, in_right;
  double h, k;
  double A, h2, k2;

  int K;
  Vector tau;
  Matrix v_local_prev;
  Matrix v_global_prev;

 public:
  Matrix u, v, diff;
  Vector x, y;
  double max_diff, max_diff_x, max_diff_y;
  double R_null, R_res;
  double accuracy;
  int n_iter;

 private:
  void SetUpGrid();
  void SetUpChebishevParameters();
  void CalculateBorder(Matrix& z);
  void CalculateTrueSolution();
  double ComputeNextValue(std::size_t i, std::size_t j, double t) const;
  void ChebishevLocalIteration(double t);
  void ChebishevGlobalIteration();
  void ChebishevMethod();
  void CalculateDiffSolutions();
  double ComputeDiscrepancyValue(std::size_t i, std::size_t j) const;
  double CalculateDiscrepancy();

 public:
  explicit Solver(int n = 8, int m = 8, double eps = 1e-6, int max_iter = 10000, int K = 4);
  void Solve();
};