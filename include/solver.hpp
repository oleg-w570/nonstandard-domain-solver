#pragma once
#include <vector>

#include "task.hpp"

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

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
  double VectorDiffNorm(const Matrix& v1, const Matrix& v2) const;
  void UpdateElement(std::size_t i, std::size_t j, double t, const Matrix& v_prev);
  void ChebishevIteration(int s, const Matrix& v_prev);
  void CopyMatrix(const Matrix& from, Matrix& to);
  void ChebishevMethod();
  void CalculateDiffSolutions();
  void UpdateDiscrepancy(std::size_t i, std::size_t j, double& R_max);
  double CalculateDiscrepancy();

 public:
  explicit Solver(int n = 8, int m = 8, double eps = 1e-6, int max_iter = 10000, int K = 4);
  void Solve();
};