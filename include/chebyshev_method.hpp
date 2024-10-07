// Copyright Zorin Oleg
#pragma once
#include <vector>

#include "grid.hpp"
#include "task.hpp"

class ChebyshevMethod {
  const Task &task;
  const Grid &grid;

  double eps;
  unsigned max_iter;

  unsigned K;
  std::vector<double> tau;
  std::vector<std::vector<double>> v_local_prev;
  std::vector<std::vector<double>> v_global_prev;

  double accuracy;
  unsigned n_iter;

  void InitializeChebyshevParameters();
  [[nodiscard]] double MaxDifference(const std::vector<std::vector<double>> &v1,
                                     const std::vector<std::vector<double>> &v2) const;

 public:
  ChebyshevMethod(const Task &task, const Grid &grid, double eps, unsigned max_iter, unsigned K);
  void run(std::vector<std::vector<double>> &v);
};
