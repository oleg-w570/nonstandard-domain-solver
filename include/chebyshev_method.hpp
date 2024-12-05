// Copyright Zorin Oleg
#pragma once
#include <vector>

#include "grid.hpp"

class ChebyshevMethod {
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
  [[nodiscard]] auto
  MaxDifference(const std::vector<std::vector<double>> &v1,
                const std::vector<std::vector<double>> &v2) const;

public:
  ChebyshevMethod(const Grid &grid, double eps, unsigned max_iter, unsigned K);
  void run(std::vector<std::vector<double>> &v,
           const std::vector<bool> &node_mask,
           const std::vector<double> &f_values);
  [[nodiscard]] double GetAccuracy() const;
  [[nodiscard]] double GetIterationCount() const;
};
