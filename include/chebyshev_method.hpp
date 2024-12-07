// Copyright Zorin Oleg
#pragma once
#include <vector>

#include "method.hpp"

class ChebyshevMethod : public Method {
  unsigned K_;
  std::vector<double> tau_;
  std::vector<std::vector<double>> v_local_prev_;
  std::vector<std::vector<double>> v_global_prev_;

  void InitializeChebyshevParameters();

 public:
  ChebyshevMethod(const Grid &grid, double eps, unsigned max_iter, unsigned K);

  virtual void run(std::vector<std::vector<double>> &v, const std::vector<bool> &node_mask,
                   const std::vector<double> &f_values) override;
};
