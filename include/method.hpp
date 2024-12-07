#pragma once
#include <vector>

#include "grid.hpp"

class Method {
 protected:
  const Grid &grid_;

  double eps_;
  unsigned max_iter_;

  double accuracy_;
  unsigned n_iter_;

 public:
  Method(const Grid &grid, double eps, unsigned max_iter)
      : grid_(grid), eps_(eps), max_iter_(max_iter), accuracy_(eps + 1.0), n_iter_(0u) {}

  virtual void run(std::vector<std::vector<double>> &v, const std::vector<bool> &node_mask,
                   const std::vector<double> &f_values) = 0;

  double accuracy() const { return accuracy_; }
  double iteration_count() const { return n_iter_; }
};