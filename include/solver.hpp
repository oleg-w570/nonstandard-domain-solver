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

  std::vector<std::vector<double>> u, v, diff;
  double max_diff, max_diff_x, max_diff_y;
  double initial_discrepancy;
  double result_discrepancy;
};
