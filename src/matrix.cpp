// Copyright Zorin Oleg
#include "matrix.hpp"

Matrix::Matrix(std::size_t n, std::size_t m) : n(n), m(m), data(n * m) {}

double& Matrix::operator()(std::size_t i, std::size_t j) { return data[i * m + j]; }

const double& Matrix::operator()(std::size_t i, std::size_t j) const { return data[i * m + j]; }

double& Matrix::operator[](std::size_t idx) { return data[idx]; }

const double& Matrix::operator[](std::size_t idx) const { return data[idx]; }

std::size_t Matrix::row(std::size_t idx) const { return idx / m; }

std::size_t Matrix::col(std::size_t idx) const { return idx % m; }

std::size_t Matrix::size() const { return n * m; }

void Matrix::clear() { data.clear(); }

double MaxDifference(const Matrix& v1, const Matrix& v2) {
  double max_diff = 0;
  for (std::size_t idx = 0; idx < v1.size(); ++idx) {
    const double curr_diff = std::abs(v1[idx] - v2[idx]);
    if (curr_diff > max_diff) max_diff = curr_diff;
  }
  return max_diff;
}
