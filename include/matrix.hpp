// Copyright Zorin Oleg
#pragma once
#include <vector>

class Matrix {
  std::size_t n;
  std::size_t m;
  std::vector<double> data;

 public:
  explicit Matrix(std::size_t n = 0, std::size_t m = 0);
  double &operator()(std::size_t i, std::size_t j);
  const double &operator()(std::size_t i, std::size_t j) const;
  double &operator[](std::size_t idx);
  const double &operator[](std::size_t idx) const;
  std::size_t row(std::size_t idx) const;
  std::size_t col(std::size_t idx) const;
  std::size_t size() const;
  void clear();
};

double MaxDifference(const Matrix &v1, const Matrix &v2);
