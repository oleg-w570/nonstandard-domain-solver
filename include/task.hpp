// Copyright Zorin Oleg
#pragma once
#include <cstddef>
#include <functional>

class Task {
 protected:
  double a_;
  double b_;
  double c_;
  double d_;

 public:
  virtual ~Task() = default;

  double a() const { return a_; }
  double b() const { return b_; }
  double c() const { return c_; }
  double d() const { return d_; }

  virtual double U(double x, double y) const = 0;
  virtual double F(double x, double y) const = 0;

  virtual std::function<bool(std::size_t, std::size_t)> BorderPredicate(std::size_t n, std::size_t m) const = 0;
  virtual std::function<bool(std::size_t, std::size_t)> SkipNodePredicate(std::size_t n, std::size_t m) const = 0;
};

class TestTask : public Task {
 public:
  TestTask();
  
  virtual double U(double x, double y) const override;
  virtual double F(double x, double y) const override;

  virtual std::function<bool(std::size_t, std::size_t)> BorderPredicate(std::size_t n, std::size_t m) const override;
  virtual std::function<bool(std::size_t, std::size_t)> SkipNodePredicate(std::size_t n, std::size_t m) const override;
};