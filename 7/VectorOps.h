#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <vector>
#include <stdexcept>
#include <omp.h>

// Overload subtraction: v1 - v2
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

// Scalar multiplication: scalar * vector
std::vector<double> operator*(double scalar, const std::vector<double>& v);

// Vector addition: v1 + v2
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

// In-place addition: v1 += v2
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b);

// In-place addition: v1 -= v2
std::vector<double>& operator-=(std::vector<double>& a, const std::vector<double>& b);

void set_zero(std::vector<std::vector<double>>& a);

#endif // VECTOR_OPS_H