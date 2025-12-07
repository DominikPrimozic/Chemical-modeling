#include <VectorOps.h>

// Vector subtraction: v1 - v2
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match.");
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) result[i] = a[i] - b[i];
    return result;
}

// Scalar multiplication: scalar * vector
std::vector<double> operator*(double scalar, const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) result[i] = scalar * v[i];
    return result;
}

// Vector addition: v1 + v2
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match.");
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) result[i] = a[i] + b[i];
    return result;
}

// In-place addition: v1 += v2
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match.");
    for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
    return a;
}

// In-place addition: v1 -= v2
std::vector<double>& operator-=(std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Vector sizes must match.");
    for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
    return a;
}

void set_zero(std::vector<std::vector<double>>& a){
    int inner_size=a[0].size();
    #pragma omp parallel for
    for (int i=0;i<a.size();i++){
        for (int j=0;j<inner_size;j++){
            a[i][j]=0;
        }
    }
}