#ifndef eigen_decomposition
#define eigen_decomposition

#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <omp.h>

class matrix{
    private:
        size_t rows_, cols_;
        std::vector<double> data;

    public:
        matrix(size_t n1, size_t n2);
        void set(const std::vector<double>& data_given);
        double& operator()(int i, int j);
        const double& operator()(int i, int j) const;
        size_t rows() const;
        size_t cols() const;
        std::vector<double>& get_data();
        matrix add(const matrix& other) const; 
        
        friend matrix operator+(const matrix& a, const matrix& b);
        friend matrix operator-(const matrix& a, const matrix& b);
        friend matrix operator*(const matrix& a, const matrix& b);
        friend matrix operator*(const matrix& mat, double scalar);
        friend matrix operator*(double scalar, const matrix& mat);

        friend std::ostream& operator<<(std::ostream& out, const matrix& mat);

        int mrow,mcol;
        double mval;
        void max_offdiag();
};

class jacobi_solver{
    private:
        const matrix* oA;
        matrix eigenvectors;
        std::vector<double> eigenvalues;

    public:
        jacobi_solver(const matrix& A);
        void solve(double tol = 1e-10, int max_iter = 1000);
        const matrix& get_eigenvectors() const;
        const std::vector<double>& get_eigenvalues() const;

        
};
std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec);

#endif // eigen_decomposition

