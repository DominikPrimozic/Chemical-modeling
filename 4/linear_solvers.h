#ifndef linear_solvers
#define linear_solvers

#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

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

        friend std::ostream& operator<<(std::ostream& out, const matrix& mat);
};

class csr{
    private:
        std::vector<double> value;
        std::vector<int> col;
        std::vector<int> row;

        void convert(const matrix& mat);

    public:
        csr(const matrix& mat);
        csr();

        void tridiagonal(int n, double sub, double diag, double up);
        void btridiagonal(int n, int n1, double sub, double diag, double up, double extra);
        void block_tridiagonal(const matrix& A, const matrix& B);
        void block_tridiagonal(const csr& A, const csr& B);
        std::vector<double> matvec(const std::vector<double>& vec) const; 

        void print_matrix(int numRows, int numCols) const;

        void set_value_at(size_t index, double val);
        void set_col_at(size_t index, int val);
        void set_row_at(size_t index, int val);
        double get_value_at(size_t index) const;
        int get_col_at(size_t index) const;
        int get_row_at(size_t index) const;
        void modify_last_diagonal(double new_value);
        double get_value_at(int row1, int col1) const;
        void modify_off_diagonal(int row1, int col1, double value_to_modify); 

};

std::vector<double> GC_solver(const csr& A,const std::vector<double>& b, double tol = 1e-6, int max_iter = 10000);

#endif // linear_solvers
