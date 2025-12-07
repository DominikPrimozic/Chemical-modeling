#include <eigensolver.h>

matrix::matrix(size_t n1, size_t n2): rows_(n1),cols_(n2){
    data.resize(n1*n2,0);

}

void matrix::set(const std::vector<double>& data_given){
    data=data_given;
}

size_t matrix::rows() const { return rows_; }
size_t matrix::cols() const { return cols_; }

double& matrix::operator()(int i, int j){
    return data[i + rows_ * j];
}


const double& matrix::operator()(int i, int j) const{
    return data[i + rows_ * j];
}

std::ostream& operator<<(std::ostream& out, const matrix& mat) {
    for (size_t i = 0; i < mat.rows_; ++i) {
        for (size_t j = 0; j < mat.cols_; ++j) {
            out << mat(i, j);
            if (j + 1 < mat.cols_) {
                out << "\t"; 
            }
        }
        out << "\n";
    }
    return out;
}

std::vector<double>& matrix::get_data() {
    return data;
}


matrix matrix::add(const matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    }

    matrix result(rows_, cols_);

    size_t total_size = rows_ * cols_;

    #pragma omp simd
    for (size_t i = 0; i < total_size; ++i) {
        result.data[i] = data[i] + other.data[i];
    }

    return result;
}

matrix operator+(const matrix& a, const matrix& b) {
    return a.add(b);  
}

matrix operator-(const matrix& a, const matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction.");
    }

    matrix result(a.rows(), a.cols());
    size_t total_size = a.rows() * a.cols();

    #pragma omp simd
    for (size_t i = 0; i < total_size; ++i) {
        result.get_data()[i] = a.data[i] - b.data[i];
    }

    return result;
}


matrix operator*(const matrix& a, const matrix& b) {
    if (a.cols() != b.rows()) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication.");
    }

    matrix result(a.rows(), b.cols());
    #pragma omp parallel for
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            double sum = 0.0;
            for (int k = 0; k < a.cols(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }

    return result;
}

matrix operator*(const matrix& mat, double scalar) {
    matrix result(mat.rows(), mat.cols());
    size_t total_size = mat.rows() * mat.cols();

    #pragma omp simd
    for (size_t i = 0; i < total_size; ++i) {
        result.get_data()[i] = mat.data[i] * scalar;
    }

    return result;
}

matrix operator*(double scalar, const matrix& mat) {
    return mat * scalar;  
}


void matrix::max_offdiag() {
    mval = 0.0;
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = i + 1; j < cols_; ++j) {
            double val = std::abs((*this)(i, j));
            if (val > mval) {
                mval = val;
                mrow = i;
                mcol = j;
            }
        }
    }
}

jacobi_solver::jacobi_solver(const matrix& A): eigenvectors(A.rows(),A.cols()), oA(&A){
}

void jacobi_solver::solve( double tol, int max_iter){
    matrix A = *oA;
    size_t n = A.rows();
    //for (size_t i = 0; i < n; ++i) eigenvectors(i, i) = 1.0;
    for (size_t i = 0; i < n; ++i){ //better because it resets them
        for (size_t j = 0; j < n; ++j){
            eigenvectors(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }    
    for (int iter = 0; iter < max_iter; ++iter){
        A.max_offdiag();
        int p=A.mrow;
        int q=A.mcol;
        if (A.mval < tol) break;
        double app = A(p, p);
        double aqq = A(q, q);
        double apq = A(p, q);
    
        // Compute rotation
        double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
        double c = std::cos(phi);
        double s = std::sin(phi);
    
        // Rotate matrix A
        for (size_t i = 0; i < n; ++i) {
            if (i != p && i != q) {
                double aip = A(i, p);
                double aiq = A(i, q);
                A(i, p) = A(p, i) = c * aip - s * aiq;
                A(i, q) = A(q, i) = s * aip + c * aiq;
            }
        }
    
        // Update diagonal elements using saved values
        A(p, p) = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        A(q, q) = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A(p, q) = A(q, p) = 0.0;  // Explicitly zero out the off-diagonal element
    
        // Rotate eigenvectors
        for (size_t i = 0; i < n; ++i) {
            double vip = eigenvectors(i, p);
            double viq = eigenvectors(i, q);
            eigenvectors(i, p) = c * vip - s * viq;
            eigenvectors(i, q) = s * vip + c * viq;
        }
    }
    eigenvalues.resize(n);
    for (size_t i = 0; i < n; ++i)
        eigenvalues[i] = A(i, i);
}

const matrix& jacobi_solver::get_eigenvectors() const{
    return eigenvectors;
}

const std::vector<double>& jacobi_solver::get_eigenvalues() const{
    return eigenvalues;
}


std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        out << vec[i];
        if (i + 1 < vec.size()) {
            out << "\n";  // You can adjust this if you want something else like ", "
        }
    }
    return out;
}