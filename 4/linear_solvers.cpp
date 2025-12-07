#include <linear_solvers.h>

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

csr::csr(const matrix& mat){
    convert(mat);
}
csr::csr(){
    
}


void csr::convert(const matrix& mat) {
    int rows = mat.rows();
    int cols = mat.cols();
    
    row.assign(cols + 1, 0);  
    value.clear();
    col.clear();

    for (int j = 0; j < cols; ++j) {  
        for (int i = 0; i < rows; ++i) {
            if (mat(i, j) != 0) {
                value.push_back(mat(i, j));
                col.push_back(i);  
            }
        }
        row[j + 1] = value.size(); 
    }
}

void csr::tridiagonal(int n, double sub, double diag, double up){
    row.push_back(0);
        for (int i = 0; i < n; ++i) {
            value.push_back(diag);
            col.push_back(i);

            // Sub-diagonal
            if (i < n - 1) {
                value.push_back(sub);
                col.push_back(i + 1);
            }

            // Sup-diagonal
            if (i > 0) {
                value.push_back(up);
                col.push_back(i - 1);
            }

            row.push_back(value.size());
        }
}

void csr::btridiagonal(int n, int n1, double sub, double diag, double up, double extra){
    row.push_back(0);
        for (int i = 0; i < n; ++i) {
            value.push_back(diag);
            col.push_back(i);

            // Sub-diagonal
            if (i < n - 1) {
                value.push_back(sub);
                col.push_back(i + 1);
            }

            // Sup-diagonal
            if (i > 0) {
                value.push_back(up);
                col.push_back(i - 1);
            }
            if (i + n1 < n) {
                value.push_back(extra);
                col.push_back(i + n1);
            }
    
            // Extra diagonal (upper every nth super-diagonal)
            if (i - n1 >= 0) {
                value.push_back(extra);
                col.push_back(i - n1);
            }

            row.push_back(value.size());
        }
}

void csr::block_tridiagonal(const matrix& A, const matrix& B){
    int block_size = A.rows(); // A is square
    int block_count = B.cols(); // B is also square
    int total_size = block_count * block_size; 

    row.clear();
    col.clear();
    value.clear();
    row.push_back(0); 

    for (int i = 0; i < block_count; ++i) {
        
        int block_start = i * block_size; 

        // A on diagonal
        for (int r = 0; r < block_size; ++r) {
            for (int c = 0; c < block_size; ++c) {
                if (A(r,c) != 0) {
                    value.push_back(A(r,c));
                    col.push_back(block_start + c);
                }
            }

            // Sub-diagonal 
            if (i > 0) {
                for (int c = 0; c < block_size; ++c) {
                    if (B(r,c) != 0) {
                        value.push_back(B(r,c));
                        col.push_back((i - 1) * block_size + c);
                    }
                }
            }

            // Super-diagonal 
            if (i < block_count - 1) {
                for (int c = 0; c < block_size; ++c) {
                    if (B(r,c) != 0) {
                        value.push_back(B(r,c));
                        col.push_back((i + 1) * block_size + c);
                    }
                }
            }
            row.push_back(value.size());
        }
    }
}

void csr::block_tridiagonal(const csr& A, const csr& B) {
    int block_size = A.row.size() - 1;  
    int block_count = B.row.size() -1;
    int total_size = block_count * block_size;  

    row.clear();
    col.clear();
    value.clear();
    row.push_back(0);

    for (int i = 0; i < block_count; ++i) {
        int block_start = i * block_size;  

        for (int r = 0; r < block_size; ++r) {
            int row_offset = block_start + r;

            // Main diagonal block (A)
            for (int j = A.row[r]; j < A.row[r + 1]; ++j) {
                value.push_back(A.value[j]);
                col.push_back(block_start + A.col[j]);
            }

            // Sub-diagonal block (B)
            if (i > 0) {
                for (int j = B.row[r]; j < B.row[r + 1]; ++j) {
                    value.push_back(B.value[j]);
                    col.push_back((i - 1) * block_size + B.col[j]);
                }
            }

            // Super-diagonal block (B)
            if (i < block_count - 1) {
                for (int j = B.row[r]; j < B.row[r + 1]; ++j) {
                    value.push_back(B.value[j]);
                    col.push_back((i + 1) * block_size + B.col[j]);
                }
            }
            row.push_back(value.size());
        }
    }
}

//thi should be parallel
std::vector<double> csr::matvec(const std::vector<double>& vec) const {
    std::vector<double> result(vec.size(), 0);
    for (size_t i = 0; i < row.size() - 1; ++i) {
        int row_start = row[i];
        int row_end = row[i + 1];
        for (int j = row_start; j < row_end; ++j) {
            result[i] += value[j] * vec[col[j]];
        }
    }
    return result;
}


std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

double dot(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

std::vector<double> GC_solver(const csr& A,const std::vector<double>& b, double tol, int max_iter){
    int n=b.size();
    std::vector<double> M(n,0);
    std::vector<double> r=b - A.matvec(M);
    std::vector<double> p=r;
    std::vector<double> rs_old=r;
    std::vector<double> rs_new;

    for (int iter=0;iter<max_iter;iter++){
        std::vector<double> Ap=A.matvec(p);
        double alpha= dot(rs_old,rs_old) / dot(p,Ap);
        M=M+p*alpha;
        r=r-Ap*alpha;
        rs_new=r;

        if (std::sqrt(dot(rs_new,rs_new))<tol) break;

        p=r+ p* (dot(rs_new,rs_new)/dot(rs_old,rs_old));
        rs_old = rs_new;
    }
    return M;
}


void csr::print_matrix(int numRows, int numCols) const {
    std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols, 0));

    for (int i = 0; i < row.size() - 1; ++i) {
        for (int j = row[i]; j < row[i + 1]; ++j) {
            matrix[i][col[j]] = value[j];
        }
    }

    // Print the matrix
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

// Modify specific element in value vector
void csr::set_value_at(size_t index, double val) {
    if (index >= value.size()) {
        throw std::out_of_range("Index out of range in set_value_at");
    }
    value[index] = val;
}

// Modify specific element in col vector
void csr::set_col_at(size_t index, int val) {
    if (index >= col.size()) {
        throw std::out_of_range("Index out of range in set_col_at");
    }
    col[index] = val;
}

// Modify specific element in row vector
void csr::set_row_at(size_t index, int val) {
    if (index >= row.size()) {
        throw std::out_of_range("Index out of range in set_row_at");
    }
    row[index] = val;
}

double csr::get_value_at(size_t index) const {
    if (index >= value.size()) throw std::out_of_range("Index out of range in get_value_at");
    return value[index];
}

// Likewise for col and row (if you want by-index access)
int csr::get_col_at(size_t index) const {
    if (index >= col.size()) throw std::out_of_range("Index out of range in get_col_at");
    return col[index];
}

int csr::get_row_at(size_t index) const {
    if (index >= row.size()) throw std::out_of_range("Index out of range in get_row_at");
    return row[index];
}

void csr::modify_last_diagonal(double new_value) {
    // Loop through rows to find the last diagonal element
    for (int i = row.size() - 2; i >= 0; --i) { // Looping backwards
        // Search for the diagonal element in the current row
        for (int j = row[i]; j < row[i + 1]; ++j) {
            if (col[j] == i) {  // Diagonal element (i, i)
                // Update the value of the diagonal element
                value[j] = new_value;
                return;
            }
        }
    }
}
double csr::get_value_at(int row1, int col1) const {
    // Find the starting index for the given row
    int row_start = row1 == 0 ? 0 : this->row[row1 - 1];
    int row_end = this->row[row1];

    // Iterate through the values in the given row
    for (int i = row_start; i < row_end; ++i) {
        if (this->col[i] == col1) {
            return this->value[i];  // Return the value if the column matches
        }
    }

    // If the column wasn't found, it means it's a zero in the matrix
    return 0.0;
}

void csr::modify_off_diagonal(int row1, int col1, double value_to_modify) {
    // Find the index of the row in the CSR structure
    int row_start = row[row1];
    int row_end = row[row1 + 1];
    
    // Loop through the row's elements to find the column index
    for (int i = row_start; i < row_end; ++i) {
        if (col[i] == col1) {
            // Found the corresponding column, modify the value
            value[i] += value_to_modify;  // Modify the value at the specified position
            return;
        }
    }
}