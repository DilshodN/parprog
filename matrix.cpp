#include  "matrix.h"
#include <thread>
#include <future>


size_t Matrix::size() const
{
    return getCols() * getRows();
}

Matrix::Matrix(size_t rank) 
{
    m_matrix.resize(rank);
    for (size_t i = 0; i < rank; i++) {
        m_matrix[i].resize(rank);
    }
}

Matrix::Matrix(size_t nRows, size_t nCols, double value)
{
    m_matrix.resize(nRows);
    for (size_t i = 0; i < nRows; i++) {
        m_matrix[i].resize(nCols);
        for (size_t j = 0; j < nCols; j++) {
            m_matrix[i][j] = value;
        }
    }
}

void Matrix::multithreadingOn() {
    m_multithread = true;
}

void Matrix::multithreadingOff() {
    m_multithread = false;
}

Matrix operator+(const Matrix &first, const Matrix &second) {
    return first.m_multithread or second.m_multithread ?
           first.fast_sum_with(second, std::thread::hardware_concurrency()) : first.sum_with(second);
}

Matrix operator*(const Matrix &first, const Matrix &second) {
    return first.m_multithread or second.m_multithread ?
           first.fast_multiply_with(second, std::thread::hardware_concurrency()) : first.multiply_with(second);
}


Matrix operator-(const Matrix &first, const Matrix &second) {
    return first.m_multithread or second.m_multithread ?
           first.fast_subtract_with(second, std::thread::hardware_concurrency()) : first.subtract_with(second);
}

double Matrix::det() const {
    if(m_multithread){
//        double det = 1;
//        fast_det(*this, std::thread::hardware_concurrency(), det);
        return fast_det(*this, std::thread::hardware_concurrency());
    }
    else{
        return fast_det(*this, 1);
    }
}

bool Matrix::operator==(const Matrix &another) const {
    if (getRows() != another.getRows() || getCols() != another.getCols()) return false;
    for (size_t i = 0; i < getRows(); i++) {
        if (m_matrix[i] != another.m_matrix[i]) return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix &another) const {
    return !(*this == another);
}


Matrix Matrix::createDiagonal(size_t rank, double value) {
    auto diagonal = Matrix(rank, rank);

    for (size_t i = 0; i < rank; i++) {
        diagonal.m_matrix[i][i] = value;
    }

    return diagonal;
}


void Matrix::fill(double value) {
    for (size_t i = 0; i < getRows(); i++) {
        for (size_t j = 0; j < getCols(); j++) {
            m_matrix[i][j] = value;
        }
    }
}

size_t Matrix::getCols() const{
    return m_matrix[0].size();
}

size_t Matrix::getRows() const {
    return m_matrix.size();
}

Matrix Matrix::sum_with(const Matrix &another) const {
    const size_t rows = getRows();
    const size_t cols = getCols();
    if (rows != another.getRows() || cols != another.getCols()) {
        return {};
    }
    Matrix sum(rows, cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            sum.m_matrix[i][j] = m_matrix[i][j] + another.m_matrix[i][j];
        }
    }
    return sum;
}

Matrix Matrix::subtract_with(const Matrix &another) const {
    const size_t rows = getRows();
    const size_t cols = getCols();
    if (rows != another.getRows() or cols != another.getCols()) {
        return {};
    }
    Matrix subtract(rows, cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            subtract.m_matrix[i][j] = m_matrix[i][j] - another.m_matrix[i][j];
        }
    }
    return subtract;
}

Matrix Matrix::multiply_with(const Matrix &another) const {
    const size_t rows = getRows();
    const size_t cols = getCols();
    Matrix result(rows, another.getCols());
    for (size_t i = 0; i < result.getRows(); i++) {
        for (size_t j = 0; j < result.getCols(); j++) {
            for (size_t k = 0; k < result.getCols(); k++) {
                result.m_matrix[i][j] += m_matrix[i][k] * another.m_matrix[k][j];
            }
        }
    }
    return result;
}


Matrix Matrix::minor(const Matrix &mat, size_t col_index) {
    auto col_idx = static_cast<long>(col_index);
    Matrix sub_mat(mat.getRows() - 1);
    for (size_t i = 1; i < mat.getCols(); i++) {
        const std::vector<double> &temp_row = mat.m_matrix[i];
        std::copy(temp_row.begin(), temp_row.begin() + col_idx, sub_mat.m_matrix[i - 1].begin());
        std::copy(temp_row.begin() + col_idx + 1, temp_row.end(), sub_mat.m_matrix[i - 1].begin() + col_idx);
    }
    return sub_mat;
}

size_t Matrix::col_max(const size_t column) const{
    double max = std::abs(m_matrix[column][column]);
    auto max_pos = column;
    for (auto i = column + 1; i < getRows(); ++i) {
        double element = std::abs(m_matrix[i][column]);
        if (element > max) {
            max = element;
            max_pos = i;
        }
    }
    return max_pos;
}

double& Matrix::at(size_t i, size_t j) {
    return m_matrix[i][j];
}

void Matrix::swap_rows(const size_t i, const size_t j) {
    std::swap(m_matrix[i], m_matrix[j]);
}

void Matrix::triangulation(Matrix& mat, const size_t current, const size_t begin, const size_t end)
{
    for (auto j = begin; j < end; ++j) {
        const auto mul = - mat.at(j, current) / mat.at(current, current);
        for (auto k = current; k < mat.getRows(); ++k) {
            mat.at(j, k) += mat.at(current, k) * mul;
        }
    }
}

void Matrix::setContraints(const Matrix::MultithreadMatrixContraints& constraints) 
{
    m_contraints = constraints;
}






