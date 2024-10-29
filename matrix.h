#pragma once

#include  <vector>

class Matrix final {

    struct MultithreadMatrixContraints {
        size_t m_maxRowsSum = 300;
        size_t m_maxRowsMult = 200;
        size_t m_maxRowsDet = 2;
        size_t m_maxRowsPerThread = 100;
    };

public:
    Matrix() = default;

    explicit Matrix(size_t rank);

    Matrix(size_t nRows, size_t nCols, double value=0);

    void multithreadingOn();

    void multithreadingOff();

    static Matrix createDiagonal(size_t rank, double value);

    void fill(double value);

    [[nodiscard]] size_t getCols() const;

    [[nodiscard]] size_t getRows() const;

    [[nodiscard]] double det() const;

    friend Matrix operator+(const Matrix &first, const Matrix &second);

    friend Matrix operator-(const Matrix &first, const Matrix &second);

    friend Matrix operator*(const Matrix &first, const Matrix &second);

    double fast_det(const Matrix &mat, size_t num_of_threads) const;

    [[nodiscard]] Matrix fast_subtract_with(const Matrix &another, size_t num_of_threads) const;

    [[nodiscard]] Matrix fast_sum_with(const Matrix &another, size_t num_of_threads) const;

    [[nodiscard]] Matrix fast_multiply_with(const Matrix &another, size_t num_of_threads) const;

    bool operator==(const Matrix &another) const;

    bool operator!=(const Matrix &another) const;

    double &at(size_t i, size_t j);

    void setContraints(const MultithreadMatrixContraints&);

private:
    bool m_multithread = false;

    MultithreadMatrixContraints m_contraints;

    std::vector<std::vector<double>> m_matrix;

    size_t size() const;

    Matrix subtract_with(const Matrix &another) const;

    Matrix sum_with(const Matrix &another) const;

    Matrix multiply_with(const Matrix &another) const;

    static Matrix minor(const Matrix &mat, size_t col_index);

    size_t col_max(const size_t column) const;

    static void triangulation(Matrix &mat, const size_t current, const size_t begin, const size_t end);

    void swap_rows(const size_t i, const size_t j);

    friend class CalculationManager;
};