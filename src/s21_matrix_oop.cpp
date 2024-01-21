#include "s21_matrix_oop.h"

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) throw std::runtime_error("Invalid rows.");
  rows_ = rows;
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) throw std::runtime_error("Invalid columns.");
  cols_ = cols;
}

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  if (rows <= 0 || cols <= 0)
    throw std::runtime_error("Invalid rows or columns.");
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; ++j) matrix_[i][j] = 0;
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(nullptr) {
  if (other.rows_ <= 0 || other.cols_ <= 0)
    throw std::runtime_error("Invalid rows or columns.");
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (matrix_[i][j] != other.matrix_[i][j]) return false;
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_ || rows_ <= 0 ||
      cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0)
    throw std::runtime_error("Invalid matrices. Cannot add.");

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] += other.matrix_[i][j];
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_ || rows_ <= 0 ||
      cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0)
    throw std::runtime_error("Invalid matrices. Cannot subtract.");

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] -= other.matrix_[i][j];
}

void S21Matrix::MulNumber(const double num) {
  if (rows_ <= 0 || cols_ <= 0)
    throw std::runtime_error("Invalid matrix. Cannot multiply by a number.");

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_ || rows_ <= 0 || cols_ <= 0 || other.rows_ <= 0 ||
      other.cols_ <= 0) {
    throw std::runtime_error("Invalid matrices. Cannot multiply.");
  }

  S21Matrix result(rows_, other.cols_);

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.cols_; j++)
      for (int k = 0; k < cols_; k++)
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];

  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() {
  if (rows_ <= 0 || cols_ <= 0 || matrix_ == nullptr)
    throw std::runtime_error("Invalid matrix. Cannot transpose.");

  S21Matrix result(cols_, rows_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }

  return result;
}

S21Matrix S21Matrix::GetMinor(int row, int column, const S21Matrix& A) {
  if (row <= 0 || column <= 0 || A.GetRows() != A.GetCols()) {
    throw std::runtime_error("Invalid matrix. Cannot calculate minor.");
  }

  S21Matrix result(A.GetRows() - 1, A.GetCols() - 1);

  int tmp_row = 0;
  for (int i = 0; i < A.GetRows(); ++i) {
    if (i == row - 1) {
      continue;
    }

    int tmp_col = 0;
    for (int j = 0; j < A.GetCols(); ++j) {
      if (j == column - 1) {
        continue;
      }
      result(tmp_row, tmp_col) = A(i, j);
      ++tmp_col;
    }
    ++tmp_row;
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::runtime_error("Non-square matrix. Cannot compute complements.");

  S21Matrix result(rows_, cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      S21Matrix minor = GetMinor(i + 1, j + 1, *this);
      double det = minor.Determinant();
      result.matrix_[i][j] = pow(-1, i + j) * det;
    }
  }

  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_)
    throw std::runtime_error("Non-square matrix. Cannot compute determinant.");

  if (rows_ == 1)
    return matrix_[0][0];
  else if (rows_ == 2)
    return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  else {
    double det = 0.0;
    for (int i = 0; i < rows_; ++i) {
      S21Matrix minor = GetMinor(1, i + 1, *this);
      double minor_det = minor.Determinant();
      det += pow(-1, i) * matrix_[0][i] * minor_det;
    }
    return det;
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (rows_ != cols_ || det == 0)
    throw std::runtime_error(
        "Matrix must be square and have a non-zero determinant for inversion.");

  S21Matrix complements = CalcComplements();
  S21Matrix transposed = complements.Transpose();
  S21Matrix inverse = transposed * (1.0 / det);

  return inverse;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  try {
    S21Matrix result = *this;
    result.SumMatrix(other);
    return result;
  } catch (const std::exception& e) {
    // Обработка исключения или его пробрасывание наверх, если необходимо
    throw e;
  }
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  try {
    S21Matrix result = *this;
    result.SubMatrix(other);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  try {
    S21Matrix result = *this;
    result.MulMatrix(other);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix S21Matrix::operator*(const double num) const {
  try {
    S21Matrix result = *this;
    result.MulNumber(num);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  try {
    S21Matrix result = *this;
    return result.EqMatrix(other);
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) return *this;

  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      if (matrix_[i] != nullptr) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = new double*[rows_];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) matrix_[i][j] = other.matrix_[i][j];
  }

  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  try {
    *this = *this + other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  try {
    *this = *this - other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  try {
    *this = *this * other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator*=(double num) {
  try {
    *this = *this * num;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("Matrix index out of range");
  return matrix_[i][j];
}

const double& S21Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("Matrix index out of range");
  return matrix_[i][j];
}