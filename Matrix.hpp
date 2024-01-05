#ifndef MATRIX_HPP
#define MATRIX_HPP

#include<tuple>

#include "DeterminantPolicies.hpp"
#include "InversePolicies.hpp"
#include "MultiplicationPolicies.hpp"
#include "LUPolicies.hpp"
#include "QRPolicies.hpp"
#include "CholeskyPolicies.hpp"
#include "EigenvaluesPolicies.hpp"
#include "SolvingPolicies.hpp"
#include "Concepts.hpp"

//Struct for Policies
template<typename T>
struct MatrixPolicies {
    using DeterminantPolicy = LaplaceExpansion<T>;
    using InversionPolicy = ClassicalAdjoint<T>;
    using MultiplicationPolicy = StandardMatrixMultiplication<T>;
    using LUPolicy = Doolittle<T>;
    using QRPolicy = Householder<T>;
    using CholeskyPolicy = Cholesky<T>;
    using EigenvaluePolicy = PowerIteration<T>;
    using SolvingPolicy = GaussianEliminationSolver<T>;
    using SolvingDecomposePolicy = QRSolver<T>;
    using SolvingIterativePolicy = GaussSeidelSolver<T>;
};

template <int M, int N, typename T, typename Policies = MatrixPolicies<T>>
class Matrix {
public:
    
//====================CONSTRUCTORS====================================

    Matrix() {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i][j] = T();
            }
        }
    }

    Matrix(const T (&initData)[M][N]) {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i][j] = initData[i][j];
            }
        }
    }

    Matrix(const std::vector<std::vector<T>>& initData) {
        if (initData.size() != M || (initData.size() > 0 && initData[0].size() != N)) {
            throw std::invalid_argument("Invalid dimensions for matrix initialization.");
        }

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i][j] = initData[i][j];
            }
        }
    }

    Matrix(const std::vector<T>& initData) {
        if (initData.size() != M || N != 1) {
            throw std::invalid_argument("Invalid dimensions for matrix initialization.");
        }

        for (int i = 0; i < M; ++i) {
            data[i][0] = initData[i];
        }
    }

//====================================METHODS=======================================================

    // Method to transpose the matrix
    Matrix<N, M, T> transpose() const {
        Matrix<N, M, T> result;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result.data[j][i] = this->data[i][j];
            }
        }
        return result;
    }

    // Method to negate the matrix
    Matrix<M, N, T> negate() const requires Negatable<T> {
        Matrix<M, N, T> result;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result.data[i][j] = -this->data[i][j];
            }
        }
        return result;
    }

    // Method to convert the fixed-size array to a vector of vectors
    std::vector<std::vector<T>> toVectorMatrix() const {
        std::vector<std::vector<T>> vecMatrix(M, std::vector<T>(N));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                vecMatrix[i][j] = data[i][j];
            }
        }
        return vecMatrix;
    }

    // Method for determinant calculation
    T determinant() const requires SquareMatrix<M, N, T> && Arithmetic<T>{
        auto vecMatrix = toVectorMatrix();
        return Policies::DeterminantPolicy::calculate(vecMatrix);
    }

    // Method for inverting a matrix
    Matrix<M, N, T> inverse() const requires SquareMatrix<M, N, T> && Arithmetic<T> {
        auto vecMatrix = toVectorMatrix();
        auto result = Policies::InversionPolicy::calculate(vecMatrix);
        return Matrix<M, N, T, Policies>(result);
    }

    // Method for matrix multiplication 
    template<int P>
    Matrix<M, P, T, Policies> multiply(const Matrix<N, P, T, Policies>& other) const requires Arithmetic<T> {
        auto vecMatrix1 = toVectorMatrix();
        auto vecMatrix2 = other.toVectorMatrix();
        auto result = Policies::MultiplicationPolicy::calculate(vecMatrix1, vecMatrix2);
        return Matrix<M, P, T, Policies>(result);
    }

    // Method for calculating trace
    T trace() const requires SquareMatrix<M, N, T> {

        T traceSum = 0;
        for (int i = 0; i < M; ++i) {
            traceSum += data[i][i];
        }
        return traceSum;
    }

    // Method for LU decomposition
    std::pair<Matrix<M, N, T, Policies>, Matrix<M, N, T, Policies>> luDecomposition() const requires SquareMatrix<M, N, T> && Arithmetic<T>  {
        auto vecMatrix = toVectorMatrix();
        auto [L, U, row, column] = Policies::LUPolicy::calculate(vecMatrix);
        return {Matrix<M, N, T, Policies>(L), Matrix<M, N, T, Policies>(U)};
    }

    // Method for QR decomposition
    std::pair<Matrix<M, N, T, Policies>, Matrix<N, N, T, Policies>> qrDecomposition() const requires Arithmetic<T> {
        auto vecMatrix = toVectorMatrix();
        auto [Q, R] = Policies::QRPolicy::calculate(vecMatrix);
        return {Matrix<M, N, T, Policies>(Q), Matrix<N, N, T, Policies>(R)};
    }

    // Method for Cholesky Decomposition
    Matrix<M, M, T, Policies> choleskyDecomposition() const requires SquareMatrix<M, N, T> && Arithmetic<T> {
        auto vecMatrix = toVectorMatrix();  
        auto L = Policies::CholeskyPolicy::calculate(vecMatrix);
        return Matrix<M, M, T, Policies>(L);
    }

    // Method for computing eigenvalue decomposition
    T eigenvalueDecomposition() const requires SquareMatrix<M, N, T> && Arithmetic<T> && ComparableWithTolerance<T>{
        auto vecMatrix = toVectorMatrix();
        auto [e, eigenvalues] = Policies::EigenvaluePolicy::calculate(vecMatrix);
        
        return e;
    }

    // Method for gaussian solving
    Matrix<M, 1, T, Policies> solve(const Matrix<M, 1, T, Policies>& b) const requires Arithmetic<T> {
        auto A = this->toVectorMatrix();
        auto vecB = b.toVectorMatrix();
        auto solution = Policies::SolvingPolicy::solve(A, vecB[0]);
        return Matrix<M, 1, T, Policies>(solution);
    }

    // Method for decomposition-based solving
    Matrix<M, 1, T, Policies> solveWithDecompose(const Matrix<M, 1, T, Policies>& b) const requires Arithmetic<T> {
        auto A = toVectorMatrix();
        auto vecB = b.toVectorMatrix();
        auto x = Policies::SolvingDecomposePolicy::solve(A, vecB[0]);
        return Matrix<M, 1, T, Policies>(x);
    }

    // Method for iterative solving
    Matrix<M, 1, T, Policies> solveIteratively(const Matrix<M, 1, T, Policies>& b, T tolerance = 1e-7, int maxIterations = 1000) const requires Arithmetic<T> && ComparableWithTolerance<T>{
        auto A = toVectorMatrix();
        auto vecB = b.toVectorMatrix();
        auto solution = Policies::SolvingIterativePolicy::solve(A, vecB[0], tolerance, maxIterations);
        return Matrix<M, 1, T, Policies>(solution);
    }

    // Method for QR decomposition
    Matrix<M, N, T, Policies> ortogonalize() const requires Arithmetic<T> {
        auto vecMatrix = toVectorMatrix();
        auto [Q, R] = Policies::QRPolicy::calculate(vecMatrix);
        return Matrix<M, N, T, Policies>(Q);
    }

    //=================================OPERATORS====================================================================================

    Matrix<M, N, T> operator*(T scalar) const {
        Matrix<M, N, T> result;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result.data[i][j] = this->data[i][j] * scalar;
            }
        }
        return result;
    }

    T& operator()(int row, int col) {
        return data[row][col];
    }

    const T& operator()(int row, int col) const {
        return data[row][col];
    }

    Matrix<M, N, T> operator-() const requires Negatable<T> {
        Matrix<M, N, T> result;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result.data[i][j] = -this->data[i][j];
            }
        }
        return result;
    }

    Matrix<M, N, T> operator+(const Matrix<M, N, T>& rhs) const requires Addable<T> {
        Matrix<M, N, T> result;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result.data[i][j] = this->data[i][j] + rhs.data[i][j];
            }
        }
        return result;
    }

    template<int P>
    Matrix<M, P, T, Policies> operator*(const Matrix<N, P, T, Policies>& other) const requires Arithmetic<T> {
        return this->multiply(other);
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix<M, N, T>& matrix) requires Streamable<T> {
        std::cout << std::endl;
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                os << matrix.data[i][j];
                if (j < N - 1) os << " ";
            }
            if (i < M - 1) os << "\n";
        }
        return os;
    }

private:
    T data[M][N];
};

#endif // MATRIX_HPP