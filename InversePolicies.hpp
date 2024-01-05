#ifndef INVERSE_POLICIES_HPP
#define INVERSE_POLICIES_HPP

#include <vector>
#include <cmath>

template<typename T>
class RowReduction {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix) {
        int M = matrix.size();
        int N = matrix[0].size();

        std::vector<std::vector<T>> augmented(M, std::vector<T>(2 * N));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                augmented[i][j] = matrix[i][j];
                augmented[i][j + N] = (i == j) ? static_cast<T>(1) : static_cast<T>(0);
            }
        }


        for (int i = 0; i < M; ++i) {
            if (augmented[i][i] == 0) {
                int swapRow = findPivotRow(augmented, i, M);
                if (swapRow < 0) {
                    throw std::runtime_error("Matrix is singular and cannot be inverted.");
                }
                std::swap(augmented[i], augmented[swapRow]);
            }

            T pivotValue = augmented[i][i];
            for (int j = 0; j < 2 * N; ++j) {
                augmented[i][j] /= pivotValue;
            }

            for (int j = 0; j < M; ++j) {
                if (j != i && augmented[j][i] != 0) {
                    T factor = augmented[j][i];
                    for (int k = 0; k < 2 * N; ++k) {
                        augmented[j][k] -= factor * augmented[i][k];
                    }
                }
            }
        }

        std::vector<std::vector<T>> inverse(M, std::vector<T>(N));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                inverse[i][j] = augmented[i][j + N];
            }
        }
        return inverse;
    }

private:
    static int findPivotRow(const std::vector<std::vector<T>>& augmented, int col, int M) {
        int maxRow = col;
        T maxVal = std::abs(augmented[col][col]);
        for (int i = col + 1; i < M; ++i) {
            T val = std::abs(augmented[i][col]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }
        return maxRow;
    }
};

template<typename T>
class ClassicalAdjoint {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix) {
        int M = matrix.size();
        int N = matrix[0].size();

        T det = determinant(matrix);
        if (det == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        std::vector<std::vector<T>> cofactors(M, std::vector<T>(N));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                cofactors[i][j] = pow(-1, i + j) * determinant(minor(matrix, i, j));
            }
        }
        std::vector<std::vector<T>> adjugate = transpose(cofactors);

        std::vector<std::vector<T>> inverse(M, std::vector<T>(N));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                inverse[i][j] = adjugate[i][j] / det;
            }
        }
        return inverse;
    }

private:
    static T determinant(const std::vector<std::vector<T>>& matrix) {
        int size = matrix.size();

        if (size == 0) {
            throw std::invalid_argument("Empty matrix has no determinant.");
        }

        for (const auto& row : matrix) {
            if (row.size() != size) {
                throw std::invalid_argument("Determinant can only be calculated for square matrices.");
            }
        }

        if (size == 1) {
            return matrix[0][0];
        }

        if (size == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }

        T det = 0;
        int sign = 1;
        for (int i = 0; i < size; ++i) {
            std::vector<std::vector<T>> subMatrix = minor(matrix, 0, i);
            det += sign * matrix[0][i] * determinant(subMatrix);
            sign = -sign;
        }

        return det;
    }


    static std::vector<std::vector<T>> minor(const std::vector<std::vector<T>>& matrix, int row, int col) {
        int M = matrix.size();
        std::vector<std::vector<T>> result(M - 1, std::vector<T>(M - 1));
        for (int i = 0, m = 0; i < M; ++i) {
            if (i == row) continue;
            for (int j = 0, n = 0; j < M; ++j) {
                if (j == col) continue;
                result[m][n] = matrix[i][j];
                n++;
            }
            m++;
        }
        return result;
    }

    static std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& matrix) {
        int M = matrix.size();
        int N = matrix[0].size();
        std::vector<std::vector<T>> result(N, std::vector<T>(M));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }
};

#endif // INVERSE_POLICIES_HPP