#ifndef DETERMINANT_POLICIES_HPP
#define DETERMINANT_POLICIES_HPP

#include <vector>
#include <cmath>

template<typename T>
class LaplaceExpansion {
public:
    static T calculate(const std::vector<std::vector<T>>& matrix) {
        return determinant(matrix);
    }

private:
    static T determinant(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        if (N == 1) {
            return matrix[0][0];
        } else {
            T det = 0;
            int sign = 1;
            for (int i = 0; i < N; ++i) {
                std::vector<std::vector<T>> subMatrix = createSubMatrix(matrix, 0, i);

                det += sign * matrix[0][i] * determinant(subMatrix);
                sign = -sign;
            }
            return det;
        }
    }

    static std::vector<std::vector<T>> createSubMatrix(const std::vector<std::vector<T>>& matrix, int excludingRow, int excludingCol) {
        int N = matrix.size();
        std::vector<std::vector<T>> subMatrix;
        for (int i = 0; i < N; ++i) {
            if (i == excludingRow) continue;
            std::vector<T> row;
            for (int j = 0; j < N; ++j) {
                if (j == excludingCol) continue;
                row.push_back(matrix[i][j]);
            }
            subMatrix.push_back(row);
        }
        return subMatrix;
    }
};

template<typename T>
class GaussianElimination {
public:
    static T calculate(std::vector<std::vector<T>> matrix) {
        int M = matrix.size();
        int N = matrix[0].size();
        T det = 1;
        for (int i = 0; i < M; ++i) {
            int maxRow = findPivotRow(matrix, i);
            if (matrix[maxRow][i] == 0) {
                return 0; 
            }

            if (i != maxRow) {
                std::swap(matrix[i], matrix[maxRow]);
                det *= -1; 
            }

            det *= matrix[i][i];
            for (int j = i + 1; j < N; ++j) {
                for (int k = M - 1; k >= i; k--) {
                    matrix[j][k] -= matrix[i][k] * (matrix[j][i] / matrix[i][i]);
                }
            }
        }
        return det;
    }

private:
    static int findPivotRow(const std::vector<std::vector<T>>& matrix, int col) {
        int maxRow = col;
        for (int i = col + 1; i < matrix.size(); ++i) {
            if (std::abs(matrix[i][col]) > std::abs(matrix[maxRow][col])) {
                maxRow = i;
            }
        }
        return maxRow;
    }
};






#endif // DETERMINANT_POLICIES_HPP
