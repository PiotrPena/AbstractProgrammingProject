#ifndef LU_POLICIES_HPP
#define LU_POLICIES_HPP

#include<tuple>

template<typename T>
class Doolittle {
public:
    static std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, std::vector<T>, std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        std::vector<std::vector<T>> L(N, std::vector<T>(N, 0));
        std::vector<std::vector<T>> U(N, std::vector<T>(N, 0));

        std::vector<T> rowPermutation(N), colPermutation(N);
        for (int i = 0; i < N; ++i) {
            rowPermutation[i] = colPermutation[i] = i;
        }

        for (int i = 0; i < N; ++i) {
            for (int k = i; k < N; ++k) {
                T sum = 0;
                for (int j = 0; j < i; ++j) {
                    sum += (L[i][j] * U[j][k]);
                }
                U[i][k] = matrix[i][k] - sum;
            }

            for (int k = i; k < N; ++k) {
                if (i == k)
                    L[i][i] = 1;
                else {
                    T sum = 0;
                    for (int j = 0; j < i; ++j) {
                        sum += (L[k][j] * U[j][i]);
                    }
                    L[k][i] = (matrix[k][i] - sum) / U[i][i];
                }
            }
        }

        return {L, U, rowPermutation, colPermutation};
    }
};

template<typename T>
class Crout {
public:
    static std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, std::vector<T>, std::vector<T>>  calculate(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        std::vector<std::vector<T>> L(N, std::vector<T>(N, 0));
        std::vector<std::vector<T>> U(N, std::vector<T>(N, 0));

        std::vector<T> rowPermutation(N), colPermutation(N);
        for (int i = 0; i < N; ++i) {
            rowPermutation[i] = colPermutation[i] = i;
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; j++) {
                T sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = matrix[i][j] - sum;
            }

            for (int j = i; j < N; j++) {
                if (i == j)
                    U[i][j] = 1;
                else {
                    T sum = 0;
                    for (int k = 0; k < i; k++) {
                        sum += L[i][k] * U[k][j];
                    }
                    U[i][j] = (matrix[i][j] - sum) / L[i][i];
                }
            }
        }

        return {L, U, rowPermutation, colPermutation};
    }
};

template<typename T>
class GaussianFullPivoting {
public:
    static std::tuple<std::vector<std::vector<T>>, std::vector<std::vector<T>>, std::vector<int>, std::vector<int>> calculate(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        std::vector<std::vector<T>> L(N, std::vector<T>(N, 0));
        std::vector<std::vector<T>> U = matrix;
        std::vector<int> rowPermutation(N), colPermutation(N);

        for (int i = 0; i < N; ++i) {
            rowPermutation[i] = i;
            colPermutation[i] = i;
        }

        for (int i = 0; i < N; ++i) {
            T maxVal = 0;
            int maxRow = i, maxCol = i;
            for (int row = i; row < N; ++row) {
                for (int col = i; col < N; ++col) {
                    if (std::abs(U[row][col]) > maxVal) {
                        maxVal = std::abs(U[row][col]);
                        maxRow = row;
                        maxCol = col;
                    }
                }
            }

            std::swap(U[i], U[maxRow]);
            std::swap(rowPermutation[i], rowPermutation[maxRow]);
            for (int k = 0; k < N; ++k) {
                std::swap(U[k][i], U[k][maxCol]);
            }
            std::swap(colPermutation[i], colPermutation[maxCol]);

            for (int j = i + 1; j < N; ++j) {
                L[j][i] = U[j][i] / U[i][i];
                for (int k = i; k < N; ++k) {
                    U[j][k] -= L[j][i] * U[i][k];
                }
            }
        }

        for (int i = 0; i < N; ++i) {
            L[i][i] = 1;
        }

        return {L, U, rowPermutation, colPermutation};
    }
};


#endif // LU_POLICIES_HPP
