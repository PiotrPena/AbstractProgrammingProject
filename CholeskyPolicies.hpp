#ifndef CHOLESKY_POLICIES_HPP
#define CHOLESKY_POLICIES_HPP

template<typename T>
class Cholesky {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        std::vector<std::vector<T>> L(N, std::vector<T>(N, 0));

        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++) {
                T sum = 0;
                if (j == i) {
                    for (int k = 0; k < j; k++) {
                        sum += L[j][k] * L[j][k];
                    }
                    L[j][j] = std::sqrt(matrix[j][j] - sum);
                } else {
                    for (int k = 0; k < j; k++) {
                        sum += L[i][k] * L[j][k];
                    }
                    L[i][j] = (matrix[i][j] - sum) / L[j][j];
                }
            }
        }

        return L;
    }
};

template<typename T>
class RecursiveCholesky {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix) {
        int N = matrix.size();
        std::vector<std::vector<T>> L(N, std::vector<T>(N, 0));
        decompose(matrix, L, 0, N);
        return L;
    }

private:
    static void decompose(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& L, int start, int N) {
        if (N <= 1) {
            L[start][start] = std::sqrt(A[start][start]);
            return;
        }

        int half = N / 2;

        decompose(A, L, start, half);

        for (int i = start + half; i < start + N; ++i) {
            for (int j = start; j < start + half; ++j) {
                T sum = 0;
                for (int k = start; k < j; ++k) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }

        for (int i = start + half; i < start + N; ++i) {
            for (int j = start + half; j <= i; ++j) {
                T sum = 0;
                for (int k = start; k < j; ++k) {
                    sum += L[i][k] * L[j][k];
                }
                if (i == j) {
                    L[i][j] = std::sqrt(A[i][j] - sum);
                } else {
                    L[i][j] = (A[i][j] - sum) / L[j][j];
                }
            }
        }
    }
};




#endif // Cholesky_POLICIES_HPP
