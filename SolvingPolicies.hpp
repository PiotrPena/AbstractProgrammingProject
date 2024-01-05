#ifndef SOLVING_POLICIES_HPP
#define SOLVING_POLICIES_HPP

#include <vector>
#include <stdexcept>
#include <cmath>


template<typename T>
class GaussianEliminationSolver {
public:
    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
        int n = A.size();
        
        std::vector<std::vector<T>> matrix = A;
        std::vector<T> vec = b;

        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(matrix[i], matrix[maxRow]);
            std::swap(vec[i], vec[maxRow]);

            if (std::abs(matrix[i][i]) < 1e-9) {
                throw std::runtime_error("Singular matrix encountered during Gaussian Elimination.");
            }
            for (int j = i + 1; j < n; ++j) {
                T factor = matrix[j][i] / matrix[i][i];
                vec[j] -= factor * vec[i];
                for (int k = i; k < n; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }

        std::vector<T> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = vec[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= matrix[i][j] * x[j];
            }
            x[i] /= matrix[i][i];
        }

        return x;
    }
};

template<typename T>
class LUDecomposition {
public:
    static std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> decompose(const std::vector<std::vector<T>>& A) {
        int n = A.size();
        std::vector<std::vector<T>> L(n, std::vector<T>(n, 0));
        std::vector<std::vector<T>> U(n, std::vector<T>(n, 0));

        for (int i = 0; i < n; i++) {
            for (int k = i; k < n; k++) {
                T sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[i][j] * U[j][k]);

                U[i][k] = A[i][k] - sum;
            }

            for (int k = i; k < n; k++) {
                if (i == k)
                    L[i][i] = 1;
                else {
                    T sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (L[k][j] * U[j][i]);

                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }
        }
        return {L, U};
    }

    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
        auto [L, U] = decompose(A);
        int n = L.size();

        std::vector<T> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = b[i];
            for (int j = 0; j < i; ++j) {
                y[i] -= L[i][j] * y[j];
            }
            y[i] /= L[i][i];
        }

        std::vector<T> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= U[i][j] * x[j];
            }
            x[i] /= U[i][i];
        }

        return x;
    }
};



template<typename T>
class CholeskySolver {
public:
    static std::vector<std::vector<T>> decompose(const std::vector<std::vector<T>>& A) {
        int n = A.size();
        std::vector<std::vector<T>> L(n, std::vector<T>(n, 0));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                T sum = 0;
                if (j == i) {
                    for (int k = 0; k < j; k++)
                        sum += L[j][k] * L[j][k];
                    L[j][j] = sqrt(A[j][j] - sum);
                } else {
                    for (int k = 0; k < j; k++)
                        sum += L[i][k] * L[j][k];
                    L[i][j] = (A[i][j] - sum) / L[j][j];
                }
            }
        }

        return L;
    }

    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
        std::cout << "Cholesky" << std::endl;
        auto L = decompose(A);
        int n = L.size();

        std::vector<T> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = b[i];
            for (int j = 0; j < i; ++j) {
                y[i] -= L[i][j] * y[j];
            }
            y[i] /= L[i][i];
        }

        std::vector<T> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= L[j][i] * x[j];
            }
            x[i] /= L[i][i];
        }

        return x;
    }
};


template<typename T>
class QRSolver {
public:
    static std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> decompose(const std::vector<std::vector<T>>& A) {
        int n = A.size();
        std::vector<std::vector<T>> Q(n, std::vector<T>(n, 0));
        std::vector<std::vector<T>> R(n, std::vector<T>(n, 0));
        std::vector<std::vector<T>> A_copy = A;

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                R[k][k] += A_copy[i][k] * A_copy[i][k];
            }
            R[k][k] = sqrt(R[k][k]);
            for (int i = 0; i < n; i++) {
                Q[i][k] = A_copy[i][k] / R[k][k];
            }
            for (int j = k + 1; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    R[k][j] += Q[i][k] * A_copy[i][j];
                }
                for (int i = 0; i < n; i++) {
                    A_copy[i][j] -= Q[i][k] * R[k][j];
                }
            }
        }

        return {Q, R};
    }

    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
        auto [Q, R] = decompose(A);
        int n = Q.size();

        std::vector<T> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = 0;
            for (int j = 0; j < n; ++j) {
                y[i] += Q[j][i] * b[j];
            }
        }

        std::vector<T> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= R[i][j] * x[j];
            }
            if (R[i][i] == 0) {
                throw std::runtime_error("Singular matrix");
            }
            x[i] /= R[i][i];
        }

        return x;
    }
};


template<typename T>
class JacobiSolver {
public:
    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b, T tolerance = 1e-7, int maxIterations = 1000) {
        int n = A.size();
        std::vector<T> x(n, 0);
        std::vector<T> x_old(n, 0);

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            for (int i = 0; i < n; ++i) {
                T sigma = 0;
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        sigma += A[i][j] * x_old[j];
                    }
                }
                x[i] = (b[i] - sigma) / A[i][i];
            }
            T error = 0;
            for (int i = 0; i < n; ++i) {
                error += std::abs(x[i] - x_old[i]);
            }

            if (error < tolerance) {
                break;
            }

            x_old = x;
        }

        return x;
    }
};

template<typename T>
class GaussSeidelSolver {
public:
    static std::vector<T> solve(const std::vector<std::vector<T>>& A, const std::vector<T>& b, T tolerance = 1e-7, int maxIterations = 1000) {
        int n = A.size();
        std::vector<T> x(n, 0);
        std::vector<T> x_old(n, 0);

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            for (int i = 0; i < n; ++i) {
                T sigma = 0;
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        sigma += A[i][j] * x[j];
                    }
                }
                x_old[i] = x[i];
                x[i] = (b[i] - sigma) / A[i][i];
            }
            T error = 0;
            for (int i = 0; i < n; ++i) {
                error += std::abs(x[i] - x_old[i]);
            }

            if (error < tolerance) {
                break;
            }
        }

        return x;
    }
};






#endif // SOLVING_POLICIES_HPP
