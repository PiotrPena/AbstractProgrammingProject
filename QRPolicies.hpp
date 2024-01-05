#ifndef QR_POLICIES_HPP
#define QR_POLICIES_HPP

#include<tuple>
#include <vector>
#include <cmath>

template<typename T>
class GramSchmidt {
public:
    static std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> calculate(const std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();

        std::vector<std::vector<T>> Q(rows, std::vector<T>(cols, 0));
        std::vector<std::vector<T>> R(cols, std::vector<T>(cols, 0));

        for (int k = 0; k < cols; ++k) {
            for (int i = 0; i < rows; ++i) {
                Q[i][k] = matrix[i][k];
            }
            for (int j = 0; j < k; ++j) {
                T dot_product = 0;
                for (int i = 0; i < rows; ++i) {
                    dot_product += Q[i][j] * matrix[i][k];
                }
                R[j][k] = dot_product;
                for (int i = 0; i < rows; ++i) {
                    Q[i][k] -= R[j][k] * Q[i][j];
                }
            }
            T norm = 0;
            for (int i = 0; i < rows; ++i) {
                norm += Q[i][k] * Q[i][k];
            }
            norm = sqrt(norm);
            R[k][k] = norm;
            if (norm != 0) {
                for (int i = 0; i < rows; ++i) {
                    Q[i][k] /= norm;
                }
            }
        }

        return {Q, R};
    }
};

template<typename T>
class Householder {
public:
    static std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> calculate(const std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();

        std::vector<std::vector<T>> Q(rows, std::vector<T>(rows, 0));
        for (int i = 0; i < rows; ++i) {
            Q[i][i] = 1;
        }

        std::vector<std::vector<T>> R = matrix;

        for (int k = 0; k < cols && k < rows - 1; ++k) {
            std::vector<T> x(rows - k, 0);
            for (int i = k; i < rows; ++i) {
                x[i - k] = R[i][k];
            }

            T norm_x = norm(x);
            if (norm_x == 0) continue;

            x[0] += (x[0] >= 0 ? norm_x : -norm_x);
            T norm_v = norm(x);

            for (T& element : x) {
                element /= norm_v;
            }

            for (int j = k; j < cols; ++j) {
                T sum = 0;
                for (int i = k; i < rows; ++i) {
                    sum += x[i - k] * R[i][j];
                }
                for (int i = k; i < rows; ++i) {
                    R[i][j] -= 2 * x[i - k] * sum;
                }
            }
            for (int j = 0; j < rows; ++j) {
                T sum = 0;
                for (int i = k; i < rows; ++i) {
                    sum += x[i - k] * Q[j][i];
                }
                for (int i = k; i < rows; ++i) {
                    Q[j][i] -= 2 * x[i - k] * sum;
                }
            }
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = i + 1; j < rows; ++j) {
                std::swap(Q[i][j], Q[j][i]);
            }
        }

        return {Q, R};
    }

private:
    static T norm(const std::vector<T>& vec) {
        T sum = 0;
        for (T val : vec) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }
};

template<typename T>
class Givens {
public:
    static std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> calculate(const std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();

        std::vector<std::vector<T>> R = matrix;
        std::vector<std::vector<T>> Q(rows, std::vector<T>(rows, 0));
        for (int i = 0; i < rows; ++i) {
            Q[i][i] = 1;
        }

        for (int j = 0; j < cols; ++j) {
            for (int i = rows - 1; i > j; --i) {
                T a = R[i - 1][j];
                T b = R[i][j];

                if (b != 0) {
                    T hypot = std::sqrt(a * a + b * b);
                    T c = a / hypot;
                    T s = -b / hypot;

                    for (int k = 0; k < cols; ++k) {
                        T tmp1 = R[i - 1][k];
                        T tmp2 = R[i][k];
                        R[i - 1][k] = c * tmp1 + s * tmp2;
                        R[i][k] = -s * tmp1 + c * tmp2;
                    }

                    for (int k = 0; k < rows; ++k) {
                        T tmp1 = Q[k][i - 1];
                        T tmp2 = Q[k][i];
                        Q[k][i - 1] = c * tmp1 - s * tmp2;
                        Q[k][i] = s * tmp1 + c * tmp2;
                    }
                }
            }
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = i + 1; j < rows; ++j) {
                std::swap(Q[i][j], Q[j][i]);
            }
        }

        return {Q, R};
    }
};




#endif // QR_POLICIES_HPP
