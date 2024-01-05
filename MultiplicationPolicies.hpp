#ifndef MULTIPLICATION_POLICIES_HPP
#define MULTIPLICATION_POLICIES_HPP

#include <vector>
#include <cmath>

template<typename T>
class StandardMatrixMultiplication {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& matrixA, const std::vector<std::vector<T>>& matrixB) {
        int rowsA = matrixA.size();
        int colsA = matrixA[0].size();
        int rowsB = matrixB.size();
        int colsB = matrixB[0].size();

        std::vector<std::vector<T>> result(rowsA, std::vector<T>(colsB, 0));
        for (int i = 0; i < rowsA; ++i) {
            for (int j = 0; j < colsB; ++j) {
                for (int k = 0; k < colsA; ++k) {
                    result[i][j] += matrixA[i][k] * matrixB[k][j];
                }
            }
        }
        return result;
    }
};

template<typename T>
class DivideAndConquerMultiplication {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n, 0));

        if (n == 1) {
            result[0][0] = A[0][0] * B[0][0];
        } else {
            int newSize = n / 2;
            std::vector<std::vector<T>> 
                a11(newSize, std::vector<T>(newSize)),
                a12(newSize, std::vector<T>(newSize)),
                a21(newSize, std::vector<T>(newSize)),
                a22(newSize, std::vector<T>(newSize)),
                b11(newSize, std::vector<T>(newSize)),
                b12(newSize, std::vector<T>(newSize)),
                b21(newSize, std::vector<T>(newSize)),
                b22(newSize, std::vector<T>(newSize));

            for (int i = 0; i < newSize; i++) {
                for (int j = 0; j < newSize; j++) {
                    a11[i][j] = A[i][j];
                    a12[i][j] = A[i][j + newSize];
                    a21[i][j] = A[i + newSize][j];
                    a22[i][j] = A[i + newSize][j + newSize];

                    b11[i][j] = B[i][j];
                    b12[i][j] = B[i][j + newSize];
                    b21[i][j] = B[i + newSize][j];
                    b22[i][j] = B[i + newSize][j + newSize];
                }
            }

            auto c11 = add(calculate(a11, b11), calculate(a12, b21));
            auto c12 = add(calculate(a11, b12), calculate(a12, b22));
            auto c21 = add(calculate(a21, b11), calculate(a22, b21));
            auto c22 = add(calculate(a21, b12), calculate(a22, b22));

            for (int i = 0; i < newSize; i++) {
                for (int j = 0; j < newSize; j++) {
                    result[i][j] = c11[i][j];
                    result[i][j + newSize] = c12[i][j];
                    result[i + newSize][j] = c21[i][j];
                    result[i + newSize][j + newSize] = c22[i][j];
                }
            }
        }
        return result;
    }

private:
    static std::vector<std::vector<T>> add(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = A[i][j] + B[i][j];
            }
        }
        return result;
    }
};

template<typename T>
class StrassenMultiplication {
public:
    static std::vector<std::vector<T>> calculate(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n, 0));

        if (n <= 2) {
            return standardMultiplication(A, B);
        } else {
            int newSize = n / 2;
            std::vector<std::vector<T>> 
                a11(newSize, std::vector<T>(newSize)),
                a12(newSize, std::vector<T>(newSize)),
                a21(newSize, std::vector<T>(newSize)),
                a22(newSize, std::vector<T>(newSize)),
                b11(newSize, std::vector<T>(newSize)),
                b12(newSize, std::vector<T>(newSize)),
                b21(newSize, std::vector<T>(newSize)),
                b22(newSize, std::vector<T>(newSize));

            splitMatrix(A, a11, 0, 0);
            splitMatrix(A, a12, 0, newSize);
            splitMatrix(A, a21, newSize, 0);
            splitMatrix(A, a22, newSize, newSize);

            splitMatrix(B, b11, 0, 0);
            splitMatrix(B, b12, 0, newSize);
            splitMatrix(B, b21, newSize, 0);
            splitMatrix(B, b22, newSize, newSize);

            auto p1 = calculate(add(a11, a22), add(b11, b22));
            auto p2 = calculate(add(a21, a22), b11);
            auto p3 = calculate(a11, subtract(b12, b22));
            auto p4 = calculate(a22, subtract(b21, b11));
            auto p5 = calculate(add(a11, a12), b22);
            auto p6 = calculate(subtract(a21, a11), add(b11, b12));
            auto p7 = calculate(subtract(a12, a22), add(b21, b22));

            auto c11 = add(subtract(add(p1, p4), p5), p7);
            auto c12 = add(p3, p5);
            auto c21 = add(p2, p4);
            auto c22 = add(subtract(add(p1, p3), p2), p6);

            joinMatrix(c11, result, 0, 0);
            joinMatrix(c12, result, 0, newSize);
            joinMatrix(c21, result, newSize, 0);
            joinMatrix(c22, result, newSize, newSize);
        }
        return result;
    }

private:
    static void splitMatrix(const std::vector<std::vector<T>>& source, std::vector<std::vector<T>>& destination, int startX, int startY) {
        int n = destination.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                destination[i][j] = source[i + startX][j + startY];
            }
        }
    }

    static void joinMatrix(const std::vector<std::vector<T>>& source, std::vector<std::vector<T>>& destination, int startX, int startY) {
        int n = source.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                destination[i + startX][j + startY] = source[i][j];
            }
        }
    }

    static std::vector<std::vector<T>> add(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = A[i][j] + B[i][j];
            }
        }
        return result;
    }

    static std::vector<std::vector<T>> subtract(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }
        return result;
    }

    static std::vector<std::vector<T>> standardMultiplication(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B) {
        int n = A.size();
        std::vector<std::vector<T>> result(n, std::vector<T>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }
};


#endif // MULTIPLICATION_POLICIES_HPP