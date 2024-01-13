#include <iostream>
#include "Matrix.hpp"

int main() {

    Matrix<2, 2, double> mat1;
    std::cout << "Matrix initialized with default constructor:\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << mat1(i, j) << " ";
        }
        std::cout << "\n";
    }

    double initData[2][2] = {{1.0, 2.0}, {3.0, 4.0}};
    Matrix<2, 2, double> mat2(initData);
    std::cout << "\nMatrix initialized with specific values:\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << mat2(i, j) << " ";
        }
        std::cout << "\n";
    }

    double scalar = 2.3;
    Matrix<2, 2, double> mat3 = mat2 * scalar;
    std::cout << "\nMatrix after multiplying by scalar " << scalar << ":\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << mat3(i, j) << " ";
        }
        std::cout << "\n";
    }

    auto transposedMat = mat2.transpose();
    std::cout << "\nTransposed Matrix:\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << transposedMat(i, j) << " ";
        }
        std::cout << "\n";
    }

    auto negatedMat = mat3.negate();
    std::cout << "\nNegated Matrix:\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << negatedMat(i, j) << " ";
        }
        std::cout << "\n";
    }

    std::cout << "Determinant (Laplace Expansion): " << mat2.determinant() << std::endl;

    std::cout << "Determinant (Laplace Expansion): " << mat3.determinant() << std::endl;

    auto inversedMat = mat3.inverse();
    std::cout << "\nInversed Matrix:\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << inversedMat(i, j) << " ";
        }
        std::cout << "\n";
    }

    std:: cout << "Sum of mat2 and mat3: " << mat2 + mat3 << std::endl;
    std:: cout << "Mul of mat2 and mat3: " << mat2 * mat3 << std::endl;
    std:: cout << "Elem mul of mat2 and mat3" << mat2.elementWiseMultiply(mat3) << std::endl;


    std::cout << "mat3" << mat3 << std::endl;
    std::cout << "L of mat3" << mat3.luDecomposition().first << std::endl;
    std::cout << "U of mat3" << mat3.luDecomposition().second << std::endl;
    std::cout << "Recomposition of mat3" << mat3.luDecomposition().first * mat3.luDecomposition().second << std::endl;

    std::cout << "mat3" << mat3 << std::endl;
    std::cout << "Q of mat3" << mat3.qrDecomposition().first << std::endl;
    std::cout << "QTQ of mat3" << mat3.qrDecomposition().first.transpose() * mat3.qrDecomposition().first << std::endl;
    std::cout << "R of mat3" << mat3.qrDecomposition().second << std::endl;
    std::cout << "Recomposition of mat3" << mat3.qrDecomposition().first * mat3.qrDecomposition().second << std::endl;

    double data4[4][4] = {{4, 1, 1, 1},{1, 4, 1, 1},{1, 1, 4, 1},{1, 1, 1, 4}};
    Matrix<4, 4, double> mat4(data4);
    std::cout << "mat4" << mat4 << std::endl;
    std::cout << "L of mat4" << mat4.choleskyDecomposition() << std::endl;
    std::cout << "Recomposition of mat4" << mat4.choleskyDecomposition() * mat4.choleskyDecomposition().transpose() << std::endl;

    std::cout << "Eigenvalue of mat4 " << mat4.eigenvalueDecomposition() << std::endl;

    std::cout << "Ortonormalization of mat4: " << mat4.ortogonalize() << std::endl;

    std::vector<std::vector<double>> matA = {
        {4, 1, 1, 1},
        {1, 4, 1, 1},
        {1, 1, 4, 1},
        {1, 1, 1, 4}
    };

    std::vector<double> vecB = {6, 7, 8, 9};

    Matrix<4, 4, double> A(matA);
    Matrix<4, 1, double> b(vecB);

    auto x = A.solve(b);

    std::cout << "Solution using Gaussian Elimination:" << std::endl;
    std::cout << x << std::endl;

    auto xIterative = A.solveIteratively(b, 1e-7, 1000);

    std::cout << "Solution using Gauss-Seidel method:" << std::endl;
    std::cout << xIterative << std::endl;

    auto xDecompose = A.solveWithDecompose(b);

    std::cout << "Solution using Decomposition:" << std::endl;
    std::cout << xDecompose << std::endl; 

    return 0;
}