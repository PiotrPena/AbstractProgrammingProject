#ifndef EIGENVALUES_POLICIES_HPP
#define EIGENVALUES_POLICIES_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>

template<typename T>
class PowerIteration {
public:
    static std::pair<T, std::vector<T>> calculate(const std::vector<std::vector<T>>& matrix, int maxIterations = 1000, T tolerance = 1e-10) {
        int n = matrix.size();
        std::vector<T> b_k(n, 1);

        T eigenvalue = 0;
        std::vector<T> b_k1;

        for (int iter = 0; iter < maxIterations; ++iter) {

            b_k1 = multiply(matrix, b_k);

            T norm = std::sqrt(std::inner_product(b_k1.begin(), b_k1.end(), b_k1.begin(), 0.0));
            std::for_each(b_k1.begin(), b_k1.end(), [norm](T& val) { val /= norm; });

            if (std::inner_product(b_k.begin(), b_k.end(), b_k1.begin(), 0.0, std::plus<>(), [](T a, T b) { return std::abs(a - b); }) < tolerance) {
                eigenvalue = norm;
                break;
            }

            b_k = b_k1;
        }

        return {eigenvalue, b_k1};
    }

private:
    static std::vector<T> multiply(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec) {
        std::vector<T> result(vec.size(), 0);
        for (size_t i = 0; i < matrix.size(); ++i)
            for (size_t j = 0; j < matrix[0].size(); ++j)
                result[i] += matrix[i][j] * vec[j];
        return result;
    }
};


#endif // EIGENVALUES_POLICIES_HPP
