#ifndef CONCEPTS_HPP
#define CONCEPTS_HPP

#include <concepts>

template<typename T>
concept Negatable = requires(T a) {
    { -a } -> std::convertible_to<T>;
};

template<typename T>
concept Multiplicable = requires(T a, T b) {
    { a * b } -> std::convertible_to<T>;
};

template<typename T>
concept Addable = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
};

template<typename T>
concept Streamable = requires(std::ostream& os, const T& t) {
    { os << t } -> std::convertible_to<std::ostream&>;
};

template<typename T>
concept Arithmetic = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};

template<int M, int N, typename T>
concept SquareMatrix = (M == N);

template<typename T>
concept ComparableWithTolerance = requires(T a, T b, T tolerance) {
    { a < b } -> std::convertible_to<bool>;
    { a > b } -> std::convertible_to<bool>;
    { tolerance >= 0 } -> std::convertible_to<bool>;
};

#endif // CONCEPTS_HPP
