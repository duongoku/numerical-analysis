#include <bits/stdc++.h>

using namespace std;

namespace single {
float fast_pow(float base, int exp) {
    if (exp == 0) return 1;
    if (exp == 1) return base;
    float result = fast_pow(base, exp / 2);
    result *= result;
    if (exp % 2 == 1) result *= base;
    return result;
}

double polynomial(vector<double> polynomial, double x) {
    double result = 0;
    double current = 1;
    for (int i = 0; i < polynomial.size(); i++) {
        result += polynomial[i] * current;
        current *= x;
    }
    return result;
}

vector<double> bisection(double (*f)(double), double lower, double upper,
                         double tolerance, int max_iter = 1000) {
    vector<double> result;
    for (int iter = 1; iter < max_iter; iter++) {
        double mid = (upper + lower) / 2;
        double value = f(mid);
        result.push_back(mid);
        if (abs((upper - lower) / 2) <= tolerance || value == 0) {
            return result;
        }
        if (value > 0) {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    result.push_back((upper + lower) / 2);
    return result;
}

vector<double> fixed_point(double (*g)(double), double x0, double tolerance,
                           int max_iter = 1000) {
    vector<double> result;
    double x = g(x0);
    for (int iter = 1; iter < max_iter; iter++) {
        result.push_back(x);
        if (abs(x - x0) <= tolerance) {
            return result;
        }
        x0 = x;
        x = g(x0);
    }
    result.push_back(x);
    return result;
}

vector<double> newtons(double (*f)(double), double (*f_prime)(double),
                       double x0, double tolerance, int max_iter = 1000) {
    vector<double> result;
    double x = x0 - f(x0) / f_prime(x0);
    for (int iter = 1; iter < max_iter; iter++) {
        result.push_back(x);
        if (abs(x - x0) <= tolerance) {
            return result;
        }
        x0 = x;
        x = x0 - f(x0) / f_prime(x0);
    }
    result.push_back(x);
    return result;
}

vector<double> secant(double (*f)(double), double x0, double x1,
                      double tolerance, int max_iter = 1000) {
    vector<double> result;
    double x = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    for (int iter = 1; iter < max_iter; iter++) {
        result.push_back(x);
        if (abs(x - x1) <= tolerance) {
            return result;
        }
        x0 = x1;
        x1 = x;
        x = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    }
    result.push_back(x);
    return result;
}

vector<double> false_position(double (*f)(double), double x0, double x1,
                              double tolerance, int max_iter = 1000) {
    vector<double> result;
    double x = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    for (int iter = 1; iter < max_iter; iter++) {
        result.push_back(x);
        if (abs(x - x1) <= tolerance || f(x) == 0) {
            return result;
        }
        if (f(x) * f(x1) < 0) {
            x0 = x1;
        }
        x1 = x;
        x = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    }
    result.push_back(x);
    return result;
}

}  // namespace single