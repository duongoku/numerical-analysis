#include <bits/stdc++.h>

#include "matrix.cpp"
#include "single.cpp"

using namespace std;

using namespace matrix_general;

#define PI 3.14159265358979323846

#define TAB setw(10)

int random() { return rand() % 10000 + 1; }

matrix generateSymmetric(int n) {}

pair<matrix, matrix> main_read_file() {
    freopen("input.txt", "r", stdin);
    int n;
    cin >> n;
    matrix a(n, row(n));
    matrix b(n, row(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> a[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> b[i][j];
        }
    }
    return make_pair(a, b);
}

pair<matrix, matrix> main_random(int n) {
    matrix a(n, row(n));
    matrix b(n, row(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = random();
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = random();
        }
    }
    return make_pair(a, b);
}

// f1(x) = x^3 - 4x^2 - 10
double f1(double x) { return x * x * x + 4 * x * x - 10; }

// g1(x) = 0.5(10 - x^3)^0.5
double g1(double x) { return 0.5 * sqrt(10 - x * x * x); }

// f2(x) = cos(x) − x
double f2(double x) { return cos(x) - x; }

// f2_prime(x) = -sin(x) - 1
double f2_prime(double x) { return -sin(x) - 1; }

void test_bisection() {
    cout << "Bisection Method" << endl;
    vector<double> result = single::bisection(f1, 1, 2, 1e-4);
    for (int i = 0; i < result.size(); i++) {
        cout << "Iteration " << i + 1 << ": " << fixed << setprecision(9)
             << result[i] << endl;
    }
    cout << endl;
}

void test_fixed_point() {
    cout << "Fixed Point Iteration" << endl;
    vector<double> result = single::fixed_point(g1, 1.5, 0);
    int interval = 13;
    for (int i = 0; i < interval; i++) {
        int temp = result.size() / interval;
        for (int j = 0; j < temp; j++) {
            int k = j * interval + i;
            if (k >= result.size()) {
                break;
            }
            cout << setw(12) << "Iteration " << k + 1 << ": " << fixed
                 << setprecision(9) << result[k];
        }
        cout << endl;
    }
    cout << endl;
}

void test_newtons() {
    cout << "Newton's Method" << endl;
    vector<double> result = single::newtons(f2, f2_prime, PI / 4, 0);
    cout << "p0: " << fixed << setprecision(10) << PI / 4 << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << "p" << i + 1 << ": " << fixed << setprecision(10) << result[i]
             << endl;
    }
    cout << endl;
}

void test_secant() {
    cout << "Secant Method" << endl;
    vector<double> result = single::secant(f2, 0.5, PI / 4, 0);
    cout << "p0: " << fixed << setprecision(10) << 0.5 << endl;
    cout << "p1: " << fixed << setprecision(10) << PI / 4 << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << "p" << i + 2 << ": " << fixed << setprecision(10) << result[i]
             << endl;
    }
    cout << endl;
}

void test_false_position() {
    cout << "False Position Method" << endl;
    vector<double> result = single::false_position(f2, 0.5, PI / 4, 0);
    cout << "p0: " << fixed << setprecision(10) << 0.5 << endl;
    cout << "p1: " << fixed << setprecision(10) << PI / 4 << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << "p" << i + 2 << ": " << fixed << setprecision(10) << result[i]
             << endl;
    }
    cout << endl;
}

void test_gauss_elimination() {
    cout << "Gauss Elimination" << endl;
    matrix A = {{1, 1, 0, 3}, {2, 1, -1, 1}, {3, -1, -1, 2}, {-1, 2, 3, -1}};
    matrix B = {{4}, {1}, {-3}, {4}};
    matrix X = matrix_general::gauss_elimination(A, B);
    matrix C = matrix_general::multiply(A, X);
    if (matrix_general::compare(B, C)) {
        cout << "Correct" << endl;
    } else {
        cout << "Incorrect" << endl;
    }
    cout << "X: " << endl;
    print_matrix(X);
}

void test_gauss_jordan() {
    cout << "Gauss Jordan" << endl;
    matrix A = {{1, 1, 0, 3}, {2, 1, -1, 1}, {3, -1, -1, 2}, {-1, 2, 3, -1}};
    matrix B = {{4}, {1}, {-3}, {4}};
    matrix X = matrix_general::gauss_jordan(A, B);
    matrix C = matrix_general::multiply(A, X);
    if (matrix_general::compare(B, C)) {
        cout << "Correct" << endl;
    } else {
        cout << "Incorrect" << endl;
    }
    cout << "X: " << endl;
    print_matrix(X);
}

void test_doolittle() {
    pair<matrix, matrix> data = main_random(4);
    cout << "Random matrix A:" << endl;
    print_matrix(data.first);
    matrix a = data.first;
    matrix b = data.second;
    pair<matrix, matrix> doo = doolittle(a);
    cout << "L:" << endl;
    print_matrix(doo.first);
    cout << "U:" << endl;
    print_matrix(doo.second);
    cout << "L * U:\n";
    print_matrix(matrix_general::multiply(doo.first, doo.second));
    cout << "Doolittle: "
         << (verifyDecomp(a, doo.first, doo.second) ? "Correct" : "Incorrect")
         << endl
         << endl;
}

void test_cholesky() {
    matrix a = randomPositiveDefiniteMatrix(4);
    cout << "Random positive definite matrix A:" << endl;
    print_matrix(a);
    matrix l = cholesky(a);
    matrix u = transpose(l);
    cout << "L:" << endl;
    print_matrix(l);
    cout << "Verify Cholesky: "
         << (verifyDecomp(a, l, u) ? "Correct" : "Incorrect") << endl;
}

void test_jacobi() {
    cout << "Jacobi" << endl;
    matrix a = {
        {10, -1, 2, 0}, {-1, 11, -1, 3}, {2, -1, 10, -1}, {0, 3, -1, 8}};
    matrix b = {{6}, {25}, {-11}, {15}};
    matrix x = {{0}, {0}, {0}, {0}};
    vector<matrix> result = jacobi(a, b, x, 1e-3);
    result.insert(result.begin(), x);
    cout << " " << TAB << "k" << TAB;
    for (int i = 0; i < result.size(); i++) {
        cout << i << TAB;
    }
    cout << endl;
    for (int i = 0; i < x.size(); i++) {
        cout << "x" << i + 1 << TAB;
        for (int j = 0; j < result.size(); j++) {
            cout << fixed << setprecision(5) << result[j][i][0] << TAB;
        }
        cout << endl;
    }
    cout << "------------" << endl;
}

void test_gauss_seidel() {
    cout << "Gauss-Seidel" << endl;
    matrix a = {
        {10, -1, 2, 0}, {-1, 11, -1, 3}, {2, -1, 10, -1}, {0, 3, -1, 8}};
    matrix b = {{6}, {25}, {-11}, {15}};
    matrix x = {{0}, {0}, {0}, {0}};
    vector<matrix> result = gauss_seidel(a, b, x, 1e-3);
    result.insert(result.begin(), x);
    cout << " " << TAB << "k" << TAB;
    for (int i = 0; i < result.size(); i++) {
        cout << i << TAB;
    }
    cout << endl;
    for (int i = 0; i < x.size(); i++) {
        cout << "x" << i + 1 << TAB;
        for (int j = 0; j < result.size(); j++) {
            cout << fixed << setprecision(5) << result[j][i][0] << TAB;
        }
        cout << endl;
    }
    cout << "------------" << endl;
}

int main() {
    srand(2);

    const int TEST = 4;

    switch (TEST) {
        case 1:
            test_gauss_seidel();
            test_jacobi();
            break;
        case 2:
            test_doolittle();
            break;
        case 3:
            test_cholesky();
            break;
        case 4:
            test_bisection();
            test_fixed_point();
            test_newtons();
            test_secant();
            test_false_position();
            break;
        case 5:
            test_gauss_elimination();
            test_gauss_jordan();
            break;
    }

    // test_gauss_seidel();
    // test_jacobi();
    // return 0;

    // test_doolittle();
    // return 0;

    // test_cholesky();
    // return 0;

    // test_bisection();
    // test_fixed_point();
    // test_newtons();
    // test_secant();
    // test_false_position();
    // test_gauss_elimination();
    // test_gauss_jordan();
    // return 0;
}
