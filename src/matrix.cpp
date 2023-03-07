#include <bits/stdc++.h>

using namespace std;

namespace matrix_general {

typedef vector<vector<double>> matrix;
typedef vector<double> row;

#define NULL_MATRIX matrix(0, row(0))

void print_matrix(matrix a) {
    int n = a.size();
    for (int i = 0; i < n; i++) {
        cout << "[";
        int m = a[i].size();
        for (int j = 0; j < m; j++) {
            cout << setw(8) << a[i][j] << (j == m - 1 ? "" : ", ");
        }
        cout << "]" << endl;
    }
    cout << endl;
}

bool compare(matrix a, matrix b, double epsilon = 0.0001) {
    if (a.size() != b.size()) return false;
    int n = a.size();
    for (int i = 0; i < n; i++) {
        if (a[i].size() != b[i].size()) return false;
        int m = a[i].size();
        for (int j = 0; j < m; j++) {
            if (abs(a[i][j] - b[i][j]) > epsilon) return false;
        }
    }
    return true;
}

double norm(matrix a) {
    double result = 0;
    int n = a.size();
    for (int i = 0; i < n; i++) {
        int m = a[i].size();
        for (int j = 0; j < m; j++) {
            result += a[i][j] * a[i][j];
        }
    }
    return sqrt(result);
}

matrix transpose(matrix a) {
    int size = a.size();
    matrix result(size, row(size));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = a[j][i];
        }
    }
    return result;
}

// Square matrices
matrix add(matrix a, matrix b, int size) {
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            a[i][j] += b[i][j];
        }
    }
    return a;
}

// Square matrices
matrix subtract(matrix a, matrix b) {
    int n = a.size();
    int m = a[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            a[i][j] -= b[i][j];
        }
    }
    return a;
}

bool checkSymmetric(matrix a) {
    int size = a.size();
    if (size != a[0].size()) return false;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (a[i][j] != a[j][i]) return false;
        }
    }
    return true;
}

bool checkSquare(matrix a) {
    int size = a.size();
    return size == a[0].size();
}

bool matrixEqual(matrix a, matrix b) {
    if (a.size() != b.size() || a[0].size() != b[0].size()) return false;

    int size = a.size();
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (a[i][j] != b[i][j]) return false;
        }
    }

    return true;
}

// Divide and conquer matrix multiplication
matrix dndMultiply(matrix n, matrix m, int n_x, int n_y, int m_x, int m_y,
                   int size) {
    if (size == 1) {
        matrix result(1, row(1));
        result[0][0] = n[n_y][n_x] * m[m_y][m_x];
        return result;
    }

    // hs = half size
    int hs = size / 2;
    matrix a = add(dndMultiply(n, m, n_x, n_y, m_x, m_y, hs),
                   dndMultiply(n, m, n_x + hs, n_y, m_x, m_y + hs, hs), hs);
    matrix b =
        add(dndMultiply(n, m, n_x, n_y, m_x + hs, m_y, hs),
            dndMultiply(n, m, n_x + hs, n_y, m_x + hs, m_y + hs, hs), hs);
    matrix c =
        add(dndMultiply(n, m, n_x, n_y + hs, m_x, m_y, hs),
            dndMultiply(n, m, n_x + hs, n_y + hs, m_x, m_y + hs, hs), hs);
    matrix d =
        add(dndMultiply(n, m, n_x, n_y + hs, m_x + hs, m_y, hs),
            dndMultiply(n, m, n_x + hs, n_y + hs, m_x + hs, m_y + hs, hs), hs);

    matrix result(size, row(size));
    for (int i = 0; i < hs; i++) {
        for (int j = 0; j < hs; j++) {
            result[i][j] = a[i][j];
            result[i][j + hs] = b[i][j];
            result[i + hs][j] = c[i][j];
            result[i + hs][j + hs] = d[i][j];
        }
    }
    return result;
}

matrix multiply(matrix a, matrix b) {
    int n = a.size();
    int m = b[0].size();
    int p = b.size();
    matrix result(n, row(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

// c = ab
void verifyWithNaiveMethod(matrix a, matrix b, int size) {
    auto naive_start = chrono::high_resolution_clock::now();
    matrix c = multiply(a, b);
    auto naive_end = chrono::high_resolution_clock::now();
    auto naive_duration =
        chrono::duration_cast<chrono::microseconds>(naive_end - naive_start);
    cout << "Naive execution time: " << naive_duration.count() / 1000.0
         << " milliseconds";
    cout << endl;

    auto dnd_start = chrono::high_resolution_clock::now();
    matrix d = dndMultiply(a, b, 0, 0, 0, 0, size);
    auto dnd_end = chrono::high_resolution_clock::now();
    auto dnd_duration =
        chrono::duration_cast<chrono::microseconds>(dnd_end - dnd_start);
    cout << "DND execution time: " << dnd_duration.count() / 1000.0
         << " milliseconds";
    cout << endl;

    bool result = matrixEqual(c, d);

    if (result) {
        cout << "Correct" << endl;
    } else {
        cout << "Incorrect" << endl;
    }
}

// Matrix a must be a symmetric matrix
matrix cholesky(matrix a) {
    if (!checkSymmetric(a)) {
        return NULL_MATRIX;
    }
    int n = a.size();
    matrix l(n, row(n));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            int temp = 0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    temp += l[j][k] * l[j][k];
                }
                l[j][j] = sqrt(a[j][j] - temp);
            } else {
                for (int k = 0; k < j; k++) {
                    temp += (l[i][k] * l[j][k]);
                }
                l[i][j] = (a[i][j] - temp) / l[j][j];
            }
        }
    }

    return l;
}

// Generate a random positive definite matrix of size nxn
matrix randomPositiveDefiniteMatrix(int n) {
    matrix a(n, row(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = rand() % 100;
        }
    }
    matrix b = transpose(a);
    matrix c = multiply(a, b);
    return c;
}

// Returns pair<l, u>
pair<matrix, matrix> doolittle(matrix a) {
    if (!checkSquare(a)) {
        return make_pair(NULL_MATRIX, NULL_MATRIX);
    }

    int n = a.size();
    matrix l(n, row(n));
    matrix u(n, row(n));

    for (int i = 0; i < n; i++) l[i][i] = 1;

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            u[i][j] = a[i][j];
            for (int k = 0; k < i; k++) {
                u[i][j] -= l[i][k] * u[k][j];
            }
        }

        for (int j = i + 1; j < n; j++) {
            l[j][i] = a[j][i];
            for (int k = 0; k < i; k++) {
                l[j][i] -= l[j][k] * u[k][i];
            }
            l[j][i] /= u[i][i];
        }
    }

    return make_pair(l, u);
}

bool verifyDecomp(matrix a, matrix l, matrix u) {
    matrix lu = multiply(l, u);
    return matrixEqual(a, lu);
}

matrix gauss_elimination(matrix A, matrix B) {
    int n = A.size();
    matrix a = A;
    matrix b = B;

    for (int i = 0; i < n; i++) {
        // Find max element in current column
        int maxElement = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(a[k][i]) > abs(a[maxElement][i])) {
                maxElement = k;
            }
        }

        // Swap found element with current row
        for (int k = i; k < n; k++) {
            double tmp = a[maxElement][k];
            a[maxElement][k] = a[i][k];
            a[i][k] = tmp;
        }
        double tmp = b[maxElement][0];
        b[maxElement][0] = b[i][0];
        b[i][0] = tmp;

        // Eliminate all elements below current row
        for (int k = i + 1; k < n; k++) {
            double c = -a[k][i] / a[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    a[k][j] = 0;
                } else {
                    a[k][j] += c * a[i][j];
                }
            }
            b[k][0] += c * b[i][0];
        }
    }

    // Solve
    matrix x(n, row(1));
    for (int i = n - 1; i >= 0; i--) {
        x[i][0] = b[i][0] / a[i][i];
        for (int k = i - 1; k >= 0; k--) {
            b[k][0] -= a[k][i] * x[i][0];
        }
    }
    return x;
}

matrix gauss_jordan(matrix A, matrix B) {
    int n = A.size();
    matrix a = A;
    matrix b = B;

    for (int i = 0; i < n; i++) {
        // Find max element in current column
        int maxElement = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(a[k][i]) > abs(a[maxElement][i])) {
                maxElement = k;
            }
        }

        // Swap found element with current row
        for (int k = i; k < n; k++) {
            double tmp = a[maxElement][k];
            a[maxElement][k] = a[i][k];
            a[i][k] = tmp;
        }
        double tmp = b[maxElement][0];
        b[maxElement][0] = b[i][0];
        b[i][0] = tmp;

        // Eliminate all elements below current row
        for (int k = i + 1; k < n; k++) {
            double c = -a[k][i] / a[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    a[k][j] = 0;
                } else {
                    a[k][j] += c * a[i][j];
                }
            }
            b[k][0] += c * b[i][0];
        }
    }

    // Eliminate all elements above current row
    for (int i = n - 1; i >= 0; i--) {
        for (int k = i - 1; k >= 0; k--) {
            double c = -a[k][i] / a[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    a[k][j] = 0;
                } else {
                    a[k][j] += c * a[i][j];
                }
            }
            b[k][0] += c * b[i][0];
        }
    }

    // Solve
    matrix x(n, row(1));
    for (int i = 0; i < n; i++) {
        x[i][0] = b[i][0] / a[i][i];
    }
    return x;
}

vector<matrix> jacobi(matrix A, matrix B, matrix x, double tolerance,
                      int max_iter = 1000) {
    vector<matrix> result;

    int n = A.size();
    matrix a = A;
    matrix b = B;

    // Normalize diagonal to 1
    for (int i = 0; i < n; i++) {
        double c = a[i][i];
        for (int j = 0; j < n; j++) {
            a[i][j] /= c;
        }
        b[i][0] /= c;
    }

    for (int i = 0; i < max_iter; i++) {
        matrix x_new(n, row(1));
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++) {
                if (k != j) {
                    sum += a[j][k] * x[k][0];
                }
            }
            x_new[j][0] = b[j][0] - sum;
        }
        result.push_back(x_new);
        if (norm(subtract(x_new, x)) < tolerance) {
            return result;
        }
        x = x_new;
    }
    result.push_back(x);
    return result;
}

vector<matrix> gauss_seidel(matrix A, matrix B, matrix x, double tolerance,
                            int max_iter = 1000) {
    vector<matrix> result;

    int n = A.size();
    matrix a = A;
    matrix b = B;

    // Normalize diagonal to 1
    for (int i = 0; i < n; i++) {
        double c = a[i][i];
        for (int j = 0; j < n; j++) {
            a[i][j] /= c;
        }
        b[i][0] /= c;
    }

    for (int i = 0; i < max_iter; i++) {
        matrix x_new(n, row(1));
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum += a[j][k] * x_new[k][0];
            }
            for (int k = j + 1; k < n; k++) {
                sum += a[j][k] * x[k][0];
            }
            x_new[j][0] = b[j][0] - sum;
        }
        result.push_back(x_new);
        if (norm(subtract(x_new, x)) < tolerance) {
            return result;
        }
        x = x_new;
    }
    result.push_back(x);
    return result;
}

}  // namespace matrix_general