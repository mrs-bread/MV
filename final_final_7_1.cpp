#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <functional>

using namespace std;
const double M_PI_VAL = 3.14159265358979323846;

double Y(double x, double y) {
    return (2.0 / sqrt(M_PI_VAL)) * exp(-x * x);
}

double Рунге_Кутт_4(double X, int N) {
    double h = X / N;
    double x = 0.0;
    double y = 0.0;
    for (int i = 0; i < N; ++i) {
        double r1 = Y(x, y);
        double r2 = Y(x + h / 2, y + h * r1 / 2);
        double r3 = Y(x + h / 2, y + h * r2 / 2);
        double r4 = Y(x + h, y + h * r3);
        y += (h / 6) * (r1 + 2 * r2 + 2 * r3 + r4);
        x += h;
    }
    return y;
}

double Симпсон(function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) ++n;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x) * ((i % 2) ? 4.0 : 2.0);
    }
    return sum * h / 3.0;
}

double интеграл_ошибки(double t) {
    return exp(-t * t);
}

double ошибка(double x) {
    const double inv_sqrt_pi = 2.0 / sqrt(M_PI_VAL);
    double I = Симпсон(интеграл_ошибки, 0.0, x, 1000);
    return inv_sqrt_pi * I;
}

int main() {
    setlocale(LC_ALL, "RUS");
    const int STEPS = 20;
    const int N = 1000;
    cout << fixed << setprecision(6);
    cout << " x          Р-К4        Симпсон     erf     |Р-К4 - Симпсон||Р-К4 - erf||Симпсон - erf|\n";

    auto t1 = chrono::high_resolution_clock::now();
    vector<double> рунге(STEPS + 1);
    for (int i = 0; i <= STEPS; ++i) {
        double x = 0.1 * i;
        рунге[i] = Рунге_Кутт_4(x, N);
    }
    auto t2 = chrono::high_resolution_clock::now();

    auto t3 = chrono::high_resolution_clock::now();
    vector<double> симпсон(STEPS + 1);
    for (int i = 0; i <= STEPS; ++i) {
        double x = 0.1 * i;
        симпсон[i] = ошибка(x);
    }
    auto t4 = chrono::high_resolution_clock::now();

    for (int i = 0; i <= STEPS; ++i) {
        double x = 0.1 * i;
        double stdv = erf(x);
        cout << setw(4) << x << "   "
            << setw(9) << рунге[i] << "   "
            << setw(9) << симпсон[i] << "   "
            << setw(9) << stdv << "   "
            << setw(9) << fabs(рунге[i] - симпсон[i]) << "   "
            << setw(9) << fabs(рунге[i] - stdv) << "   "
            << setw(9) << fabs(симпсон[i] - stdv) << "\n";
    }

    auto t_rk4 = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
    auto t_simp = chrono::duration_cast<chrono::microseconds>(t4 - t3).count();
    cout << "\nВремя выполнения:\n";
    cout << " Рунге-Кутт 4-го порядка (N = " << N << "):     " << t_rk4 << " мксек\n";
    cout << " Симпсон (N = 1000): " << t_simp << " мксек\n";

    return 0;
}
