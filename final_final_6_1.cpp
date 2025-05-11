#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;
const double M_PI = 3.1415926535;

double Прямоугольник(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += f(a + i * h);
    }
    return sum * h;
}

double Трапеция(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return sum * h;
}

double Симпсон(double (*f)(double), double a, double b, int n) {
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

double Функция_Ошибки(double x) {
    const double inv_sqrt_pi = 2.0 / sqrt(M_PI);
    double I = Симпсон(интеграл_ошибки, 0.0, x, 1000);
    return inv_sqrt_pi * I;
}

double Функция_2(double x) {
    return 4.0 / (1.0 + x * x);
}

double Функция_3(double x) {
    if (x <= 2.0)
        return exp(x * x);
    else
        return 1.0 / (4.0 - sin(16.0 * M_PI * x));
}

struct Spline {
    vector<double> x, a, b, c, d;
    Spline(const vector<double>& xs, const vector<double>& ys) {
        int n = xs.size() - 1;
        x = xs; a = ys;
        c.assign(n + 1, 0.0);
        vector<double> alpha(n), beta(n);
        for (int i = 1; i < n; ++i) {
            double hi = xs[i] - xs[i - 1];
            double hip = xs[i + 1] - xs[i];
            double A = hi, C = 2 * (hi + hip), B = hip;
            double F = 3 * ((ys[i + 1] - ys[i]) / hip - (ys[i] - ys[i - 1]) / hi);
            double z = A * alpha[i - 1] + C;
            alpha[i] = -B / z;
            beta[i] = (F - A * beta[i - 1]) / z;
        }
        for (int i = n - 1; i >= 1; --i) {
            c[i] = alpha[i] * c[i + 1] + beta[i];
        }
        b.assign(n, 0.0);
        d.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double h = xs[i + 1] - xs[i];
            d[i] = (c[i + 1] - c[i]) / (3 * h);
            b[i] = (ys[i + 1] - ys[i]) / h - h * (2 * c[i] + c[i + 1]) / 3;
        }
    }
    double Вычисление() const {
        double I = 0.0;
        int n = a.size() - 1;
        for (int i = 0; i < n; ++i) {
            double h = x[i + 1] - x[i];
            I += a[i] * h
                + b[i] * h * h / 2
                + c[i] * h * h * h / 3
                + d[i] * h * h * h * h / 4;
        }
        return I;
    }
};

int main() {
    setlocale(LC_ALL, "RUS");
    cout << fixed << setprecision(6);
    cout << " x           Симпсон      Таблица      Ошибка\n";
    for (int i = 0; i <= 20; ++i) {
        double x = 0.1 * i;
        double en = Функция_Ошибки(x);
        double es = erf(x);
        cout << setw(4) << x << "   "
            << setw(10) << en << "   "
            << setw(10) << es << "   "
            << setw(10) << fabs(en - es) << "\n";
    }
    cout << "\n";

    cout << " n    Прямоугольник  Трапеция     Сплайн       Ошибка_Пр   Ошибка_Тр   Ошибка_Сп\n";
    for (int n : {8, 32, 128}) {
        double pi_rect = Прямоугольник(Функция_2, 0.0, 1.0, n);
        double pi_trap = Трапеция(Функция_2, 0.0, 1.0, n);

        vector<double> xs(n + 1), ys(n + 1);
        for (int i = 0; i <= n; ++i) {
            xs[i] = double(i) / n;
            ys[i] = Функция_2(xs[i]);
        }
        Spline sp(xs, ys);
        double pi_spline = sp.Вычисление();

        cout << setw(3) << n << "   "
            << setw(10) << pi_rect << "   "
            << setw(10) << pi_trap << "   "
            << setw(10) << pi_spline
            << "   " << setw(10) << fabs(pi_rect - M_PI)
            << setw(12) << fabs(pi_trap - M_PI)
            << setw(12) << fabs(pi_spline - M_PI)
            << "\n";
    }
    cout << "\n";

    double I3 = Симпсон(Функция_3, 0.0, 4.0, 1000000);
    cout << "Интеграл (задание 3) = " << setprecision(10) << I3 << "\n";

    return 0;
}
