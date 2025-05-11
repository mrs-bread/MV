#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
const double M_PI = 3.1415926535;
enum class Rule { Прямоугольник, Трапеция, Симпсон };

double ЧисленноеИнтегрирование(const function<double(double)>& f, double a, double b, int n, Rule rule)
{
    double h = (b - a) / n;
    double sum = 0.0;
    switch (rule) {
    case Rule::Прямоугольник:
        for (int i = 0; i < n; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    case Rule::Трапеция:
        sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    case Rule::Симпсон:
        if (n % 2 != 0) { ++n; h = (b - a) / n; }
        sum = f(a) + f(b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            sum += f(x) * ((i % 2) ? 4.0 : 2.0);
        }
        return sum * h / 3.0;
    }
    return 0.0;
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
        b.assign(n, 0.0); d.assign(n, 0.0);
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
            I += a[i] * h + b[i] * h * h / 2 + c[i] * h * h * h / 3 + d[i] * h * h * h * h / 4;
        }
        return I;
    }
};

int main() {
    setlocale(LC_ALL, "RUS");
    cout << fixed << setprecision(6);
    const double inv_sqrt_pi = 2.0 / sqrt(M_PI);
    auto Функция_Ошибки = [&](double x) {
        auto integrand = [&](double t) { return exp(-t * t); };
        double I = ЧисленноеИнтегрирование(integrand, 0.0, x,1000, Rule::Симпсон);
        return inv_sqrt_pi * I;
        };
    cout << " x           Симпсон      Таблица      Ошибка\n";
    for (int i = 0; i <= 20; ++i) {
        double x = 0.1 * i;
        double en = Функция_Ошибки(x), es = erf(x);
        cout << setw(4) << x << "   "
            << setw(10) << en << "   "
            << setw(10) << es << "   "
            << setw(10) << fabs(en - es) << "\n";
    }
    cout << "\n";

    auto Функция_2 = [&](double x) { return 4.0 / (1.0 + x * x); };
    cout << " n    Прямоугольник  Трапеция     Сплайн       Ошибка_Пр   Ошибка_Тр   Ошибка_Сп\n";
    for (int n : {8, 32, 128}) {
        double pi_rect = ЧисленноеИнтегрирование(Функция_2, 0.0, 1.0, n, Rule::Прямоугольник);
        double pi_trap = ЧисленноеИнтегрирование(Функция_2, 0.0, 1.0, n, Rule::Трапеция);
        vector<double> xs(n + 1), ys(n + 1);
        for (int i = 0; i <= n; ++i) {
            xs[i] = double(i) / n;
            ys[i] = Функция_2(xs[i]);
        }
        Spline sp(xs, ys);
        double pi_spline = sp.Вычисление();
        double err_rect = fabs(pi_rect - M_PI);
        double err_trap = fabs(pi_trap - M_PI);
        double err_spline = fabs(pi_spline - M_PI);

        cout << setw(3) << n << "   "
            << setw(10) << pi_rect << "   "
            << setw(10) << pi_trap << "   "
            << setw(10) << pi_spline << "   "
            << setw(10) << err_rect
            << setw(12) << err_trap << setw(12)
            << err_spline << "\n";
    }
    cout << "\n";

    auto Функция_3 = [&](double x) {
        if (x <= 2.0)       return exp(x * x);
        else return 1.0 / (4.0 - sin(16.0 * M_PI * x));
        };
    double I = ЧисленноеИнтегрирование(Функция_3, 0.0, 4.0,100000, Rule::Симпсон);
    cout << "Интеграл (задание 3) = " << setprecision(10) << I << "\n";

    return 0;
}
