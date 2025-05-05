#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
const double M_PI = 3.1415926535;
enum class Rule { Rectangle, Trapezoid, Simpson };

// --- Общая composite–интеграция (для части 1 и 2A) ---
double integrateComposite(const function<double(double)>& f,
    double a, double b, int n, Rule rule)
{
    double h = (b - a) / n;
    double sum = 0.0;
    switch (rule) {
    case Rule::Rectangle:
        for (int i = 0; i < n; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    case Rule::Trapezoid:
        sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    case Rule::Simpson:
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

// --- Адаптивный метод Симпсона (для части 1 и 3) ---
double adaptiveSimpsonAux(const function<double(double)>& f,
    double a, double b, double eps, double S, double fa, double fb, double fm, int depth)
{
    double m = 0.5 * (a + b), h = 0.5 * (b - a);
    double f1 = f(a + 0.5 * h), f2 = f(b - 0.5 * h);
    double Sleft = (fa + 4 * f1 + fm) * (h / 6);
    double Sright = (fm + 4 * f2 + fb) * (h / 6);
    double S2 = Sleft + Sright;
    if (depth <= 0 || fabs(S2 - S) < 15 * eps)
        return S2 + (S2 - S) / 15.0;
    return adaptiveSimpsonAux(f, a, m, eps / 2, Sleft, fa, fm, f1, depth - 1)
        + adaptiveSimpsonAux(f, m, b, eps / 2, Sright, fm, fb, f2, depth - 1);
}

double integrateAdaptiveSimpson(const function<double(double)>& f,
    double a, double b, double eps = 1e-8)
{
    double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
    double S = (fa + 4 * fm + fb) * (b - a) / 6.0;
    return adaptiveSimpsonAux(f, a, b, eps, S, fa, fb, fm, 20);
}

// --- Кубический сплайн и сплайн-квадратура (для части 2B) ---
struct Spline {
    vector<double> x, a, b, c, d;
    Spline(const vector<double>& xs, const vector<double>& ys) {
        int n = xs.size() - 1;
        x = xs; a = ys;
        c.assign(n + 1, 0.0);
        vector<double> alpha(n), beta(n);
        // Натуральный сплайн: c[0]=c[n]=0
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
    double integrate() const {
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
    cout << fixed << setprecision(6);

    // 1) Табуляция erf(x) через адаптивный Симпсон
    const double inv_sqrt_pi = 2.0 / sqrt(M_PI);
    auto erf_num = [&](double x) {
        auto integrand = [&](double t) { return exp(-t * t); };
        double I = integrateAdaptiveSimpson(integrand, 0.0, x, 1e-8);
        return inv_sqrt_pi * I;
        };
    cout << " x      erf_num     erf_std     error\n";
    for (int i = 0; i <= 20; ++i) {
        double x = 0.1 * i;
        double en = erf_num(x), es = erf(x);
        cout << setw(4) << x << "   "
            << setw(10) << en << "   "
            << setw(10) << es << "   "
            << setw(10) << fabs(en - es) << "\n";
    }
    cout << "\n";

    // 2) Приближение pi: прямоуг., трапеций и сплайн-квадратура
    auto g = [&](double x) { return 4.0 / (1.0 + x * x); };
    cout << " n    Rect        Trap        Spline      Err_Trap    Err_Spline\n";
    for (int n : {8, 32, 128}) {
        double pi_rect = integrateComposite(g, 0.0, 1.0, n, Rule::Rectangle);
        double pi_trap = integrateComposite(g, 0.0, 1.0, n, Rule::Trapezoid);
        vector<double> xs(n + 1), ys(n + 1);
        for (int i = 0; i <= n; ++i) {
            xs[i] = double(i) / n;
            ys[i] = g(xs[i]);
        }
        Spline sp(xs, ys);
        double pi_spline = sp.integrate();
        double err_trap = fabs(pi_trap - M_PI);
        double err_spline = fabs(pi_spline - M_PI);

        cout << setw(3) << n << "   "
            << setw(10) << pi_rect << "   "
            << setw(10) << pi_trap << "   "
            << setw(10) << pi_spline << "   "
            << setw(10) << err_trap << setw(12) << err_spline << "\n";
    }
    cout << "\n";

    // 3) Кусочно-заданная f(x) и адаптивный Симпсон
    auto f = [&](double x) {
        if (x <= 2.0)       return exp(x * x);
        else /*2 < x <= 4*/ return 1.0 / (4.0 - sin(16.0 * M_PI * x));
        };
    double I = integrateAdaptiveSimpson(f, 0.0, 4.0, 1e-8);
    cout << "Integral of f(x) over [0,4] = " << setprecision(10) << I << "\n";

    return 0;
}
