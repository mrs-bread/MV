#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
using namespace std;
void Kvad(double a, double b, double c) {
    cout << fixed << setprecision(5);
    cout << "Решаем уравнение: " << a << "x^2 + " << b << "x + " << c << " = 0" << endl;
    constexpr double machine_zero = numeric_limits<double>::min();
    constexpr double dbl_max = numeric_limits<double>::max();
    if (abs(a) <= machine_zero) {
        cout << "Коэффициент 'a' близок к нулю"<< endl;
        if (abs(b) <= machine_zero) {
            cout << "Коэффициент 'b' также близок к нулю" << endl;
            if (abs(c) <= machine_zero) {
                cout << "Коэффициент 'c' также близок к нулю" << endl;
                cout << "Уравнение вырождается: Бесконечно много решений." << endl;
            }
            else {
                cout << "Коэффициент 'c' не близок к нулю" << endl;
                cout << "Уравнение вырождается: a = 0, b = 0, c != 0. Нет решений." << endl;
            }
        }
        else {
            cout << "Коэффициент 'b' не близок к нулю" << endl;
            double x = -c / b;
            cout << "Уравнение является линейным: bx + c = 0." << endl << " Решение: x = " << x << endl;
        }
        return;
    }
    bool over = false;
    if (abs(a)*4*abs(c) >= dbl_max) {
        cout << endl << "Предупреждение: Возможно переполнение при вычислении 4*a*c." << endl;
        over = true;
    }
    if (abs(b) > sqrt(dbl_max)) {
        cout << endl << "Предупреждение: Возможно переполнение при вычислении b*b." << endl;
        over = true;
    }
    if (over) {
        cout <<"Внимание: Возможны переполнения при вычислениях дискриминанта." << endl;
    }
    long double D = (long double)b * b - 4.0L * (long double)a * (long double)c;
    if (D < 0) {
        cout << "Дискриминант отрицательный: D = " << D << endl << "Комплексные корни:" << endl;
        double rPart = -b / (2 * a);
        double imPart = sqrt(-D) / (2 * a);
        cout << "x1 = " << rPart << " + " << imPart << "i" << endl;
        cout << "x2 = " << rPart << " - " << imPart << "i" << endl;
    }
    else {
        cout << "Дискриминант неотрицательный: D = " << D << "." << endl;
        double x1, x2;
        if (b >= 0) {
            x1 = (-b - sqrt(D)) / (2 * a);
            x2 = c / (a * x1);
        }
        else {
            x1 = (-b + sqrt(D)) / (2 * a);
            x2 = c / (a * x1);
        }
        cout << "Действительные корни:" << endl;
        cout << "x1 = " << x1 << endl;
        cout << "x2 = " << x2 << endl;
    }
}
int main() {
    setlocale(LC_ALL, "RUS");
    cout << "Пример 1: a = 1, b = -10^8, c = 1" << endl;
    Kvad(1.0, -1e8, 1.0);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 2: a = 6*10^30, b = 5*10^30, c = -4*10^30" << endl;
    Kvad(6e30, 5e30, -4e30);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 3: a = 1.0, b = -4.0, c = 3.9999999" << endl;
    Kvad(1.0, -4.0, 3.9999999);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 4: a = 2, b = 6, c = 5" << endl;
    Kvad(2.0, 6.0, 5.0);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 5: a = 10^-10, b = -10^30, c = 10^30" << endl;
    Kvad(1e-10, -1e30, 1e30);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 6: a = " << numeric_limits<double>::max() / 4 << ", b = 1, c = 1" << endl;
    Kvad(numeric_limits<double>::max() / 4, 1.0, 1.0);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 7: a = " << numeric_limits<double>::min() / 4 << ", b = 1, c = 1" << endl;
    Kvad(numeric_limits<double>::min() / 4, 1.0, 1.0);
    cout << endl << "----------------------------------------" << endl;
    cout << "Пример 8: a = " << 0 << ", b = 0, c = 0" << endl;
    Kvad(0.0, 0.0, 0.0);
    cout << endl << "----------------------------------------" << endl;
    return 0;
}
