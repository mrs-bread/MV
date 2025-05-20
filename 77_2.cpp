#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// ----------------------------------------------------------------------------------
// Решение системы ODE методом Рунге–Кутты 4-го порядка (RK4).
// Система:
//    dr/dt = 2r - a * r * f
//    df/dt = -f + a * r * f
//
// Параметры:
//   a   – константа взаимодействия (здесь 0.01)
//   r0  – начальное число кроликов
//   f0  – начальное число лис
//   dt  – шаг интегрирования по времени
//   T_max – максимально допустимое время (выход по времени, если не произошло вымирание)
//
// Функции:
//   deriv_r(r, f)   – правая часть d r / dt
//   deriv_f(r, f)   – правая часть d f / dt
//   rk4_step(r, f)  – один шаг RK4, возвращает пару (r_next, f_next)
// ----------------------------------------------------------------------------------

static const double a_const = 0.01; // константа a = 0.01

// Правая часть для dr/dt
double deriv_r(double r, double f) {
    return 2.0 * r - a_const * r * f;
}

// Правая часть для df/dt
double deriv_f(double r, double f) {
    return -1.0 * f + a_const * r * f;
}

// Обычный RK4-шаг для системы (r, f)
void rk4_step(double r, double f, double dt, double& r_next, double& f_next) {
    // k1
    double k1_r = deriv_r(r, f);
    double k1_f = deriv_f(r, f);

    // k2
    double r_mid = r + 0.5 * dt * k1_r;
    double f_mid = f + 0.5 * dt * k1_f;
    double k2_r = deriv_r(r_mid, f_mid);
    double k2_f = deriv_f(r_mid, f_mid);

    // k3
    r_mid = r + 0.5 * dt * k2_r;
    f_mid = f + 0.5 * dt * k2_f;
    double k3_r = deriv_r(r_mid, f_mid);
    double k3_f = deriv_f(r_mid, f_mid);

    // k4
    double r_end = r + dt * k3_r;
    double f_end = f + dt * k3_f;
    double k4_r = deriv_r(r_end, f_end);
    double k4_f = deriv_f(r_end, f_end);

    // Итоговый шаг
    r_next = r + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r);
    f_next = f + (dt / 6.0) * (k1_f + 2.0 * k2_f + 2.0 * k3_f + k4_f);
}

// ------------------------------------------------------------------------------------------------
// Функция simulate_one_pair(r0, f0, dt, T_max)
//   – имитирует эволюцию популяций (r, f) с начальных r0, f0 до момента вымирания или T_max.
//   – Возвращает кортеж: (time_extinction, final_r, final_f, status),
//     где status = 1, если вымерли кролики (r < 1), status = 2, если вымерли лисы (f < 1),
//     status = 3, если вымерли оба (в один и тот же шаг), status = 0, если ни одна
//     популяция не вымерла к T_max.
// ------------------------------------------------------------------------------------------------
struct Result {
    double t_ext;   // время, когда обнаружено вымирание (в секундах)
    double r_final; // r(t_ext) или r(T_max)
    double f_final; // f(t_ext) или f(T_max)
    int status;     // 0 – нет вымирания, 1 – кролики, 2 – лисы, 3 – оба сразу
};

Result simulate_one_pair(double r0, double f0, double dt, double T_max) {
    double r = r0;
    double f = f0;
    double t = 0.0;

    // Если сразу одно из значений < 1, считаем мгновенным вымиранием
    if (r < 1.0 && f < 1.0) {
        return { 0.0, r, f, 3 };
    }
    if (r < 1.0) {
        return { 0.0, r, f, 1 };
    }
    if (f < 1.0) {
        return { 0.0, r, f, 2 };
    }

    while (t < T_max) {
        double r_next, f_next;
        rk4_step(r, f, dt, r_next, f_next);
        t += dt;
        r = r_next;
        f = f_next;

        bool rabbit_dead = (r < 1.0);
        bool fox_dead = (f < 1.0);

        if (rabbit_dead && fox_dead) {
            return { t, r, f, 3 };
        }
        else if (rabbit_dead) {
            return { t, r, f, 1 };
        }
        else if (fox_dead) {
            return { t, r, f, 2 };
        }
        // Иначе продолжаем
    }

    // Если дошли до T_max без вымирания
    return { T_max, r, f, 0 };
}

int main() {
    setlocale(LC_ALL, "RUS");
    cout << fixed << setprecision(6);
    
    // -----------------------
    // Часть Б. Конкретный случай r0 = 15, f0 = 22
    // -----------------------
    {
        double r0 = 15.0;
        double f0 = 22.0;
        double dt = 0.01;     // шаг интегрирования по времени
        double T_max = 1000;  // максимально допустимое время (довольно большое, чтобы увидеть динамику)

        Result res = simulate_one_pair(r0, f0, dt, T_max);
        cout << "=== Симуляция для r0 = " << r0 << ", f0 = " << f0 << " ===\n";
        if (res.status == 0) {
            cout << "  Кролики и лисы не вымерли к t = " << res.t_ext
                << ". r(" << res.t_ext << ") = " << res.r_final
                << ", f(" << res.t_ext << ") = " << res.f_final << "\n";
        }
        else {
            const char* who =
                (res.status == 1 ? "Кролики" :
                    res.status == 2 ? "Лисы" :
                    "Оба вида");
            cout << "  Время вымирания (" << who << "): t = " << res.t_ext << "\n";
            cout << "  Значения в момент вымирания: r = " << res.r_final
                << ", f = " << res.f_final << "\n";
        }
        cout << endl;
    }

    // -----------------------
    // Часть А. Исследование поведения при a = 0.01 и разных r0, f0 от 2–3 до нескольких тысяч.
    // Можно снять профили популяций или просто определить, что происходит (рост/вымирание).
    // Ниже приведён пример: перебираем r0 и f0 с шагом 50 от 50 до 300.
    // Для каждого (r0,f0) выводим, что вымерло первым (или ни одно).
    // -----------------------
    {
        cout << "=== Перебор начальных условий (r0, f0) для a=0.01 ===\n";
        double dt = 0.01;
        double T_max = 1000;  // ограничим время, чтобы программа не «висела» слишком долго

        cout << " r0 \\ f0 | ";
        for (int f0 = 1; f0 <= 30; f0 += 1) {
            cout << setw(8) << f0;
        }
        cout << "\n";
        cout << string(10 + 8 * ((300 - 50) / 50 + 1), '-') << "\n";

        for (int r0 = 1; r0 <= 30; r0 += 1) {
            cout << setw(8) << r0 << " | ";
            for (int f0 = 1; f0 <= 30; f0 += 1) {
                Result res = simulate_one_pair(r0, f0, dt, T_max);
                char c;
                if (res.status == 0)       c = 'N'; // No extinctions (ни один вид не вымер)
                else if (res.status == 1)  c = 'R'; // Rabbits died first (кролики)
                else if (res.status == 2)  c = 'F'; // Foxes died first (лисы)
                else                       c = 'B'; // Both died simultaneously
                cout << setw(8) << c;
            }
            cout << "\n";
        }
        cout << "(N = ни один, R = кролики, F = лисы, B = оба)\n\n";
    }

    // -----------------------
    // Часть Б (дополнительно). Найти:
    //  - начальные условия, при которых кролики вымирают (r < 1)
    //  - начальные условия, при которых лисы вымирают (f < 1)
    //  - начальные условия с r0 = f0, при которых вымирают оба
    //
    // Ниже пример: перебираем r0 = f0 от 2 до 500 и ищем случай, когда оба умирают.
    // -----------------------
    {
        cout << "=== Пары (r0 = f0), где оба вида умирают одновременно ===\n";
        double dt = 0.01;
        double T_max = 2000;

        vector<int> both_die;
        for (int x = 2; x <= 5000; ++x) {
            Result res = simulate_one_pair((double)x, (double)x, dt, T_max);
            if (res.status == 3) {
                both_die.push_back(x);
            }
        }

        if (!both_die.empty()) {
            for (int x : both_die) {
                cout << "  r0 = f0 = " << x << "\n";
            }
        }
        else {
            cout << "  Ни одного x ∈ [2..300], где оба вида умирают одновременно, не найдено.\n";
        }
        cout << endl;
    }



    // -----------------------
    // Дополнительно. Найти пару (r0, f0), при которой вымирают исключительно кролики,
    // и пару, при которой вымирают исключительно лисы. Например, пробегаем r0 от 2..500,
    // f0 от 2..500, и выводим первые попавшиеся.
    // -----------------------
    {
        cout << "=== Все пары (r0, f0) для вымирания только одного вида ===\n";
        double dt = 0.01;
        double T_max = 500;

        vector<pair<int, int>> rabbits_die;
        vector<pair<int, int>> foxes_die;

        for (int r0 = 2; r0 <= 10; ++r0) {
            for (int f0 = 2; f0 <= 10; ++f0) {
                Result res = simulate_one_pair((double)r0, (double)f0, dt, T_max);
                if (res.status == 1) {
                    rabbits_die.emplace_back(r0, f0);
                }
                if (res.status == 2) {
                    foxes_die.emplace_back(r0, f0);
                }
            }
        }

        if (!foxes_die.empty()) {
            cout << "-- Пары, где лисы умирают первыми: (r0, f0) ---\n";
            for (auto& p : foxes_die) {
                cout << "  (" << p.first << ", " << p.second << ")\n";
            }
        }
        else {
            cout << "-- Ни одной пары, где лисы умирают первыми, не найдено.\n";
        }

        if (!rabbits_die.empty()) {
            cout << "-- Пары, где кролики умирают первыми: (r0, f0) ---\n";
            for (auto& p : rabbits_die) {
                cout << "  (" << p.first << ", " << p.second << ")\n";
            }
        }
        else {
            cout << "-- Ни одной пары, где кролики умирают первыми, не найдено.\n";
        }
    }
    return 0;
}
