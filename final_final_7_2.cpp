#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
struct Состояние {
    double r;
    double f;
};

struct РезультатСимуляции {
    bool extinct_r;
    bool extinct_f;
    double t_ext;
};

Состояние ОДУ(const Состояние& s, double a) {
    Состояние ds;
    ds.r = 2.0 * s.r - a * s.r * s.f;
    ds.f = -1.0 * s.f + a * s.r * s.f;
    return ds;
}

Состояние РК4_шаг(const Состояние& s, double a, double dt) {
    Состояние k1 = ОДУ(s, a);
    Состояние s2{ s.r + 0.5 * dt * k1.r, s.f + 0.5 * dt * k1.f };
    Состояние k2 = ОДУ(s2, a);
    Состояние s3{ s.r + 0.5 * dt * k2.r, s.f + 0.5 * dt * k2.f };
    Состояние k3 = ОДУ(s3, a);
    Состояние s4{ s.r + dt * k3.r, s.f + dt * k3.f };
    Состояние k4 = ОДУ(s4, a);

    Состояние s_new;
    s_new.r = s.r + (dt / 6.0) * (k1.r + 2 * k2.r + 2 * k3.r + k4.r);
    s_new.f = s.f + (dt / 6.0) * (k1.f + 2 * k2.f + 2 * k3.f + k4.f);
    return s_new;
}

РезультатСимуляции Симуляция(double r0, double f0, double a, double dt, double t_max) {
    Состояние s{ r0, f0 };
    double t = 0.0;
    РезультатСимуляции result{ false, false, t_max };

    while (t < t_max) {
        if (!result.extinct_r && s.r < 1.0) {
            result.extinct_r = true;
            result.t_ext = t;
        }
        if (!result.extinct_f && s.f < 1.0) {
            result.extinct_f = true;
            result.t_ext = min(result.t_ext, t);
        }
        if (result.extinct_r && result.extinct_f) break;

        s = РК4_шаг(s, a, dt);
        t += dt;
    }
    return result;
}

int main() {
    setlocale(LC_ALL, "RUS");
    const double a = 0.01;
    const double dt = 0.01;
    const double t_max = 1000.0;

    cout << "Часть А: При a = 0.01" << endl;
    cout << "r_0\tf_0\tкролики\tлисы\tвремя" << endl;
    for (double r0 : {2, 10, 100, 200, 300, 400, 500, 1000}) {
        for (double f0 : {2, 10, 100, 200, 300, 400, 500, 1000}) {
            РезультатСимуляции res = Симуляция(r0, f0, a, dt, t_max);
            cout << r0 << '\t' << f0 << '\t'
                << (res.extinct_r ? "да" : "нет") << '\t'
                << (res.extinct_f ? "да" : "нет") << '\t';
            if (res.t_ext != 1000)
                cout << res.t_ext << endl;
            else cout << "Выжили" << endl;
        }
    }
    cout << "--------------------------------------------------------------------------";
    cout << "\nЧасть Б: r_0 = 15, f_0 = 22" << endl;
    РезультатСимуляции res15 = Симуляция(15, 22, a, dt, t_max);
    cout << "Кролики вымерли: " << (res15.extinct_r ? "да" : "нет")
        << ", лисы вымерли: " << (res15.extinct_f ? "да" : "нет")
        << ", время первого вымирания: " << res15.t_ext << endl;
    cout << "--------------------------------------------------------------------------";
    cout << "\nНачальные условия, где вымирают только лисы:" << endl;
    for (int v = 30000; v <= 39750; v += 1) {
        РезультатСимуляции r = Симуляция(v, 2, a, dt, t_max);
        if (!r.extinct_r && r.extinct_f) {
            cout << "r_0 = " << v << ", f_0 = 2, лисы вымерли: да, кролики вымерли: нет, t = " << r.t_ext << endl;
        }
    }
    cout << "--------------------------------------------------------------------------";
    cout << "\nНачальные r_0 = f_0, где оба вида вымирают: " << endl;
    for (int v = 2; v <= 850; v += 1) {
        РезультатСимуляции r = Симуляция(v, v, a, dt, t_max);
        if (r.extinct_r && r.extinct_f) {
            cout << "r_0 = f_0 = " << v << ", кролики: да, лисы: да, t = " << r.t_ext << endl;
        }
    }
    cout <<"--------------------------------------------------------------------------";
    return 0;
}
