#include <iostream>
#include <vector>
#include <queue>
#include <utility>   // для pair, если понадобится
#include <algorithm> // для max, swap

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;

    // Читаем матрицу
    vector<vector<int>> A(n, vector<int>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> A[i][j];

    // Построим список смежности и массив степеней
    vector<vector<int>> adj(n);
    vector<int> deg(n, 0);
    for (int i = 0; i < n; i++) {
        if (A[i][i] != 0) {
            cout << "NO\n";
            return 0;
        }
        for (int j = 0; j < n; j++) {
            if (A[i][j]) {
                if (A[j][i] == 0) {
                    cout << "NO\n";
                    return 0;
                }
                adj[i].push_back(j);
                deg[i]++;
            }
        }
    }

    // Проверка связности (BFS от 0)
    vector<char> vis(n, 0);
    queue<int> q;
    vis[0] = 1;
    q.push(0);
    int seen = 1;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (!vis[v]) {
                vis[v] = 1;
                q.push(v);
                seen++;
            }
        }
    }
    if (seen != n) {
        cout << "NO\n";
        return 0;
    }

    // Собираем угловые вершины (deg == 2)
    vector<int> corner;
    for (int i = 0; i < n; i++)
        if (deg[i] == 2)
            corner.push_back(i);

    if (corner.size() != 4) {
        cout << "NO\n";
        return 0;
    }

    // Функция BFS для расчёта расстояний
    auto bfs = [&](int start) {
        const int INF = 1e9;
        vector<int> dist(n, INF);
        queue<int> qq;
        dist[start] = 0;
        qq.push(start);
        while (!qq.empty()) {
            int u = qq.front(); qq.pop();
            for (int v : adj[u]) {
                if (dist[v] == INF) {
                    dist[v] = dist[u] + 1;
                    qq.push(v);
                }
            }
        }
        return dist;
    };

    // От первого угла считаем все расстояния
    vector<int> d0 = bfs(corner[0]);

    // Находим противоположный угол (макс. расстояние)
    int opp = corner[1];
    int maxd = d0[opp];
    for (int i = 2; i < 4; i++) {
        if (d0[corner[i]] > maxd) {
            maxd = d0[corner[i]];
            opp = corner[i];
        }
    }

    // Два оставшихся — на сторонах
    vector<int> side;
    for (int v : corner)
        if (v != corner[0] && v != opp)
            side.push_back(v);

    int d1 = d0[side[0]];
    int d2 = d0[side[1]];
    int a = d1 + 1;
    int b = d2 + 1;

    // Проверка на произведение
    if (1LL * a * b != n) {
        swap(a, b);
        if (1LL * a * b != n) {
            cout << "NO\n";
            return 0;
        }
    }

    // Проверяем распределение степеней вершин
    int cnt2 = 0, cnt3 = 0, cnt4 = 0;
    for (int x : deg) {
        if (x == 2)      cnt2++;
        else if (x == 3) cnt3++;
        else if (x == 4) cnt4++;
        else {
            cout << "NO\n";
            return 0;
        }
    }

    int exp2 = 4;
    int exp3 = 2 * max(0, a - 2) + 2 * max(0, b - 2);
    int exp4 = max(0, a - 2) * max(0, b - 2);

    if (cnt2 == exp2 && cnt3 == exp3 && cnt4 == exp4) {
        cout << "YES\n" << a << " " << b << "\n";
    } else {
        cout << "NO\n";
    }

    return 0;
}
