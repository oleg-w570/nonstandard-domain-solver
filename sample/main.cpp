// Copyright Zorin Oleg
#include <chrono>
#include <iostream>

#include "solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
//   const int n = atoi(argv[1]);
  const int n = 256;
  const int m = n;
  const double eps = 1e-6;
  const int max_iter = 999999;
  const int K = 16;

  Solver S(n, m, eps, max_iter, K);

  const auto begin = chrono::high_resolution_clock::now();
  S.Solve();
  const auto end = chrono::high_resolution_clock::now();
  const chrono::duration<double> elapsed = end - begin;

  cout << "Время исполнения: " << elapsed.count() << " сек." << endl << endl;

  cout << "Для решения тестовой задачи использована сетка-основа с числом "
          "разбиений по x n=«"
       << n
       << "» "
          "и числом разбиений по y m=«"
       << m << "»." << endl;

  cout << "Метод Чебышева (" << K << "), критерии остановки по точности εмет=«" << eps << "» "
       << "и по числу итераций Nmax=«" << max_iter << "»." << endl
       << endl;

  cout << "На решение схемы (СЛАУ) затрачено итераций N=«" << S.n_iter << "» "
       << "и достигнута точность итерационного метода ε(N)=«" << S.accuracy << "»." << endl;

  cout << "Схема (СЛАУ) решена с невязкой ||R(N)|| = «" << S.R_res << "» "
       << "для невязки СЛАУ использована норма «max»." << endl
       << endl;

  cout << "Тестовая задача должна быть решена с погрешностью не более ε = "
          "0.5⋅10^–6; "
       << "задача решена с погрешностью ε1=«" << S.max_diff << "»." << endl;

  cout << "Максимальное отклонение точного и численного решений наблюдается в узле x=«"
       << S.max_diff_x << "», y=«" << S.max_diff_y << "»." << endl
       << endl;

  cout << "В качестве начального приближения использовано «нулевое»." << endl;

  cout << "Невязка СЛАУ на начальном приближении ||R(0)|| = «" << S.R_null << "» (max)." << endl;

  return 0;
}
