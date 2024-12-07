#include <chrono>
#include <iostream>
#include <string>

#include "solver.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <number_of_partitions>" << endl;
    return 1;
  }
  const auto n = stoull(argv[1]);
  const auto m = n;
  const auto eps = 1e-6;
  const auto max_iter = 999999999u;
  const auto K = 16u;

  Solver S(n, m, eps, max_iter, K);

  const auto begin = chrono::high_resolution_clock::now();
  S.Solve();
  const auto end = chrono::high_resolution_clock::now();
  const chrono::duration<double> elapsed = end - begin;

  cout << "ВРЕМЯ ИСПОЛНЕНИЯ: " << elapsed.count() << " сек." << endl << endl;

  cout << "СПРАВКА О РЕШЕНИИ:" << endl;
  cout << "Для решения тестовой задачи использована сетка-основа с числом разбиений по x n=«" << n
       << "» и числом разбиений по y m=«" << m << "».\n";
  cout << "Метод Чебышева (" << K << "), критерии остановки по точности εмет=«" << eps << "» и по числу итераций Nmax=«"
       << max_iter << "».\n\n";

  cout << "На решение схемы (СЛАУ) затрачено итераций N=«" << S.GetIterationCount()
       << "» и достигнута точность итерационного метода ε(N)=«" << S.GetAccuracy() << "».\n";
  cout << "Схема (СЛАУ) решена с невязкой ||R(N)|| = «" << S.GetResultDiscrepancy()
       << "» для невязки СЛАУ использована норма «max».\n\n";

  constexpr auto error_threshold = 0.5e-6;
  cout << "Тестовая задача должна быть решена с погрешностью не более ε = " << error_threshold
       << "; задача решена с погрешностью ε1=«" << S.GetMaxDiff() << "».\n\n";

  cout << "Максимальное отклонение точного и численного решений наблюдается в узле x=«" << S.GetMaxDiffX() << "», y=«"
       << S.GetMaxDiffY() << "».\n\n";

  cout << "В качестве начального приближения использовано «нулевое».\n";
  cout << "Невязка СЛАУ на начальном приближении ||R(0)|| = «" << S.GetInitialDiscrepancy() << "» (max).\n";

  return 0;
}
