#include <ctime>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <time.h>

int main() {
  float x;
  std::cin >> x;
  int accuracy;
  std::cin >> accuracy;
  clock_t ts = clock();
  float toadd = x;
  float ans = 0;
  for (int a = 1; a <= accuracy; a++) {
    ans += toadd;
    toadd *= x;
    toadd *= x;
    toadd /= (2 * a);
    toadd /= (2 * a + 1);
    toadd *= -1;
  }
  std::cout << "Answer: " << std::fixed << std::setprecision(20) << ans << "\n";
  std::cout << "Time: " << (double)(clock() - ts) / CLOCKS_PER_SEC;
  return 0;
}
