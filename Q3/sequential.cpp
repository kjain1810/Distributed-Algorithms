#include <iostream>
#include <iomanip>
#include <math.h>

int main()
{
  float x;
  std::cin >> x;
  int accuracy;
  std::cin >> accuracy;
  float toadd = x;
  float ans = 0;
  for(int a = 1; a <= accuracy; a++)
  {
    ans += toadd;
    toadd *= x;
    toadd *= x;
    toadd /= (2 * a);
    toadd /= (2 * a + 1);
    toadd *= -1;
  }
  std::cout << std::fixed << std::setprecision(20) << ans << "\n";
  return 0; 
}
