#include <bits/stdc++.h>
#include <random>

using namespace std;

int main() {
  int n;
  cin >> n;
  cout << n << endl;
  mt19937 gen(time(NULL));
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++)
      cout << gen() << " ";
    cout << endl;
  }
  return 0;
}
