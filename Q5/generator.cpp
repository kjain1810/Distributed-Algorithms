#include <bits/stdc++.h>
#include <random>

using namespace std;

int main() {
  mt19937 gen(time(NULL));
  int n;
  cin >> n;
  vector<vector<int>> adjMat(n, vector<int>(n, 0));
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      adjMat[i][j] = adjMat[j][i] = gen() % 2;
    }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      cout << adjMat[i][j] << " ";
    cout << "\n";
  }
  return 0;
}
