#include <iostream>
#include <vector>

int main() {
  int n;
  std::cin >> n;
  std::vector<std::vector<int>> A(n, std::vector<int>(n, 0));
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++)
      std::cin >> A[a][b];

  // B = A * A
  std::vector<std::vector<int>> B(n, std::vector<int>(n, 0));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        B[i][k] += A[i][j] * A[j][k];

  // C = B * A
  std::vector<std::vector<int>> C(n, std::vector<int>(n, 0));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        C[i][k] += B[i][j] * A[j][k];
  int ans = 0;
  for (int i = 0; i < n; i++)
    ans += C[i][i];
  std::cout << ans / 2 << std::endl;

  ans = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++)
        if (A[a][b] && A[a][c] && A[b][c])
          ans++;
  std::cout << ans << std::endl;
  return 0;
}
