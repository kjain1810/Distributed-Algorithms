#include <iostream>
#include <vector>

int main() {
  int n;
  std::cin >> n;
  std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++)
      std::cin >> adj[a][b];
  int ans = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++)
        if (adj[a][b] && adj[a][c] && adj[b][c])
          ans++;
  std::cout << ans << std::endl;
  return 0;
}
