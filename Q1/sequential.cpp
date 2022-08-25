#include <iostream>
#include <vector>

int main() {
  int n;
  std::cin >> n;
  std::vector<std::vector<float>> A(n, std::vector<float>(n));
  /* std::vector<std::vector<float>> I(n, std::vector<float>(n)); */
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++) {
      std::cin >> A[a][b];
      /* if (a == b) */
      /*   I[a][b] = 1; */
      /* else */
      /*   I[a][b] = 0; */
    }
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++)
      A[a].push_back((int)(a == b));
  int h = 0;
  int k = 0;
  while (h < n && k < n) {
    // find pivot
    int i_max = -1;
    float val_max = -1e9;
    for (int i = h; i < n; i++)
      if (std::abs(A[i][k]) >= val_max) {
        val_max = std::abs(A[i][k]);
        i_max = i;
      }
    if (A[i_max][k] == 0) {
      // no pivot for this column
      k++;
      continue;
    }
    float pivot = val_max;
    // swap hth and i_maxth row
    for (int i = 0; i < 2 * n; i++) {
      std::swap(A[i_max][i], A[h][i]);
    }
    // set pivot
    // for all rows below the pivot:
    for (int i = h + 1; i < n; i++) {
      float f = A[i][k] / pivot;
      A[i][k] = 0;
      /* I[i][k] = 0; */
      for (int j = k + 1; j < 2 * n; j++) {
        A[i][j] = A[i][j] - A[h][j] * f;
        /* I[i][j] = I[i][j] - I[h][j] * f; */
      }
    }
    h++;
    k++;

    std::cout << h << " " << k << "\n";
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 2 * n; j++) {
        if (j == n)
          std::cout << "| ";
        std::cout << A[i][j] << " ";
      }
      std::cout << "\n";
    }
  }
  for (int i = 0; i < n; i++) {
    float here = A[i][i];
    for (int j = 0; j < 2 * n; j++)
      A[i][j] /= here;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2 * n; j++) {
      if (j == n)
        std::cout << "| ";
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }

  /* // BACK SUBSTITUTION PHASE */
  for (int a = n - 1; a >= 1; a--) {
    for (int b = a - 1; b >= 0; b--) {
      float factor = A[b][a];
      for (int c = 0; c < 2 * n; c++) {
        A[b][c] = A[b][c] - factor * A[a][c];
        /* I[b][c] = I[b][c] - factor * I[a][c]; */
      }
    }
  }

  std::cout << "Successful\n";
  for (int a = 0; a < n; a++) {
    for (int b = n; b < 2 * n; b++)
      std::cout << A[a][b] << " ";
    std::cout << "\n";
  }
  return 0;
}
