#include <ctime>
#include <ios>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  double start_time;
  int numprocs, myid;
  int dims[2] = {0, 0};
  int n_here;
  int n;
  int left, right, up, down;
  int periods[2];
  MPI_Comm cannon_comm;
  long long int *A;
  long long int *A_copy;
  long long int *A_square;
  long long int *A_cube;
  long long int *A_four;
  long long int *buf;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  MPI_Dims_create(numprocs, 2, dims);
  if (dims[0] != dims[1]) {
    if (myid == 0)
      std::cout << "Need a square number of processors!\n";
    MPI_Finalize();
    return 0;
  }
  freopen("matrix10000.txt", "r", stdin);
  std::cin >> n;
  n_here = n / dims[0];

  int i_idx = myid / dims[0];
  int j_idx = myid % dims[0];

  int i_start = i_idx * n_here;
  int i_stop = (i_idx + 1) * n_here - 1;
  int j_start = j_idx * n_here;
  int j_stop = (j_idx + 1) * n_here - 1;

  A = (long long int *)malloc(n_here * n_here * sizeof(long long int));
  A_copy = (long long int *)malloc(n_here * n_here * sizeof(long long int));
  A_square = (long long int *)malloc(n_here * n_here * sizeof(long long int));
  A_cube = (long long int *)malloc(n_here * n_here * sizeof(long long int));
  A_four = (long long int *)malloc(n_here * n_here * sizeof(long long int));
  buf = (long long int *)malloc(n_here * n_here * sizeof(long long int));

  int idx = 0;
  for (int i = 0; i < n_here; i++) {
    for (int j = 0; j < n_here; j++) {
      long long x;
      std::cin >> x;
      if (i >= i_start && i <= i_stop && j >= j_start && j <= j_stop) {
        A_square[idx] = 0;
        A_cube[idx] = 0;
        A_four[idx] = 0;
        A_copy[idx] = x;
        A[idx++] = x;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  start_time = clock();

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cannon_comm);
  MPI_Cart_shift(cannon_comm, 0, 1, &left, &right);
  MPI_Cart_shift(cannon_comm, 1, 1, &up, &down);

  long long int *B =
      (long long int *)malloc(n_here * n_here * sizeof(long long int));
  for (int i = 0; i < n_here * n_here; i++)
    B[i] = A[i];

  // Computing A_square
  for (int shift = 0; shift < dims[0]; shift++) {
    for (int i = 0; i < n_here; i++)
      for (int j = 0; j < n_here; j++)
        for (int k = 0; k < n_here; k++)
          A_square[i * n_here + j] += A[i * n_here + k] * B[k * n_here + j];
    if (shift == dims[0] - 1)
      break;
    MPI_Status status;
    MPI_Sendrecv(A, n_here * n_here, MPI_LONG_LONG_INT, left, 1, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, right, 1, cannon_comm,
                 &status);
    long long int *tmp = buf;
    buf = A;
    A = tmp;
    MPI_Sendrecv(B, n_here * n_here, MPI_LONG_LONG_INT, up, 2, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, down, 2, cannon_comm,
                 &status);
    tmp = buf;
    buf = B;
    B = tmp;
  }

  long long M =
      0; // Sum of matrix A_square, needed for computing cycles of length 4

  for (int i = 0; i < n_here * n_here; i++)
    M += A_square[i];

  // Now, A is the same, so we can copy it from A_copy
  for (int i = 0; i < n_here * n_here; i++)
    A[i] = A_copy[i];
  // B is A_square
  for (int i = 0; i < n_here * n_here; i++)
    B[i] = A_square[i];

  // compute A_cube now
  for (int shift = 0; shift < dims[0]; shift++) {
    for (int i = 0; i < n_here; i++)
      for (int j = 0; j < n_here; j++)
        for (int k = 0; k < n_here; k++)
          A_cube[i * n_here + j] += A[i * n_here + k] * B[k * n_here + j];
    if (shift == dims[0] - 1)
      break;
    MPI_Status status;
    MPI_Sendrecv(A, n_here * n_here, MPI_LONG_LONG_INT, left, 1, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, right, 1, cannon_comm,
                 &status);
    std::swap(A, buf);
    MPI_Sendrecv(B, n_here * n_here, MPI_LONG_LONG_INT, up, 2, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, down, 2, cannon_comm,
                 &status);
    std::swap(B, buf);
  }

  long long diagonal_element_sum = 0; // sum of diagonal elements of A^3
  for (int i = 0; i < n_here; i++)
    for (int j = 0; j <= n_here; j++)
      if (i + i_start == j + j_start)
        diagonal_element_sum += A_cube[i * n_here + j];

  // Now, A is the same, so we can copy it from A_copy
  for (int i = 0; i < n_here * n_here; i++)
    A[i] = A_copy[i];
  // B is A_cube
  for (int i = 0; i < n_here * n_here; i++)
    B[i] = A_cube[i];

  // compute A_four now
  for (int shift = 0; shift < dims[0]; shift++) {
    for (int i = 0; i < n_here; i++)
      for (int j = 0; j < n_here; j++)
        for (int k = 0; k < n_here; k++)
          A_four[i * n_here + j] += A[i * n_here + k] * B[k * n_here + j];
    if (shift == dims[0] - 1)
      break;
    MPI_Status status;
    MPI_Sendrecv(A, n_here * n_here, MPI_LONG_LONG_INT, left, 1, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, right, 1, cannon_comm,
                 &status);
    std::swap(A, buf);
    MPI_Sendrecv(B, n_here * n_here, MPI_LONG_LONG_INT, up, 2, buf,
                 n_here * n_here, MPI_LONG_LONG_INT, down, 2, cannon_comm,
                 &status);
    std::swap(B, buf);
  }

  long long int S = 0; // sum of diagonal elements of A^4, required for
                       // calculating cycles of length 4
  for (int i = 0; i < n_here; i++)
    for (int j = 0; j < n_here; j++)
      if (i + i_start == j + j_start)
        S += A_four[i * n_here + j];

  MPI_Barrier(cannon_comm);

  // collect S, M, diagonal_element_sum at process 0
  long long int sumS, sumM, sum_diagonal_elemet_sum;
  MPI_Reduce(&S, &sumS, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&M, &sumM, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&diagonal_element_sum, &sum_diagonal_elemet_sum, 1,
             MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid == 0) {
    std::cout << "C3: " << diagonal_element_sum << std::endl;
    std::cout << "Length 4 cycles: " << (S - M) / 4 << std::endl;
    std::cout << "Runtime: " << (clock() - start_time) / CLOCKS_PER_SEC
              << std::endl;
  }
  MPI_Finalize();
  return 0;
}
