#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <vector>

struct pivotMesg {
  long double val;
  int idx;
};

void swapping(int i_max_overall, int h, int myid, int num_procs, int n,
              long double **arr) {
  MPI_Status status;
  if (i_max_overall % num_procs == myid && h % num_procs == myid) {
    // both to swap are here
    int idx1 = i_max_overall / num_procs;
    int idx2 = h / num_procs;
    for (int i = 0; i < n; i++)
      std::swap(arr[idx1][i], arr[idx2][i]);
  } else if (i_max_overall % num_procs == myid) {
    // send over the i_max_overall and recieve hth
    int idx1 = i_max_overall / num_procs;
    int dest = h % num_procs;
    MPI_Send(arr[idx1], n, MPI_LONG_DOUBLE, dest, 1, MPI_COMM_WORLD);
    MPI_Recv(arr[idx1], n, MPI_LONG_DOUBLE, dest, 1, MPI_COMM_WORLD, &status);
  } else if (h % num_procs == myid) {
    // recieve i_maxth and send over hth
    long double new_row[n];
    int src = i_max_overall % num_procs;
    int idx1 = h / num_procs;
    MPI_Recv(new_row, n, MPI_LONG_DOUBLE, src, 1, MPI_COMM_WORLD, &status);
    MPI_Send(arr[idx1], n, MPI_LONG_DOUBLE, src, 1, MPI_COMM_WORLD);
    for (int i = 0; i < n; i++)
      arr[idx1][i] = new_row[i];
    /* MPI_Bcast(arr[idx1], n, MPI_LONG_DOUBLE, myid, MPI_COMM_WORLD); */
  }
}

void get_hth_row_and_reduce(int n, int num_procs, int num_rows, int myid, int h,
                            int k, long double pivot, long double **arr) {
  // get the hth row now
  long double vec[n];
  int root = h % num_procs;
  if (root == myid)
    for (int i = 0; i < n; i++)
      vec[i] = arr[h / num_procs][i];
  MPI_Bcast(vec, n, MPI_LONG_DOUBLE, root, MPI_COMM_WORLD);
  // reduce the row with the pivot and the hth row
  std::cout << "on iteration " << k << " proc " << myid << " got: ";
  for (int i = 0; i < n; i++)
    std::cout << vec[i] << " ";
  std::cout << std::endl;
  for (int i = myid; i < n; i += num_procs) {
    if (i <= h)
      continue;
    int idx = i / num_procs;
    long double f = arr[idx][k] / pivot;
    arr[idx][k] = 0;
    for (int j = k + 1; j < n; j++)
      arr[idx][j] -= vec[j] * f;
  }
}

int main(int argc, char *argv[]) {
  int n;
  int num_rows;
  int num_procs;
  int myid;
  long double **arr;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int num_types = 2;
  int blockSize[] = {1, 1};
  MPI_Datatype types[] = {MPI_LONG_DOUBLE, MPI_INT};
  MPI_Aint offsets[] = {offsetof(pivotMesg, val), offsetof(pivotMesg, idx)};
  MPI_Datatype pivotMesg_mpi;

  MPI_Type_create_struct(num_types, blockSize, offsets, types, &pivotMesg_mpi);
  MPI_Type_commit(&pivotMesg_mpi);

  if (myid == 0) {
    std::cin >> n;
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  num_rows = n / num_procs;
  if (myid < n % num_procs)
    num_rows++;
  arr = (long double **)malloc(num_rows * sizeof(long double *));
  for (int i = 0; i < n; i++)
    arr[i] = (long double *)malloc(n * sizeof(long double));
  /* std::cout << myid << ": " << n << " " << num_rows << std::endl; */
  if (myid == 0) {
    // Take array as input
    long double *vec = (long double *)malloc(n * sizeof(long double));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        std::cin >> vec[j];
      int target = i % num_procs;
      // If it belongs to some other process, send it over
      if (target != 0)
        MPI_Send(vec, n, MPI_LONG_DOUBLE, target, i, MPI_COMM_WORLD);
      // Else store it here
      else {
        int row_idx_here = i / num_procs;
        for (int j = 0; j < n; j++)
          arr[row_idx_here][j] = vec[j];
      }
    }
  } else {
    // Recieve the rows of this process from the process 0
    long double *vec = (long double *)malloc(n * sizeof(long double));
    for (int i = 0; i < num_rows; i++) {
      MPI_Status status;
      MPI_Recv(vec, n, MPI_LONG_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
               &status);
      int row_idx = status.MPI_TAG;
      int row_idx_here = row_idx / num_procs;
      for (int j = 0; j < n; j++)
        arr[row_idx_here][j] = vec[j];
    }
  }
  /* std::cout << myid << " reached barrier!" << std::endl; */
  /* for (int i = 0; i < num_rows; i++) { */
  /*   std::cout << num_procs * i + myid << ": "; */
  /*   for (int j = 0; j < n; j++) */
  /*     std::cout << arr[i][j] << " "; */
  /*   std::cout << std::endl; */
  /* } */
  // Wait for everyone to get their arrays
  MPI_Barrier(MPI_COMM_WORLD);
  // do gauss jordan elimination
  //
  // MPI Communication Tags:
  // 0 - for communicating the pivot max
  // 1 - for swapping the rows

  int h = 0;
  for (int k = 0; k < n; k++) {
    std::cout << myid << " : " << k << std::endl;
    int i_max = -1;
    long double val_max = -1;
    for (int i = myid; i < n; i += num_procs) {
      int idx = i / num_procs;
      if (i >= h && std::abs(arr[idx][k]) > val_max) {
        val_max = std::abs(arr[idx][k]);
        i_max = i;
      }
    }
    pivotMesg mesg;
    long double pivot;
    if (myid == 0) {
      MPI_Status status;
      for (int i = 1; i < num_procs; i++) {
        MPI_Recv(&mesg, 1, pivotMesg_mpi, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
                 &status);
        if (mesg.val > val_max) {
          val_max = mesg.val;
          i_max = mesg.idx;
        }
      }
      pivot = val_max;
    } else {
      mesg.val = val_max;
      mesg.idx = i_max;
      MPI_Send(&mesg, 1, pivotMesg_mpi, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&pivot, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    if (pivot == 0)
      continue;
    MPI_Bcast(&i_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << "proc " << myid << " with pivot " << pivot << " of row "
              << i_max << std::endl;
    swapping(i_max, h, myid, num_procs, n, arr);
    get_hth_row_and_reduce(n, num_procs, num_rows, myid, h, k, pivot, arr);
    h++;
    std::cout << "after iteration " << k << ": " << std::endl;
    for (int i = myid; i < n; i += num_procs) {
      int idx = i / num_procs;
      std::cout << i << ": ";
      for (int j = 0; j < n; j++)
        std::cout << arr[idx][j] << " ";
      std::cout << std::endl;
    }
  }
  for (int i = myid; i < n; i += num_procs) {
    int idx = i / num_procs;
    long double here = arr[idx][i];
    for (int j = 0; j < n; j++)
      arr[idx][j] /= here;
  }
  // wait for everyone to finish with gauss jordan
  MPI_Barrier(MPI_COMM_WORLD);
  // start back substitution
  for (int a = n - 1; a >= 0; a--) {
  }
  MPI_Finalize();
  return 0;
}
