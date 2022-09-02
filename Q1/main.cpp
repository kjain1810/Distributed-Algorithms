#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <vector>

struct pivotMesg {
  float val;
  int idx;
};

void swapping(int i_max_overall, int h, int myid, int num_procs, int n,
              float **arr) {
  MPI_Status status;
  if (i_max_overall % num_procs == myid && h % num_procs == myid) {
    // both to swap are here
    int idx1 = i_max_overall / num_procs;
    int idx2 = h / num_procs;
    for (int i = 0; i < 2 * n; i++)
      std::swap(arr[idx1][i], arr[idx2][i]);
  } else if (i_max_overall % num_procs == myid) {
    // send over the i_max_overall and recieve hth
    int idx1 = i_max_overall / num_procs;
    int dest = h % num_procs;
    MPI_Send(arr[idx1], 2 * n, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Recv(arr[idx1], 2 * n, MPI_FLOAT, dest, 1, MPI_COMM_WORLD, &status);
  } else if (h % num_procs == myid) {
    // recieve i_maxth and send over hth
    float new_row[2 * n];
    int src = i_max_overall % num_procs;
    int idx1 = h / num_procs;
    MPI_Recv(new_row, 2 * n, MPI_FLOAT, src, 1, MPI_COMM_WORLD, &status);
    MPI_Send(arr[idx1], 2 * n, MPI_FLOAT, src, 1, MPI_COMM_WORLD);
    for (int i = 0; i < 2 * n; i++)
      arr[idx1][i] = new_row[i];
    /* MPI_Bcast(arr[idx1], n, MPI_FLOAT, myid, MPI_COMM_WORLD); */
  }
}

void get_hth_row_and_reduce(int n, int num_procs, int num_rows, int myid, int h,
                            int k, float pivot, float **arr) {
  // get the hth row now
  float vec[2 * n];
  int root = h % num_procs;
  if (root == myid)
    for (int i = 0; i < 2 * n; i++)
      vec[i] = arr[h / num_procs][i];
  MPI_Bcast(vec, 2 * n, MPI_FLOAT, root, MPI_COMM_WORLD);
  // reduce the row with the pivot and the hth row
  /* std::cout << "on iteration " << k << " proc " << myid << " got: "; */
  /* for (int i = 0; i < 2 * n; i++) */
  /*   std::cout << vec[i] << " "; */
  /* std::cout << "\n"; */
  for (int i = myid; i < n; i += num_procs) {
    if (i <= h)
      continue;
    int idx = i / num_procs;
    float f = arr[idx][k] / pivot;
    arr[idx][k] = 0;
    for (int j = k + 1; j < 2 * n; j++)
      arr[idx][j] -= vec[j] * f;
  }
}

int main(int argc, char *argv[]) {
  int n;
  int num_rows;
  int num_procs;
  int myid;
  float **arr;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int num_types = 2;
  int blockSize[] = {1, 1};
  MPI_Datatype types[] = {MPI_FLOAT, MPI_INT};
  MPI_Aint offsets[] = {offsetof(pivotMesg, val), offsetof(pivotMesg, idx)};
  MPI_Datatype pivotMesg_mpi;

  MPI_Type_create_struct(num_types, blockSize, offsets, types, &pivotMesg_mpi);
  MPI_Type_commit(&pivotMesg_mpi);

  std::freopen("./input1000.txt", "r", stdin);
  std::cin >> n;

  num_rows = n / num_procs;
  /* if (num_rows == 0) { */
  /*   MPI_Finalize(); */
  /*   return 0; */
  /* } */
  if (myid < n % num_procs)
    num_rows++;
  arr = (float **)malloc((num_rows + 1) * sizeof(float *));
  for (int i = 0; i < num_rows; i++)
    arr[i] = (float *)malloc((2 * n + 1) * sizeof(float));

  /* std::cout << myid << ": " << n << " " << num_rows << "\n"; */
  MPI_Barrier(MPI_COMM_WORLD);
  int idxHere = 0;
  for (int i = 0; i < n; i++) {
    /* std::cout << myid << " " << i << "\n"; */
    if (i % num_procs == myid) {
      for (int j = 0; j < n; j++) {
        std::cin >> arr[idxHere][j];
        if (j == i)
          arr[idxHere][n + j] = 1.0;
        else
          arr[idxHere][n + j] = 0.0;
      }
      idxHere++;
    } else {
      float tmp;
      for (int j = 0; j < n; j++)
        std::cin >> tmp;
    }
    /* if (idxHere > num_rows) */
    /*   std::cout << "wut ze fuk?\n"; */
  }
  /* std::cout << myid << " reached input barrier!" */
  /*           << "\n"; */
  /* for (int i = 0; i < num_rows; i++) { */
  /*   std::cout << num_procs * i + myid << ": "; */
  /*   for (int j = 0; j < n; j++) */
  /*     std::cout << arr[i][j] << " "; */
  /*   std::cout << "\n"; */
  /* } */
  // Wait for everyone to get their arrays
  MPI_Barrier(MPI_COMM_WORLD);
  // do gauss jordan elimination
  //
  // MPI Communication Tags:
  // 0 - for communicating the pivot max
  // 1 - for swapping the rows
  double start_time = MPI_Wtime();

  int h = 0;
  for (int k = 0; k < n; k++) {
    /* std::cout << myid << " : " << k << "\n"; */

    // finding the local pivot
    int i_max = -1;
    float val_max = -1;
    for (int i = myid; i < n; i += num_procs) {
      int idx = i / num_procs;
      if (i >= h && std::abs(arr[idx][k]) > val_max) {
        val_max = std::abs(arr[idx][k]);
        i_max = i;
      }
    }

    // collect the pivot at process 0
    pivotMesg mesg;
    float pivot;
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
    // broadcast the pivot
    MPI_Bcast(&pivot, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // if no pivot, move forward
    if (pivot == 0)
      continue;
    // broadcast the pivot row
    MPI_Bcast(&i_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* std::cout << "proc " << myid << " with pivot " << pivot << " of row " */
    /*           << i_max << "\n"; */
    // swap current and pivot rows
    swapping(i_max, h, myid, num_procs, n, arr);
    // do operation using pivot row
    get_hth_row_and_reduce(n, num_procs, num_rows, myid, h, k, pivot, arr);
    h++;
    /* std::cout << "after iteration " << k << ": " << "\n"; */
    /* for (int i = myid; i < n; i += num_procs) { */
    /*   int idx = i / num_procs; */
    /*   std::cout << i << ": "; */
    /*   for (int j = 0; j < 2 * n; j++) */
    /*     std::cout << arr[idx][j] << " "; */
    /*   std::cout << "\n"; */
    /* } */
  }
  for (int i = myid; i < n; i += num_procs) {
    int idx = i / num_procs;
    float here = arr[idx][i];
    for (int j = 0; j < 2 * n; j++)
      arr[idx][j] /= here;
  }
  // wait for everyone to finish with gauss jordan
  MPI_Barrier(MPI_COMM_WORLD);

  // start back substitution
  for (int a = n - 1; a >= 1; a--) {
    int src_proc = a % num_procs;
    float vec[2 * n];
    if (src_proc == myid) {
      int idx = a / num_procs;
      for (int i = 0; i < 2 * n; i++)
        vec[i] = arr[idx][i];
    }
    MPI_Bcast(vec, 2 * n, MPI_FLOAT, src_proc, MPI_COMM_WORLD);
    for (int i = myid; i < a; i += num_procs) {
      int idx = i / num_procs;
      float factor = arr[idx][a];
      for (int j = 0; j < 2 * n; j++)
        arr[idx][j] = arr[idx][j] - factor * vec[j];
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  double cur_time = MPI_Wtime();
  double end_time;
  MPI_Reduce(&cur_time, &end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myid == 0) {
    std::cout << "Runtime: " << (end_time - start_time) << "\n";
  }

  std::string outMatrixName = "./outMatrix.txt";
  freopen(outMatrixName.c_str(), "w+", stdout);
  for (int i = 0; i < n; i++) {
    if (i % num_procs == myid) {
      for (int j = 0; j < 2 * n; j++)
        std::cout << arr[i / num_procs][j] << " ";
      std::cout << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  /* for (int a = myid; a < n; a += num_procs) { */
  /*   std::cout << "Final " << a << ": "; */
  /*   int idx = a / num_procs; */
  /*   for (int b = 0; b < 2 * n; b++) */
  /*     std::cout << arr[idx][b] << " "; */
  /*   std::cout << "\n"; */
  /* } */
  if (myid == 0)
    std::cout << n << " " << num_procs << "\n";
  MPI_Finalize();
  return 0;
}
