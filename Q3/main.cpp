#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char *argv[]) {
  double start_time = MPI_Wtime();
  long double x;
  int accuracy;
  int pow_begin;
  int pow_end;
  int iter_start;
  int iter_end;
  long double toadd;
  long double res_here;

  int numprocs, myid;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    std::cin >> x;
    std::cin >> accuracy;
    /* x = 8; */
    /* accuracy = 100; */
  }
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&accuracy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  toadd = -1;
  res_here = 0;
  // power of first term
  iter_start = (accuracy / numprocs) * myid;
  pow_begin = iter_start * 2 + 1;
  // power of last term
  iter_end = (accuracy / numprocs) * (myid + 1) - 1;
  if (myid == numprocs - 1)
    iter_end = accuracy - 1;
  pow_end = iter_end * 2 + 1;
  if ((pow_begin + 1) % 4 == 0)
    toadd = 1;
  // take (x ^ (pow_begin - 2)) / (pow_begin - 2)! common from all terms
  for (int a = pow_begin; a <= pow_end; a += 2) {
    if (a - 1 > 0)
      toadd *= (x / (a - 1));
    toadd *= (x / a);
    toadd *= -1;
    res_here += toadd;
    /* std::cout << myid << " " << a << " " << toadd << "\n"; */
  }
  if (myid == 0) {
    // collect the result and last term from each process
    std::vector<long double> res_proc(numprocs);
    std::vector<long double> last_term(numprocs);
    res_proc[0] = res_here;
    last_term[0] = toadd;
    /* std::cout << "Have: " << res_proc[0] << " 0 0\n"; */
    /* std::cout << "Have: " << last_term[0] << "0 1\n"; */
    for (int a = 2; a < 2 * numprocs; a++) {
      MPI_Status stat;
      float data;
      MPI_Recv(&data, 1, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &stat);
      int procid = stat.MPI_SOURCE;
      int tag = stat.MPI_TAG;
      /* std::cout << "Recieved: " << data << " " << procid << " " << tag <<
       * "\n"; */
      if (tag == 0)
        res_proc[procid] = data;
      else
        last_term[procid] = data;
    }
    // iteratate over the result of each process and multiply
    float ans = res_proc[0];
    for (int a = 1; a < numprocs; a++) {
      ans += res_proc[a] * abs(last_term[a - 1]);
      /* last_term[a] *= abs(last_term[a - 1]); */
    }
    std::cout << "Answer for x = " << x << " with " << accuracy
              << " iteration: " << ans << "\n";
  } else {
    // send the result from this process to the main process
    MPI_Send(&res_here, 1, MPI_LONG_DOUBLE, 0, 0,
             MPI_COMM_WORLD); // res_here sent with tag 0
    MPI_Send(&toadd, 1, MPI_LONG_DOUBLE, 0, 1,
             MPI_COMM_WORLD); // toadd sent with tag 1
  }
  MPI_Finalize();
  if (myid == 0) {
    double end_time = MPI_Wtime();
    std::cout << "Runtime: " << (end_time - start_time);
  }
  return 0;
}
