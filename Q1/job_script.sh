#!/bin/bash
#SBATCH -N 2
#SBATCH --mem-per-cpu=1024

# module load hpcx-2.7.0/hpcx-ompi

mpic++ main.cpp

# rm out
# touch out

mpirun --use-hwthread-cpus a.out > out

# echo "Running 100" >> out
# mpirun --use-hwthread-cpus a.out < input100.txt >> out

# echo "Running 1000" >> out
# mpirun --use-hwthread-cpus a.out < input1000.txt >> out

# echo "Running 10000" >> out
# mpirun --use-hwthread-cpus a.out < input10000.txt >> out
