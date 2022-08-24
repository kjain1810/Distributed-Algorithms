#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem-per-cpu=1024

# module load hpcx-2.7.0/hpcx-ompi

mpic++ main.cpp

mpirun --use-hwthread-cpus a.out < input.txt > out
