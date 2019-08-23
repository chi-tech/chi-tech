#!/bin/bash

#SBATCH --job-name=Test.i
#SBATCH --output=Test.i.log
#SBATCH --error=Test.i.error
#SBATCH --ntasks=20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=janv@tamu.edu

module load mpi/openmpi-x86_64

mpirun -np 4 ./bin/ChiTech CHI_TEST/Transport3D_3BlockPoly.lua 
