#!/bin/bash 
#                                                    
                                  
#SBATCH --job-name=janv                                      
#SBATCH --time=01:59:00                                            
#SBATCH --ntasks=192                                  
#SBATCH --ntasks-per-node=48                                                                       
#SBATCH --output=ZOutGrace_192.txt
#SBATCH --account=132780049402
#SBATCH --mem=360G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janv@tamu.edu

## module load GCC/11.3.0
## module load iimpi/2022.00
## module load CMake/3.22.1

## YOUR COMMANDS BELOW  
CHI_TECH=/scratch/group/tees-class/gcc-11.3.0/ChiTech/naktakala-chi-tech/bin/ChiTech

mpiexec -np 192 $CHI_TECH X_Sim1.lua
