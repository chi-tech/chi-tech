#include <LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h>
#include <petscksp.h>



int LBSMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax);
