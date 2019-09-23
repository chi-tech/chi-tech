#include <LinearBoltzmanSolver/lbs_linear_boltzman_solver.h>
#include <petscksp.h>



int NPTMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax);
int NPTMatrixAction_Ax_Cycles(Mat matrix, Vec krylov_vector, Vec Ax);