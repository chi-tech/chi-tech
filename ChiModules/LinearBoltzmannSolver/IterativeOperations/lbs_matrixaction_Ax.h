#ifndef LBS_MATRIXACTION_AX_H
#define LBS_MATRIXACTION_AX_H

#include <LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h>
#include <petscksp.h>

namespace LinearBoltzmann
{
int LBSMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax);
}


#endif