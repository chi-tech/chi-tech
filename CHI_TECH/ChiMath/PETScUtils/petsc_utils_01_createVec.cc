#include "petsc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a general vector.
 *
This is a macro for:
\code
Vec x;
VecCreate(PETSC_COMM_WORLD,&x);
VecSetType(x,VECMPI);
VecSetSizes(x, local_size, global_size);
VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

return x;
\endcode*/
Vec chi_math::PETScUtils::
CreateVector(int local_size, int global_size)
{
  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetType(x,VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

  return x;
}

//###################################################################
/**Creates a general vector.
 *
This is a macro for:
\code
VecCreate(PETSC_COMM_WORLD,&x);
VecSetType(x,VECMPI);
VecSetSizes(x, local_size, global_size);
VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
\endcode*/
void chi_math::PETScUtils::
CreateVector(Vec& x, int local_size, int global_size)
{
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetType(x,VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
}

//###################################################################
/**Creates a general vector with ghost value support.
 *
This is a macro for:
\code
Vec x;
VecCreateGhost(PETSC_COMM_WORLD,
               local_size,
               global_size,
               nghosts,
               ghost_indices.data(),
               &x);

VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

return x;
\endcode*/
Vec chi_math::PETScUtils::
CreateVectorWithGhosts(int local_size, int global_size,
                       int nghosts,
                       std::vector<int>& ghost_indices)
{
  Vec x;
  VecCreateGhost(PETSC_COMM_WORLD,
                 local_size,
                 global_size,
                 nghosts,
                 (ghost_indices.empty())? NULL : ghost_indices.data(),
                 &x);

  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

  return x;
}