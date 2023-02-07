#ifndef CHITECH_LBSMIP_SHELL_OPERATIONS_H
#define CHITECH_LBSMIP_SHELL_OPERATIONS_H

#include <petscksp.h>

namespace lbs
{
int MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);

}//namespace lbs

#endif //CHITECH_LBSMIP_SHELL_OPERATIONS_H
