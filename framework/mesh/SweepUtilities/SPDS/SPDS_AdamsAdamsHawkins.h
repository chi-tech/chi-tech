#ifndef CHITECH_SPDS_ADAMSADAMSHAWKINS_H
#define CHITECH_SPDS_ADAMSADAMSHAWKINS_H

#include "SPDS.h"

namespace chi_mesh::sweep_management
{

class SPDS_AdamsAdamsHawkins : public SPDS
{
public:
  SPDS_AdamsAdamsHawkins(const chi_mesh::Vector3& omega,
                         const chi_mesh::MeshContinuum& grid,
                         bool cycle_allowance_flag);
};

}

#endif // CHITECH_SPDS_ADAMSADAMSHAWKINS_H
