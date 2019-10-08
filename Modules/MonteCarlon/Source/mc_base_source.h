#ifndef _montecarlon_source_h
#define _montecarlon_source_h

#include"../chi_montecarlon.h"
#include "mc_base_source.h"
#include"../chi_montecarlon_particle.h"

#include <ChiMesh/chi_mesh.h>
#include <FiniteVolume/fv.h>

#define MC_BASE_SRC     0
#define MC_POINT_SRC    1
#define MC_BNDRY_SRC    2
#define    MC_ALL_BOUNDARIES -1
#define MC_LOGICVOL_SRC 3
#define MC_RESID_SRC    4
#define MC_RESID_SRC_SU 5
#define MC_RESID_MOC    6
#define MC_RESID_MOC_SU 7

//######################################################### Class Def
/**Parent Monte carlo source.
This source is an isotropic source at [0 0 0] with energy of 4 MeV.*/
class chi_montecarlon::Source
{
public:
  int particles_C;
  int particles_L;
  int particles_R;
  double weights_L;
  double weights_R;

public:
  chi_mesh::MeshContinuum* grid;
  SpatialDiscretization_FV*   fv_sdm;
  int type_index;

public:
  //00
          Source();
  //01
  virtual void Initialize(chi_mesh::MeshContinuum* ref_grid,
                          SpatialDiscretization_FV*   ref_fv_sdm);
  virtual chi_montecarlon::Particle
          CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);

};

#endif