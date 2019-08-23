#ifndef _montecarlon_source_h
#define _montecarlon_source_h

#include"../chi_montecarlon.h"
#include "mc_base_source.h"
#include"../chi_montecarlon_particle.h"

#include <CHI_MESH/chi_mesh.h>
#include <CHI_DISCRETIZATION_FV/fv.h>

#define MC_BASE_SRC     0
#define MC_POINT_SRC    1
#define MC_BNDRY_SRC    2
#define    MC_ALL_BOUNDARIES -1
#define MC_LOGICVOL_SRC 3

//######################################################### Class Def
/**Parent Monte carlo source.
This source is an isotropic source at [0 0 0] with energy of 4 MeV.*/
class chi_montecarlon::Source
{
public:
  chi_mesh::MeshContinuum* grid;
  CHI_DISCRETIZATION_FV*   fv_sdm;
  int type_index;

public:
  //00
          Source();
  //01
  virtual void Initialize(chi_mesh::MeshContinuum* ref_grid,
                          CHI_DISCRETIZATION_FV*   ref_fv_sdm);
  virtual chi_montecarlon::Particle
          CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);

};

#endif