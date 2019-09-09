#include "solver_montecarlon.h"

#include <ChiMesh/Cell/cell.h>
#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

#include <chi_log.h>

extern ChiPhysics chi_physics_handler;
extern ChiLog chi_log;

#include<iostream>
#include<math.h>
#include<typeinfo>

//###################################################################
/**The default raytracing algorithm.*/
void chi_montecarlon::Solver::Raytrace(chi_montecarlon::Particle* prtcl)
{
  //======================================== Get total and scat xs
  auto cell = grid->cells[prtcl->cur_cell_ind];
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];

  chi_physics::TransportCrossSections* xs =
    (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double sigt = xs->sigma_tg[prtcl->egrp];
  double sigs = sigt - xs->sigma_ag[prtcl->egrp];

  if (sigs>sigt)
  {
    chi_log.Log(LOG_ALLERROR)
    << "Error in ray trace material";
    exit(EXIT_FAILURE);
  }


  //======================================== Distance to event
  double d_to_intract = -1.0*log(1.0-rng0.Rand())/sigt;
  double d_to_surface = 1.0e15;

  chi_mesh::Vector posf = prtcl->pos;
  chi_mesh::Vector dirf = prtcl->dir;
  int                ef = prtcl->egrp;
  int         auxinfo[] = {3,-1,-1};
  chi_mesh::RayTrace(grid, cell, prtcl->pos, prtcl->dir,
                     d_to_surface, posf, auxinfo);

  if (isnan(posf.x))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Posf corruption after surface tracking.";
    exit(EXIT_FAILURE);
  }

  if (isnan(d_to_intract))
  {
    chi_log.Log(LOG_ALLERROR)
      << "d_to_interact corrupt.";
    exit(EXIT_FAILURE);
  }
  if (isnan(d_to_surface))
  {
    chi_log.Log(LOG_ALLERROR)
      << "d_to_surface corrupt.";
    exit(EXIT_FAILURE);
  }

//  chi_log.Log(LOG_0)
//  << "d2s=" << d_to_surface << " "
//  << "d2i=" << d_to_intract;

  //======================================== Process interaction
  if (d_to_intract < d_to_surface)
  {
    posf = prtcl->pos + prtcl->dir*d_to_intract;

    if (rng0.Rand() < (sigs/sigt))
    {

      std::pair<int,chi_mesh::Vector> e_dir =
        ProcessScattering(prtcl,xs);
      ef   = e_dir.first;
      dirf = e_dir.second;

      if (mono_energy && (ef != prtcl->egrp))
        prtcl->alive = false;
    }
    else
      prtcl->alive = false;

    ContributeTally(prtcl,posf);
  }
  //======================================== Process surface
  else
  {
    ContributeTally(prtcl,posf);
    if (auxinfo[2] < 0)
    {
      bool reflecting = false;
      if (!reflecting)
        prtcl->alive = false;
      else
      {

      }
    }//if bndry
    else
      prtcl->cur_cell_ind = auxinfo[2];
  }

  prtcl->pos = posf;
  prtcl->dir = dirf;
  prtcl->egrp = ef;

}
