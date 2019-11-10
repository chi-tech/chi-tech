#include "solver_montecarlon.h"

#include <ChiMesh/Cell/cell.h>
#include <ChiMesh/Raytrace/raytracing.h>
#include <ChiPhysics/chi_physics.h>
#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

#include <chi_log.h>

extern ChiPhysics chi_physics_handler;
extern ChiLog chi_log;

#include<iostream>
#include<cmath>

//###################################################################
/**The default raytracing algorithm.*/
void chi_montecarlon::Solver::Raytrace(chi_montecarlon::Particle* prtcl)
{
  //======================================== Get total and scat xs
  auto cell = grid->cells[prtcl->cur_cell_ind];
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];

  auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double sigt = xs->sigma_tg[prtcl->egrp];
  double sigs = sigt - xs->sigma_ag[prtcl->egrp];

  //======================================== Distance to event
  double d_to_intract = -1.0*log(1.0-rng0.Rand())/sigt;
  double d_to_surface = 1.0e15;

  chi_mesh::Vector posf = prtcl->pos;
  chi_mesh::Vector dirf = prtcl->dir;
  int                ef = prtcl->egrp;
  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(grid, cell, prtcl->pos, prtcl->dir,
                       d_to_surface, posf);


  //======================================== Process interaction
  if (d_to_intract < d_to_surface)
  {
    posf = prtcl->pos + prtcl->dir*d_to_intract;

    if (rng0.Rand() < (sigs/sigt))
    {

      auto energy_dir = ProcessScattering(prtcl,xs);
      ef   = energy_dir.first;
      dirf = energy_dir.second;
//      prtcl->alive = false;  // TODO: Remove before flight

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
    //posf set in call to RayTrace
    ContributeTally(prtcl,posf);
    if (ray_dest_info.destination_face_neighbor < 0)
    {
      bool reflecting = false;
      if (!reflecting)
        prtcl->alive = false;
      else
      {

      }
    }//if bndry
    else
      prtcl->cur_cell_ind = ray_dest_info.destination_face_neighbor;
  }

  prtcl->pos = posf;
  prtcl->dir = dirf;
  prtcl->egrp = ef;

}
