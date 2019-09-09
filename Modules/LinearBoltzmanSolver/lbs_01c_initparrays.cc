#include "lbs_linear_boltzman_solver.h"
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>
#include <PiecewiseLinear/pwl.h>
#include <PiecewiseLinear/CellViews/pwl_slab.h>
#include <PiecewiseLinear/CellViews/pwl_polygon.h>
#include <PiecewiseLinear/CellViews/pwl_polyhedron.h>
#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

#include <iomanip>
#include "ChiConsole/chi_console.h"

extern ChiConsole chi_console;

//###################################################################
/**Initializes p_arrays.\n
The question arises of what datatype can store the total amount of
unknowns. For now we will say we want to be
designing for 100 billion cells with
an assumed shape of a truncated octahedron which has 24 vertices.
We will also assume that we will be able to do 2000 energy groups
and finally we will assume we will do scattering orders up to 16
which requires 289 moments.
   DOFS per truncated octahedron = 24\n
   Energy groups                 = 2000\n
   Moments                       = 289\n
   # of cells                    =   100,000,000,000\n
   Total DOFS                    = 2,400,000,000,000\n
   Unknowns per cell             =        13,872,000\n
   Total Unknowns                = A crap ton\n
\n
It is easy to see here that this is a hell of a lot so how about we think about
something more modest. Like 200 energy groups scattering order 5 (36 moments)
and 2 billion cells.\n
   Energy groups                 = 200\n
   Moments                       = 36\n
   # of cells                    =     2,000,000,000\n
   Total DOFS                    =    48,000,000,000\n
   Unknowns per cell             =             7,200\n
   Total Unknowns                = 1.44xe13\n
\n
A long int only supports up to 4.29e9. This obviously requires
unsigned long long int which can hold up to 2x2e63.\n
\n
Another interesting aspect is what it will take to get to exascale. For a
discrete ordinates this will undoubtly be evident in the amount of angular flux
unknowns. 1 billion cells, 24 vertices, 200 groups, 48 azimuthal angles per
octant, 8 polar angles per octant (3072) angles. 1.47456e16. Just a factor 68
away from exascale.
   */
int LinearBoltzmanSolver::InitializeParrays()
{
//  PetscErrorCode ierr;
//  PetscMPIInt    size;

  SpatialDiscretization_PWL* pwl_discretization =
    (SpatialDiscretization_PWL*)discretization;

  //================================================== Compute local # of dof
  local_dof_count=0;
  if (typeid(*discretization) == typeid(SpatialDiscretization_PWL))
  {
    size_t num_cell_views = pwl_discretization->cell_fe_views.size();
    for (int c=0; c<num_cell_views; c++)
    {
      local_dof_count += pwl_discretization->cell_fe_views[c]->dofs;
    }
  }
  chi_log.Log(LOG_ALLVERBOSE_2) << "Local DOF count = " << local_dof_count;

  //================================================== Compute global # of dof
  MPI_Allreduce(&local_dof_count,&glob_dof_count,1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM,MPI_COMM_WORLD);
  chi_log.Log(LOG_ALLVERBOSE_2) << "Global DOF count = " << glob_dof_count;

  //================================================== Compute num of unknowns
  int G = groups.size();
  int M = num_moments;
  unsigned long long local_unknown_count = local_dof_count*G*M;
  unsigned long long glob_unknown_count  = glob_dof_count*G*M;

  chi_log.Log(LOG_ALLVERBOSE_2) << "Local Unknown count = "
                                      << local_unknown_count;
  chi_log.Log(LOG_ALLVERBOSE_2) << "Globl Unknown count = "
                                      << glob_unknown_count;

  //================================================== Size fixed source
//  q_fixed_local.resize(local_unknown_count,0.0);
  q_moments_local.resize(local_unknown_count,0.0);
  phi_old_local.resize(local_unknown_count,0.0);
  phi_new_local.resize(local_unknown_count,0.0);

//  //================================================== Initialize local indices
//  for (int i=0; i<local_unknown_count; i++)
//    local_indices.push_back(i);

  //================================================== Set local to global mapping
  //This has been deferred to solvers

  //================================================== Initialize default
  //                                                   incident boundary
  std::vector<double> zero_boundary(G,0.0);
  incident_P0_mg_boundaries.push_back(zero_boundary);

  //================================================== Initialize transport views
  int num_grps = groups.size();
  int block_MG_counter = 0;       //Counts the strides of moment and group
  int block_counter = 0;          //Counts the base stride

  chi_mesh::Vector ihat(1.0,0.0,0.0);
  chi_mesh::Vector jhat(0.0,1.0,0.0);
  chi_mesh::Vector khat(0.0,0.0,1.0);
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell = grid->cells[cell_g_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)cell;
      SlabFEView* slab_fe_view =
        (SlabFEView*)pwl_discretization->MapFeView(cell_g_index);
      LBSCellViewFull* full_cell_view =
        (LBSCellViewFull*)cell_transport_views[slab_cell->cell_local_id];

      int mat_id = cell->material_id;

      full_cell_view->xs_id = matid_to_xs_map[mat_id];

      full_cell_view->dof_phi_map_start = block_MG_counter;
      block_MG_counter += slab_fe_view->dofs*num_grps*num_moments;

      //Init face upwind flags and adj_partition_id
      int num_faces = 2;
      full_cell_view->face_f_upwind_flag.resize(num_faces,false);
      for (int f=0; f<num_faces; f++)
      {
        if (slab_cell->edges[f] >= 0)
        {
          int adj_g_index = slab_cell->edges[f];
          auto adj_cell = grid->cells[adj_g_index];

          full_cell_view->face_f_adj_part_id.push_back(
            adj_cell->partition_id);
        }//if not bndry
        else
        {
          full_cell_view->face_f_adj_part_id.push_back(
            slab_cell->edges[f]);

          chi_mesh::Vector& face_norm = slab_cell->face_normals[f];

          if      (face_norm.Dot(ihat)>0.999)
            full_cell_view->face_boundary_id.push_back(0);
          else if (face_norm.Dot(ihat)<-0.999)
            full_cell_view->face_boundary_id.push_back(1);
          else if (face_norm.Dot(jhat)>0.999)
            full_cell_view->face_boundary_id.push_back(2);
          else if (face_norm.Dot(jhat)<-0.999)
            full_cell_view->face_boundary_id.push_back(3);
          else if (face_norm.Dot(khat)>0.999)
            full_cell_view->face_boundary_id.push_back(4);
          else if (face_norm.Dot(khat)<-0.999)
            full_cell_view->face_boundary_id.push_back(5);
        }//if bndry
      }//for f

      //Add address
      local_cell_phi_dof_array_address.push_back(full_cell_view->dof_phi_map_start);
      local_cell_dof_array_address.push_back(block_counter);
      block_counter += slab_fe_view->dofs;
    }//slab

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)cell;
      PolygonFEView* poly_fe_view =
        (PolygonFEView*)pwl_discretization->MapFeView(cell_g_index);
      LBSCellViewFull* full_cell_view =
        (LBSCellViewFull*)cell_transport_views[poly_cell->cell_local_id];

      int mat_id = cell->material_id;

      full_cell_view->xs_id = matid_to_xs_map[mat_id];

      full_cell_view->dof_phi_map_start = block_MG_counter;
      block_MG_counter += poly_fe_view->dofs*num_grps*num_moments;

      //Init face upwind flags and adj_partition_id
      int num_faces = poly_cell->edges.size();
      full_cell_view->face_f_upwind_flag.resize(num_faces,false);
      for (int f=0; f<num_faces; f++)
      {
        if (poly_cell->edges[f][EDGE_NEIGHBOR] >= 0)
        {
          int adj_g_index = poly_cell->edges[f][EDGE_NEIGHBOR];
          auto adj_cell = grid->cells[adj_g_index];

          full_cell_view->face_f_adj_part_id.push_back(
            adj_cell->partition_id);
        }//if not bndry
        else
        {
          full_cell_view->face_f_adj_part_id.push_back(
            poly_cell->edges[f][EDGE_NEIGHBOR]);

          chi_mesh::Vector& face_norm = poly_cell->edgenormals[f];

          if      (face_norm.Dot(ihat)>0.999)
            full_cell_view->face_boundary_id.push_back(0);
          else if (face_norm.Dot(ihat)<-0.999)
            full_cell_view->face_boundary_id.push_back(1);
          else if (face_norm.Dot(jhat)>0.999)
            full_cell_view->face_boundary_id.push_back(2);
          else if (face_norm.Dot(jhat)<-0.999)
            full_cell_view->face_boundary_id.push_back(3);
          else if (face_norm.Dot(khat)>0.999)
            full_cell_view->face_boundary_id.push_back(4);
          else if (face_norm.Dot(khat)<-0.999)
            full_cell_view->face_boundary_id.push_back(5);
        }//if bndry
      }//for f

      //Add address
      local_cell_phi_dof_array_address.push_back(full_cell_view->dof_phi_map_start);
      local_cell_dof_array_address.push_back(block_counter);
      block_counter += poly_fe_view->dofs;
    }//polygon

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)cell;
      PolyhedronFEView* polyh_fe_view =
        (PolyhedronFEView*)pwl_discretization->MapFeView(cell_g_index);
      LBSCellViewFull* full_cell_view =
        (LBSCellViewFull*)cell_transport_views[polyh_cell->cell_local_id];

      int mat_id = cell->material_id;

      full_cell_view->xs_id = matid_to_xs_map[mat_id];

      full_cell_view->dof_phi_map_start = block_MG_counter;
      block_MG_counter += polyh_fe_view->dofs*num_grps*num_moments;

      //Init face upwind flags and adj_partition_id
      int num_faces = polyh_cell->faces.size();
      full_cell_view->face_f_upwind_flag.resize(num_faces,false);
      for (int f=0; f<num_faces; f++)
      {
        if (polyh_cell->faces[f]->face_indices[0] >= 0)
        {
          int adj_g_index = polyh_cell->faces[f]->face_indices[0];
          auto adj_cell = grid->cells[adj_g_index];

          full_cell_view->face_f_adj_part_id.push_back(
            adj_cell->partition_id);
        }//if not bndry
        else
        {
          full_cell_view->face_f_adj_part_id.push_back(
            polyh_cell->faces[f]->face_indices[0]);

          chi_mesh::Vector& face_norm = polyh_cell->faces[f]->geometric_normal;

          if      (face_norm.Dot(ihat)>0.999)
            full_cell_view->face_boundary_id.push_back(0);
          else if (face_norm.Dot(ihat)<-0.999)
            full_cell_view->face_boundary_id.push_back(1);
          else if (face_norm.Dot(jhat)>0.999)
            full_cell_view->face_boundary_id.push_back(2);
          else if (face_norm.Dot(jhat)<-0.999)
            full_cell_view->face_boundary_id.push_back(3);
          else if (face_norm.Dot(khat)>0.999)
            full_cell_view->face_boundary_id.push_back(4);
          else if (face_norm.Dot(khat)<-0.999)
            full_cell_view->face_boundary_id.push_back(5);
        }//if bndry
      }//for f

      //Add address
      local_cell_phi_dof_array_address.push_back(full_cell_view->dof_phi_map_start);
      local_cell_dof_array_address.push_back(block_counter);
      block_counter += polyh_fe_view->dofs;
    }//polyhedron
  }//for local cell

  //================================================== Initialize Field Functions
  for (int g=0; g<groups.size(); g++)
  {
    for (int m=0; m<num_moments; m++)
    {
      chi_physics::FieldFunction* group_ff =
        new chi_physics::FieldFunction;
      group_ff->text_name = std::string("Flux_g") +
                            std::to_string(g) +
                            std::string("_m") + std::to_string(m);
      group_ff->grid = grid;
      group_ff->spatial_discretization = discretization;
      group_ff->id = chi_physics_handler.fieldfunc_stack.size();

      group_ff->type = FF_SDM_PWLD;
      group_ff->num_grps = groups.size();
      group_ff->num_moms = num_moments;
      group_ff->grp = g;
      group_ff->mom = m;
      group_ff->field_vector_local = &phi_old_local;
      group_ff->local_cell_dof_array_address = &local_cell_phi_dof_array_address;

      chi_physics_handler.fieldfunc_stack.push_back(group_ff);
      field_functions.push_back(group_ff);
    }//for m
  }//for g




  return 0;
}