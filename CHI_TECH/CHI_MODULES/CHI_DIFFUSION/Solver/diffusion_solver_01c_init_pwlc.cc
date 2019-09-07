#include "diffusion_solver.h"

#include <CHI_MESH/CHI_CELL/cell_slab.h>
#include <CHI_MESH/CHI_CELL/cell_triangle.h>
#include <CHI_MESH/CHI_CELL/cell_polygon.h>
#include <CHI_MESH/CHI_CELL/cell_polyhedron.h>
#include <CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h>

#include <ChiTimer/chi_timer.h>

#include <chi_mpi.h>
#include <chi_log.h>
#include <CHI_PHYSICS/chi_physics.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;
extern CHI_PHYSICS chi_physics_handler;

#include<fstream>
#include <unistd.h>

PetscErrorCode
DiffusionConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                            KSPConvergedReason* convergedReason,
                            void *monitordestroy);

//###################################################################
/**Initializes Piecewise Linear FEM for diffusion solver.*/
int chi_diffusion::Solver::InitializePWLC(bool verbose)
{
  //Right now I am only doing one region at a time.
  //Later I want to support multiple regions with interfaces.
  chi_mesh::Region*              aregion = this->regions.back();
  grid = aregion->volume_mesh_continua.back();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  int n = grid->nodes.size();

  //==================================================
  chi_log.Log(LOG_0) << "Computing nodal reorderings for CFEM";
  ChiTimer t_reorder; t_reorder.Reset();
  this->ReorderNodesPWLC();


  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0) << "Time taken during nodal reordering "
                     << t_reorder.GetTime()/1000.0;


  //================================================== Initialize x and b
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,local_rows_to - local_rows_from+1,n);CHKERRQ(ierr);
  ierr = VecSetType(x,VECMPI);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  VecSet(x,0.0);
  VecSet(b,0.0);



  //================================================== Initialize field function
  //                                                   if empty
  if (field_functions.size() == 0)
  {
    chi_physics::FieldFunction* initial_field_function =
      new chi_physics::FieldFunction;
    initial_field_function->text_name = std::string("phi");
    initial_field_function->grid = grid;
    initial_field_function->spatial_discretization = discretization;
    initial_field_function->id = chi_physics_handler.fieldfunc_stack.size();

    initial_field_function->field_vector = x;

    field_functions.push_back(initial_field_function);
    chi_physics_handler.fieldfunc_stack.push_back(initial_field_function);
  }
  else
  {
    size_t num_ff = field_functions.size();
    for (int ff=0; ff<num_ff; ff++)
    {
      chi_physics::FieldFunction* cur_ff = field_functions[ff];
      cur_ff->grid                   = grid;
      cur_ff->spatial_discretization = discretization;
      cur_ff->field_vector           = x;
    }
  }



  //################################################## Create matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,local_rows_to - local_rows_from+1,
                     local_rows_to - local_rows_from+1,
                     n,n);CHKERRQ(ierr);
  ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  //================================================== Setup timer
  chi_log.Log(LOG_0) << "Determining nodal connections";
  ChiTimer t_connect; t_connect.Reset();
  double t0 = 0.0;


  //================================================== Initialize nodal DOF
  //                                                   and connection info
  for (int i=0; i<grid->nodes.size(); i++)
  {
    std::vector<int>* new_node_links = new std::vector<int>;
    nodal_connections.push_back(new_node_links);
    new_node_links = new std::vector<int>;
    nodal_cell_connections.push_back(new_node_links);

    nodal_boundary_numbers.push_back(0);
    nodal_nnz_in_diag.push_back(0);
    nodal_nnz_off_diag.push_back(0);
  }



  //================================================== Determine nodal DOF
  size_t num_local_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc < num_local_cells; lc++)
  {
    int glob_index = grid->local_cell_glob_indices[lc];

    auto cell = grid->cells[glob_index];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      chi_mesh::CellSlab* slab_cell =
        (chi_mesh::CellSlab*)(cell);
      for (int i=0; i<2; i++)
      {
        int ir =  mesher->MapNode(slab_cell->v_indices[i]);

        //================================== Check if i is on boundary
        for (int e=0; e<2; e++)
        {
          if (slab_cell->edges[e]<0)
          {
            int v0_index =
              mesher->MapNode(slab_cell->v_indices[e]);
            if (ir == v0_index)
            {
              int boundary_type =
                boundaries[abs(slab_cell->edges[e])-1]->type;
              if (boundary_type == DIFFUSION_DIRICHLET)
              {
                nodal_boundary_numbers[ir]=slab_cell->edges[e];
              }
              break;
            } //if ir part of face
          }//if boundary
        }//for edge

        //======================================= Set nodal connections
        std::vector<int>* node_links = nodal_connections[ir];
        for (int j=0; j<2; j++)
        {
          int jr = mesher->MapNode(slab_cell->v_indices[j]);

          //====================== Check for duplicates
          bool already_there = false;
          for (int k=0; k<node_links->size(); k++)
          {
            if ((*node_links)[k] == jr)
            {already_there = true; break;}
          }
          if (!already_there)
          {
            (*node_links).push_back(jr);
            if ((jr>=local_rows_from) && (jr<=local_rows_to))
            {
              nodal_nnz_in_diag[ir]+=1;
            } else
            {
              nodal_nnz_off_diag[ir]+=1;
            }
          }
        }//for j
      }//for i
    }//if slab


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      chi_mesh::CellPolygon* poly_cell =
        (chi_mesh::CellPolygon*)(cell);
      for (int i=0; i<poly_cell->v_indices.size(); i++)
      {
        int ir =  mesher->MapNode(poly_cell->v_indices[i]);

        //================================== Check if i is on boundary
        for (int e=0; e<poly_cell->edges.size(); e++)
        {
          if (poly_cell->edges[e][2]<0)
          {
            int v0_index =
              mesher->MapNode(poly_cell->edges[e][0]);
            int v1_index =
              mesher->MapNode(poly_cell->edges[e][1]);
            if ((ir == v0_index) || (ir == v1_index))
            {
              int boundary_type =
                boundaries[abs(poly_cell->edges[e][2])-1]->type;
              if (boundary_type == DIFFUSION_DIRICHLET)
              {
                nodal_boundary_numbers[ir]=poly_cell->edges[e][2];
              }
              break;
            } //if ir part of face
          }
        }

        //======================================= Set nodal connections
        std::vector<int>* node_links = nodal_connections[ir];
        for (int j=0; j<poly_cell->v_indices.size(); j++)
        {
          int jr = mesher->MapNode(poly_cell->v_indices[j]);

          //====================== Check for duplicates
          bool already_there = false;
          for (int k=0; k<node_links->size(); k++)
          {
            if ((*node_links)[k] == jr)
            {already_there = true; break;}
          }
          if (!already_there)
          {
            (*node_links).push_back(jr);
            if ((jr>=local_rows_from) && (jr<=local_rows_to))
            {
              nodal_nnz_in_diag[ir]+=1;
            } else
            {
              nodal_nnz_off_diag[ir]+=1;
            }
          }
        }//for j
      }//for vertices i
    }//if Polygon
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% If polyhedral item_id
    else if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      chi_mesh::CellPolyhedron* polyh_cell =
        (chi_mesh::CellPolyhedron*)(cell);

      for (int i=0; i<polyh_cell->v_indices.size(); i++)
      {
        int ir =  mesher->MapNode(polyh_cell->v_indices[i]);

        if (ir<0)
        {
          chi_log.Log(LOG_ALLERROR)
            << "ir Mapping error node " << polyh_cell->v_indices[i];
          exit(EXIT_FAILURE);
        }

        //================================== Check if i is on boundary
        for (int f=0; f<polyh_cell->faces.size(); f++)
        {
          if (polyh_cell->faces[f]->face_indices[NEIGHBOR] < 0)
          {
            for (int e=0; e<polyh_cell->faces[f]->edges.size(); e++)
            {
              int v0_index =
                mesher->MapNode(polyh_cell->faces[f]->edges[e][0]);

              if (v0_index<0)
              {
                chi_log.Log(LOG_ALLERROR)
                  << "v0 Mapping error node " << polyh_cell->faces[f]->edges[e][0];
                exit(EXIT_FAILURE);
              }

              if ((ir == v0_index))
              {
                //================= Processing boundary
                int boundary_type =
                  boundaries[abs(polyh_cell->faces[f]->face_indices[NEIGHBOR])-1]->type;
                if (boundary_type == DIFFUSION_DIRICHLET)
                {
                  nodal_boundary_numbers[ir]=
                    polyh_cell->faces[f]->face_indices[0];
                }
                else if (polyh_cell->faces[f]->face_indices[0]>=0)
                {
                  int adj_cell_index = polyh_cell->faces[f]->face_indices[0];
                  auto adj_cell = grid->cells[adj_cell_index];
                  if (adj_cell->partition_id != chi_mpi.location_id)
                  {
                    nodal_boundary_numbers[ir]=PROCESS_BOUNDARY;
                  }
                }
                break;
              } //if ir part of face
            }
          }
        }//for f

        //======================================= Set nodal connections
        std::vector<int>* node_links = nodal_connections[ir];
        for (int j=0; j<polyh_cell->v_indices.size(); j++)
        {
          int jr = mesher->MapNode(polyh_cell->v_indices[j]);

          //====================== Check for duplicates
          bool already_there = false;
          for (int k=0; k<node_links->size(); k++)
          {
            if ((*node_links)[k] == jr)
            {already_there = true; break;}
          }
          if (!already_there)
          {
            (*node_links).push_back(jr);
            if ((jr>=local_rows_from) && (jr<=local_rows_to))
            {
              nodal_nnz_in_diag[ir]+=1;
            } else
            {
              nodal_nnz_off_diag[ir]+=1;
            }
          }
        }//for j
      }//for i

    } //if typeid

  }
  chi_log.Log(LOG_0) << "Time taken during nodal connection "
                     << t_connect.GetTime()/1000.0;

  //================================================== Allocate matrix memory
  chi_log.Log(LOG_0) << "Setting matrix preallocation.";
  MatMPIAIJSetPreallocation(A,0,&nodal_nnz_in_diag[local_rows_from],
                            0,&nodal_nnz_off_diag[local_rows_from]);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
  MatSetUp(A);

  //================================================== Set up solver
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetOperators(ksp,A,A);
  ierr = KSPSetType(ksp,KSPCG);

  ierr = KSPGetPC(ksp,&pc);
  ierr = PCSetType(pc,PCGAMG);

  PetscOptionsInsertString(NULL,options_string.c_str());

  //=================================== Set up monitor
  if (verbose)
    ierr = KSPMonitorSet(ksp,&chi_diffusion::KSPMonitorAChiTech,NULL,NULL);

  KSPSetConvergenceTest(ksp,&DiffusionConvergenceTestNPT,NULL,NULL);

  //=================================== Setup verbose viewer
  if (chi_log.GetVerbosity()>= LOG_0VERBOSE_2)
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

  ierr = KSPSetTolerances(ksp,1.e-50,residual_tolerance,1.0e50,max_iters);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  return 0;
}