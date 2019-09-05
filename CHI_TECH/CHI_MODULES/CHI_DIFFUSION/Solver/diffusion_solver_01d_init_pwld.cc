#include "diffusion_solver.h"

#include <CHI_MESH/CHI_CELL/cell_slab.h>
#include <CHI_MESH/CHI_CELL/cell_polygon.h>
#include <CHI_MESH/CHI_CELL/cell_polyhedron.h>




#include <CHI_TIMER/chi_timer.h>

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
int chi_diffusion::Solver::InitializePWLD(bool verbose)
{
  //Right now I am only doing one region at a time.
  //Later I want to support multiple regions with interfaces.
  chi_mesh::Region*     aregion = this->regions.back();
  grid = aregion->volume_mesh_continua.back();


  //================================================== Reorder nodes
  if (verbose)
    chi_log.Log(LOG_0) << "Computing nodal reorderings for PWLD";
  CHI_TIMER t_reorder; t_reorder.Reset();
  this->ReorderNodesPWLD();

  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose)
  {
    chi_log.Log(LOG_0) << "Time taken during nodal reordering "
                       << t_reorder.GetTime()/1000.0;
  }


  //================================================== Initialize field function
  //                                                   if empty
  pwld_phi_local.resize(pwld_local_dof_count);
  if (field_functions.size() == 0)
  {
    chi_physics::FieldFunction* initial_field_function =
      new chi_physics::FieldFunction;
    initial_field_function->text_name = std::string("phi");
    initial_field_function->grid = grid;
    initial_field_function->spatial_discretization = discretization;
    initial_field_function->id = chi_physics_handler.fieldfunc_stack.size();

    initial_field_function->type = FF_SDM_PWLD;
    initial_field_function->num_grps = 1;
    initial_field_function->num_moms = 1;
    initial_field_function->grp = 0;
    initial_field_function->mom = 0;
    initial_field_function->field_vector_local = &pwld_phi_local;
    initial_field_function->local_cell_dof_array_address =
      &pwld_cell_dof_array_address;

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

      cur_ff->type = FF_SDM_PWLD;
      cur_ff->num_grps = 1;
      cur_ff->num_moms = 1;
      cur_ff->grp = 0;
      cur_ff->mom = 0;
      cur_ff->field_vector_local = &pwld_phi_local;
      cur_ff->local_cell_dof_array_address = &pwld_cell_dof_array_address;
    }
  }

  //================================================== Setup timers
  if (verbose)
    chi_log.Log(LOG_0) << "Determining nodal connections";
  CHI_TIMER t_connect; t_connect.Reset();
  double t0 = 0.0;


  //================================================== Initialize nodal DOF
  //                                                   and connection info
  for (int i=0; i<pwld_local_dof_count; i++)
  {
    std::vector<int>* new_node_links = new std::vector<int>;
    nodal_connections.push_back(new_node_links);
    new_node_links = new std::vector<int>;
    nodal_cell_connections.push_back(new_node_links);

    nodal_nnz_in_diag.push_back(0);
    nodal_nnz_off_diag.push_back(1);
  }

  nodal_boundary_numbers.resize(grid->nodes.size(),0);
  int total_nnz = 0;



  //================================================== First pass store pwld view
  int num_loc_cells = grid->local_cell_glob_indices.size();
  int dof_count = 0;
  std::set<int> local_border_cells;
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;

      DIFFUSION_IP_VIEW* ip_view = new DIFFUSION_IP_VIEW;
      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
      pwld_cell_dof_array_address.push_back(dof_count);
      ip_cell_views.push_back(ip_view);

      nodal_nnz_in_diag[dof_count]   = 4;
      nodal_nnz_in_diag[dof_count+1] = 4;

      if (slab_cell->edges[0]>=0)
      {
        if (grid->IsCellLocal(slab_cell->edges[0]))
          nodal_nnz_in_diag[dof_count]  += 2;
        else
        {
          nodal_nnz_off_diag[dof_count] += 2;
          local_border_cells.insert(lc);
        }
      }
      else
        nodal_boundary_numbers[slab_cell->v_indices[0]] = slab_cell->edges[0];

      if (slab_cell->edges[1]>=0)
      {
        if (grid->IsCellLocal(slab_cell->edges[1]))
          nodal_nnz_in_diag[dof_count+1]  += 2;
        else
        {
          nodal_nnz_off_diag[dof_count+1] += 2;
          local_border_cells.insert(lc);
        }
      }
      else
        nodal_boundary_numbers[slab_cell->v_indices[1]] = slab_cell->edges[1];

      dof_count += 2;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;

      DIFFUSION_IP_VIEW* ip_view = new DIFFUSION_IP_VIEW;
      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
      pwld_cell_dof_array_address.push_back(dof_count);
      ip_cell_views.push_back(ip_view);

      for (int v=0; v<poly_cell->v_indices.size(); v++)
      {
        nodal_nnz_in_diag[dof_count] = poly_cell->v_indices.size();

        total_nnz+=poly_cell->v_indices.size();
        for (int e=0; e<poly_cell->edges.size(); e++)
        {
          if (poly_cell->edges[e][2]>=0) //Not boundary
          {
            if (grid->IsCellLocal(poly_cell->edges[e][2]))
            {
              int adj_cell_glob_index = poly_cell->edges[e][2];
              auto adj_cell =
                (chi_mesh::CellPolygon*)grid->cells[adj_cell_glob_index];
              nodal_nnz_in_diag[dof_count] += adj_cell->v_indices.size();
            }
            else
            {
              //Since we have no information about the non-local cell,
              //we can make a good assumption that it has the same amount
              //of dofs than does the current cell.
              //nodal_nnz_off_diag[dof_count] += poly_cell->v_indices.size()/2;
              local_border_cells.insert(lc);
            }
          }//if not bndry
        }//for edge
        dof_count++;
      }//for vi

      //==================================== Boundary numbers
      for (int e=0; e<poly_cell->edges.size(); e++)
      {
        if (poly_cell->edges[e][2]<0)
        {
          nodal_boundary_numbers[poly_cell->edges[e][0]]=
            poly_cell->edges[e][2];
          nodal_boundary_numbers[poly_cell->edges[e][1]]=
            poly_cell->edges[e][2];
        }
      }

    }//polygon

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;

      DIFFUSION_IP_VIEW* ip_view = new DIFFUSION_IP_VIEW;
      ip_view->cell_dof_start = dof_count + pwld_local_dof_start;
      pwld_cell_dof_array_address.push_back(dof_count);
      ip_cell_views.push_back(ip_view);

      for (int v=0; v<polyh_cell->v_indices.size(); v++)
      {
        nodal_nnz_in_diag[dof_count] = polyh_cell->v_indices.size();

        for (int f=0; f<polyh_cell->faces.size(); f++)
        {
          if (polyh_cell->faces[f]->face_indices[NEIGHBOR]>=0) //Not bndry
          {
            bool is_local =
              grid->IsCellLocal(polyh_cell->faces[f]->face_indices[NEIGHBOR]);

            if (is_local)
            {
              int adj_cell_glob_index =
                polyh_cell->faces[f]->face_indices[NEIGHBOR];
              auto adj_cell =
                (chi_mesh::CellPolyhedron*)grid->cells[adj_cell_glob_index];
              nodal_nnz_in_diag[dof_count] += adj_cell->v_indices.size();
            }
            else
            {
              //Since we have no information about the non-local cell,
              //we can make a good assumption that it has the same amount
              //of dofs than does the current cell.
              nodal_nnz_off_diag[dof_count] += polyh_cell->v_indices.size();
              local_border_cells.insert(lc);
            }
          }
        }
        dof_count++;
      }

      //==================================== Boundary numbers
      for (int f=0; f<polyh_cell->faces.size(); f++)
      {
        if (polyh_cell->faces[f]->face_indices[NEIGHBOR]<0)
        {
          for (int fv=0; fv<polyh_cell->faces[f]->v_indices.size(); fv++)
          {
            int fvi = polyh_cell->faces[f]->v_indices[fv];
            nodal_boundary_numbers[fvi] =
              polyh_cell->faces[f]->face_indices[NEIGHBOR];
          }//for fv
        }//if bndry
      }//for face v's
    }//if polyhedron

  }//for local cell
  if (verbose)
  {
    chi_log.Log(LOG_0) << "Time taken during nodal connection "
                       << t_connect.GetTime()/1000.0;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (verbose)
    chi_log.Log(LOG_0) << "Communicating border cell information.";
  //================================================== Serialize local cells
  // The vectorized values will be as follows
  // - cell_glob_index
  // - cell_dof_start
  // - cell_type
  // - cell_mat_id
  // - cell_dof_count
  // - cell_face_count
  //
  // - dof 0 glob_index
  //     to
  // - dof N glob_index
  //
  // - face_0 dof_count
  // - face_0 dof 0 glob_index
  //     to
  // - face_0 dof fN glob_index
  //
  // - repeat all face info
  std::vector<int> border_cell_info;

  //============================================= Loop over set
  std::set<int>::iterator local_cell;
  for (local_cell  = local_border_cells.begin();
       local_cell != local_border_cells.end();
       local_cell++)
  {
    int local_cell_index = *local_cell;
    int cell_glob_index = grid->local_cell_glob_indices[local_cell_index];

    auto cell = grid->cells[cell_glob_index];
    DIFFUSION_IP_VIEW* ip_view = ip_cell_views[local_cell_index];

    border_cell_info.push_back(cell_glob_index);         //cell_glob_index
    border_cell_info.push_back(ip_view->cell_dof_start); //cell_dof_start

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (typeid(*cell) == typeid(chi_mesh::CellSlab))
    {
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      border_cell_info.push_back(0);                     //cell_type
      border_cell_info.push_back(slab_cell->material_id);//cell_mat_id
      border_cell_info.push_back(2);                     //cell_dof_count
      border_cell_info.push_back(2);                     //cell_face_count

      for (int v=0; v<2; v++)
        border_cell_info.push_back(slab_cell->v_indices[v]); //dof 0 to N

      for (int f=0; f<2; f++)
      {
        border_cell_info.push_back(1);                   //face dof_count
        border_cell_info.push_back(slab_cell->v_indices[f]); //face dof 0 to fN
      }
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
    {
      auto poly_cell = (chi_mesh::CellPolygon*)cell;
      border_cell_info.push_back(1);                     //cell_type
      border_cell_info.push_back(poly_cell->material_id);//cell_mat_id
      border_cell_info.push_back(poly_cell->v_indices.size());//cell_dof_count
      border_cell_info.push_back(poly_cell->edges.size());//cell_face_count

      for (int v=0; v<poly_cell->v_indices.size(); v++)
        border_cell_info.push_back(poly_cell->v_indices[v]);//dof 0 to N

      for (int f=0; f<poly_cell->edges.size(); f++)
      {
        border_cell_info.push_back(2);                       //face dof_count
        border_cell_info.push_back(poly_cell->edges[f][0]); //face dof 0 to fN
        border_cell_info.push_back(poly_cell->edges[f][1]); //face dof 0 to fN
      }
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      border_cell_info.push_back(2);                     //cell_type
      border_cell_info.push_back(polyh_cell->material_id);//cell_mat_id
      border_cell_info.push_back(polyh_cell->v_indices.size());//cell_dof_count
      border_cell_info.push_back(polyh_cell->faces.size());//cell_face_count

      for (int v=0; v<polyh_cell->v_indices.size(); v++)
        border_cell_info.push_back(polyh_cell->v_indices[v]);//dof 0 to N

      for (int f=0; f<polyh_cell->faces.size(); f++)
      {
        int face_dof_count = polyh_cell->faces[f]->v_indices.size();
        border_cell_info.push_back(face_dof_count);         //face dof_count
        for (int fv=0; fv<face_dof_count; fv++)
          border_cell_info.push_back(polyh_cell->faces[f]->v_indices[fv]);
                                                            //face dof 0 to fN
      }
    }
  }//for local cell

  //================================================== Distribute border info
  std::vector<int> locI_info_size;
  std::vector<std::vector<int>> locI_border_cell_info;

  locI_info_size.resize(chi_mpi.process_count);
  locI_border_cell_info.resize(chi_mpi.process_count);

  //======================================== Collect sizes
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    if (locI == chi_mpi.location_id)
    {
      locI_info_size[locI] = border_cell_info.size();
    }
    MPI_Bcast(&locI_info_size[locI],1,MPI_INT,locI,MPI_COMM_WORLD);
  }

  //======================================== Collect info
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    if (locI == chi_mpi.location_id)
    {
      std::copy(border_cell_info.begin(),
                border_cell_info.end(),
                std::back_inserter(locI_border_cell_info[locI]));
    }
    else
      locI_border_cell_info[locI].resize(locI_info_size[locI]);

    MPI_Bcast(locI_border_cell_info[locI].data(),
              locI_info_size[locI],MPI_INT,locI,MPI_COMM_WORLD);
  }

  //================================================== Deserialize border info
  // The vectorized values will be as follows
  // - cell_glob_index
  // - cell_dof_start
  // - cell_type
  // - cell_mat_id
  // - cell_dof_count
  // - cell_face_count
  //
  // - dof 0 glob_index
  //     to
  // - dof N glob_index
  //
  // - face_0 dof_count
  // - face_0 dof 0 glob_index
  //     to
  // - face_0 dof fN glob_index
  //
  // - repeat all face info
  ip_locI_bordercell_info.resize(chi_mpi.process_count);
  ip_locI_bordercells.resize(chi_mpi.process_count);
  ip_locI_borderfeviews.resize(chi_mpi.process_count);
  ip_locI_borderipviews.resize(chi_mpi.process_count);
  for (int locI=0; locI<chi_mpi.process_count; locI++)
  {
    int k=0;
    while (k<locI_info_size[locI])
    {
      DIFFUSION_IP_BORDERCELL* border_cell = new DIFFUSION_IP_BORDERCELL;
      border_cell->cell_glob_index = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_dof_start  = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_type       = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_mat_id     = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_dof_count  = locI_border_cell_info[locI][k]; k++;
      border_cell->cell_face_count = locI_border_cell_info[locI][k]; k++;

      int dof_count = border_cell->cell_dof_count;

      for (int v=0; v<dof_count; v++)
      {
        border_cell->v_indices.push_back(locI_border_cell_info[locI][k]);
        k++;
      }

      int face_count = border_cell->cell_face_count;

      for (int f=0; f<face_count; f++)
      {
        int face_dof_count = locI_border_cell_info[locI][k]; k++;
        border_cell->face_v_indices.push_back(std::vector<int>());
        for (int fv=0; fv<face_dof_count; fv++)
        {
          int vgi = locI_border_cell_info[locI][k]; k++;
          border_cell->face_v_indices[f].push_back(vgi);
        }
      }

      ip_locI_bordercell_info[locI].push_back(border_cell);
    }//while less than buffersize

    int locI_num_bordercells  = ip_locI_bordercell_info[locI].size();
    ip_locI_bordercells[locI].resize(locI_num_bordercells,nullptr);
    ip_locI_borderfeviews[locI].resize(locI_num_bordercells,nullptr);
    ip_locI_borderipviews[locI].resize(locI_num_bordercells,nullptr);
  }


  //================================================== Initialize x and b
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,pwld_local_dof_count,
                     pwld_global_dof_count);CHKERRQ(ierr);
  ierr = VecSetType(x,VECMPI);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  VecSet(x,0.0);
  VecSet(b,0.0);

  //################################################## Create matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,pwld_local_dof_count,
                     pwld_local_dof_count,
                     pwld_global_dof_count,pwld_global_dof_count);CHKERRQ(ierr);
  ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);


  //================================================== Allocate matrix memory
  if (verbose)
  {
    chi_log.Log(LOG_0)
      << "Setting matrix preallocation. Total non-zeros: " << total_nnz;
  }


  MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                              0,nodal_nnz_off_diag.data());
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
  MatSetUp(A);


  //================================================== Set up solver
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetOperators(ksp,A,A);
  ierr = KSPSetType(ksp,KSPCG);

  ierr = KSPGetPC(ksp,&pc);
  PCSetType(pc,PCHYPRE);

  PCHYPRESetType(pc,"boomeramg");

  //================================================== Setting Hypre parameters
  //The default HYPRE parameters used for polyhedra
  //seemed to have caused a lot of trouble for Slab
  //geometries. This section makes some custom options
  //per cell type
  int first_cell_g_index = grid->local_cell_glob_indices[0];
  auto first_cell = grid->cells[first_cell_g_index];

  if (typeid(*first_cell) == typeid(chi_mesh::CellSlab))
  {
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_P_max 4");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");
  }

  if (typeid(*first_cell) == typeid(chi_mesh::CellPolyhedron))
  {
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_strong_threshold 0.8");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_P_max 4");

    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(NULL,"-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(NULL,"-options_left");

  }
  PetscOptionsInsertString(NULL,options_string.c_str());
  PCSetFromOptions(pc);

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