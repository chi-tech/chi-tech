#include "mc_rmc_source.h"

#include "../../RandomNumberGenerator/montecarlon_rng.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_slab.h>
#include <FiniteVolume/CellViews/fv_polygon.h>

#include <PiecewiseLinear/pwl.h>
#include <PiecewiseLinear/CellViews/pwl_slab.h>
#include <PiecewiseLinear/CellViews/pwl_polygon.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualSource::ResidualSource(
  chi_physics::FieldFunction *in_resid_ff, bool use_uniform_sampling) :
  sample_uniformly(use_uniform_sampling)
{
  resid_ff = in_resid_ff;
  particles_L = 0;
  particles_R = 0;

  weights_L = 0.0;
  weights_R = 0.0;

  total_residual = 0.0;
  total_Cresidual = 0.0;
  total_Lresidual = 0.0;
  total_Rresidual = 0.0;
}

//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void chi_montecarlon::ResidualSource::
  Initialize(chi_mesh::MeshContinuum *ref_grid,
             SpatialDiscretization_FV *ref_fv_sdm)
{
  chi_log.Log(LOG_0) << "Initializing RMC Source";
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;

  //================================================== Assert same grid
  //                                                   for source and
  //                                                   field function
  if (resid_ff->grid != grid)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon: Dissimilar grids encountered in "
      << "call to ResidualSource::Initialize. "
      << " The grid for which the residual will be computed must be"
         " the same as the grid on which the residual source is to be"
         " sampled.";
    exit(EXIT_FAILURE);
  }

  SpatialDiscretization_PWL* resid_sdm_pwl = nullptr;
  if (typeid(*resid_ff->spatial_discretization) ==
      typeid(SpatialDiscretization_PWL))
  {
    resid_sdm_pwl =
      (SpatialDiscretization_PWL*)resid_ff->spatial_discretization;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon: Unsupported spatial discretization encountered in "
      << "call to ResidualSource::Initialize.";
    exit(EXIT_FAILURE);
  }


  typedef chi_physics::TransportCrossSections TrXS;
  typedef chi_physics::IsotropicMultiGrpSource TrQ;

  chi_mesh::FieldFunctionInterpolation ff_interp;
  ff_interp.grid_view = grid;

  std::vector<double>& field = *resid_ff->field_vector_local;


  //================================================== Loop over local cells
  chi_log.Log(LOG_0) << "Computing cell residuals";
  int num_local_cells = grid->local_cell_glob_indices.size();
  std::vector<double> cell_averages(num_local_cells,0.0);
  cell_residuals.resize(num_local_cells,0.0);
  cell_interior_residual.resize(num_local_cells,0.0);
  cell_surface_residualL.resize(num_local_cells,0.0);
  cell_surface_residualR.resize(num_local_cells,0.0);
  for (int lc=0; lc<num_local_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

    chi_physics::Material* cell_mat =
      chi_physics_handler.material_stack[cell->material_id];

    TrXS* xs = nullptr;
    TrQ* source = nullptr;

    for (int p=0; p<cell_mat->properties.size(); p++)
    {
      if (typeid(*cell_mat->properties[p]) == typeid(TrXS))
        xs = (TrXS*)cell_mat->properties[p];

      if (typeid(*cell_mat->properties[p]) == typeid(TrQ))
        source = (TrQ*)cell_mat->properties[p];
    }

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      chi_log.Log(LOG_0VERBOSE_1) << "**************** Cell " << cell_glob_index;
      auto slab_cell = (chi_mesh::CellSlab*)cell;
      auto cell_fe_view =
        static_cast<SlabFEView*>(
        resid_sdm_pwl->MapFeView(cell_glob_index));

      //==================================== Creating current cell dof-mapping
      std::vector<int> dofs_to_map(cell_fe_view->dofs);
      std::vector<int> cells_to_map;
      std::vector<int> cur_cell_mapping;
      for (int i=0; i<cell_fe_view->dofs; i++)
        dofs_to_map[i] = i;
      cells_to_map.resize(cell_fe_view->dofs,cell_glob_index);

      chi_log.Log(LOG_0VERBOSE_1) << "Mapping cell dofs ";
      ff_interp.CreatePWLDMapping(
        resid_ff->num_grps,resid_ff->num_moms,0,0,
        dofs_to_map,cells_to_map,
        *resid_ff->local_cell_dof_array_address,&cur_cell_mapping);

      chi_log.Log(LOG_0VERBOSE_1)
       << "dof 0 phi=" << field[cur_cell_mapping[0]] << "\n"
       << "dof 1 phi=" << field[cur_cell_mapping[1]];

      //==================================== Creating adj cell dof-mapping
      chi_log.Log(LOG_0VERBOSE_1) << "Mapping adj cell dofs";
      std::vector<std::vector<int>> adj_cell_mapping_f;

      int num_faces = 2;
      for (int f=0; f<num_faces; f++)
      {
        int adj_cell_index = slab_cell->edges[f];

        std::vector<int> adj_mapping;

        if (adj_cell_index >= 0)
        {
          auto adj_cell = (chi_mesh::CellSlab*)grid->cells[adj_cell_index];

          int adj_num_dofs = 2;
          for (int i=0; i<adj_num_dofs; i++)
          {
            dofs_to_map[i] = i;
            cells_to_map[i] = adj_cell_index;
          }

          ff_interp.CreatePWLDMapping(
            resid_ff->num_grps,resid_ff->num_moms,0,0,
            dofs_to_map,cells_to_map,
            *resid_ff->local_cell_dof_array_address,&adj_mapping);

          chi_log.Log(LOG_0VERBOSE_1)
            << "adj_cell " << adj_cell_index << "\n"
            << "dof 0 phi=" << field[adj_mapping[0]] << "\n"
            << "dof 1 phi=" << field[adj_mapping[1]];
        }//if not bndry

        adj_cell_mapping_f.push_back(adj_mapping);
      }//for faces

      double phi     = 0.5*field[cur_cell_mapping[0]] +
                       0.5*field[cur_cell_mapping[1]];

//      phi=0.0;

      double q = source->source_value_g[0];

      //==================================== Computing interior residual
      chi_log.Log(LOG_0VERBOSE_1) << "Computing cell interior residual";
      for (int i=0; i<cell_fe_view->dofs; i++)
      {
        cell_interior_residual[lc] += q*
                                      cell_fe_view->IntV_shapeI[i];
        cell_interior_residual[lc] -= xs->sigma_rg[0]*
                                      cell_fe_view->IntV_shapeI[i]*phi;
        chi_log.Log(LOG_0VERBOSE_1) << "Vol_i[" << i << "]=" << cell_fe_view->IntV_shapeI[i];
      }
      chi_log.Log(LOG_0VERBOSE_1) << "Sigma a=" << xs->sigma_rg[0];

      //==================================== Computing surface residuals
      cell_averages[lc] = phi;
      chi_log.Log(LOG_0VERBOSE_1) << "Computing cell surface residual";
      for (int f=0; f<num_faces; f++)
      {
        int num_face_verts = 1;
        for (int fi=0; fi<num_face_verts; fi++)
        {
          double phi_adj = 0.0;
          if (slab_cell->edges[f] >= 0)
            phi_adj = 0.5*field[adj_cell_mapping_f[f][0]] +
                      0.5*field[adj_cell_mapping_f[f][1]];
          if (slab_cell->edges[f] == -1)
            phi_adj = 1.0;
//          else
//            phi_adj = 0.0;

//          phi_adj=0.0;

          chi_log.Log(LOG_0VERBOSE_1)
          << "Face " << f << "\n"
          << "phi=" << phi << " adj_phi="<<phi_adj
          << " neigbor=" <<slab_cell->edges[f];

          if (f == 0)
            cell_surface_residualL[lc] =-0.25*(phi - phi_adj);
          else
            cell_surface_residualR[lc] =-0.25*(phi - phi_adj);

        }//for face verts
      }//for face

      cell_residuals[lc] =
        std::fabs(cell_interior_residual[lc]) +
        std::fabs(cell_surface_residualL[lc]) +
        std::fabs(cell_surface_residualR[lc]);

      chi_log.Log(LOG_0VERBOSE_1)
      << "Cell residuals "
      << cell_interior_residual[lc] << " "
      << cell_surface_residualL[lc] << " "
      << cell_surface_residualR[lc] << " "
      << cell_residuals[lc];

      total_Cresidual += std::fabs(cell_interior_residual[lc]);
      total_Lresidual += std::fabs(cell_surface_residualL[lc]);
      total_Rresidual += std::fabs(cell_surface_residualR[lc]);
    }//slab
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "chi_montecarlon: Unsupported cell type encountered in "
        << "call to ResidualSource::Initialize.";
      exit(EXIT_FAILURE);
    }
  }//for local cells

  chi_log.Log(LOG_0) << "Computing total residual";
  double total_resid = 0.0;
  for (int lc=0; lc<num_local_cells; lc++)
    total_resid += cell_residuals[lc];

  total_residual = total_resid;

  chi_log.Log(LOG_0) << "Total residual = " << total_residual;
  chi_log.Log(LOG_0) << "Total cint residual = " << total_Cresidual;
  chi_log.Log(LOG_0) << "Total left residual = " << total_Lresidual;
  chi_log.Log(LOG_0) << "Total rite residual = " << total_Rresidual;


  cell_residual_cdf.resize(num_local_cells,0.0);
  double intgl = 0.0;
  for (int lc=0; lc<num_local_cells; lc++)
  {
    intgl += cell_residuals[lc];
    cell_residual_cdf[lc] = intgl/total_resid;
  }

  residual_sampler = new chi_math::CDFSampler(cell_residual_cdf);

  for (int lc=0; lc<num_local_cells; lc++)
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "Cell total residual "
      << std::setw(3) << lc
      << std::setw(12) << cell_residuals[lc]
      << " phi_avg=" << cell_averages[lc];
  }
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource::
CreateParticle(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;
  if (sample_uniformly)
    new_particle = UniformSampling(rng);
  else
    new_particle = DirectSampling(rng);

  return new_particle;
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource::
  DirectSampling(chi_montecarlon::RandomNumberGenerator* rng)
{
  int num_local_cells = grid->local_cell_glob_indices.size();
  int lc = 0;
  double sampling_normalization;
  if (sample_uniformly) {
    lc = std::floor( rng->Rand()*(num_local_cells) );
    sampling_normalization = (cell_residuals[lc]/total_residual)*
                             num_local_cells*total_residual;
  }
  else {
    lc = residual_sampler->Sample(rng->Rand());
    sampling_normalization = total_residual;
  }

  int cell_glob_index = grid->local_cell_glob_indices[lc];
  auto cell = grid->cells[cell_glob_index];

  chi_montecarlon::Particle new_particle;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_cell = (chi_mesh::CellSlab*)cell;

    int v0i = slab_cell->v_indices[0];
    int v1i = slab_cell->v_indices[1];

    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    double rn = rng->Rand();

    double center_res = std::fabs(cell_interior_residual[lc]);
    double surfL_res  = std::fabs(cell_surface_residualL[lc]);
    double surfR_res  = std::fabs(cell_surface_residualR[lc]);

    double cell_R = center_res + surfL_res + surfR_res;

    if (cell_R < 1.0e-16) {
      new_particle.alive = false;
      return new_particle;
    }


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTER
    if (rn < (center_res/cell_R))
    {
      //====================================== Sample direction
      double costheta = 2.0*rng->Rand() - 1.0;
      double theta    = acos(costheta);
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Sample position
      double w = rng->Rand();
      new_particle.pos = v0*w + v1*(1.0-w);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
        ((cell_interior_residual[lc]<0.0) ? -1.0:1.0);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_C++;

      return new_particle;
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT
    else if (rn < ((center_res + surfL_res)/cell_R))
//    if (rn < (0.5))
    {
      //====================================== Sample direction
      double costheta = rng->Rand();
      double theta    = acos(sqrt(costheta));
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Set position
      new_particle.pos = v0+chi_mesh::Vector(0,0,1.0e-8);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
        ((cell_surface_residualL[lc]<0.0) ? -1.0:1.0);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_L++;
      weights_L += new_particle.w;

      return new_particle;
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT
    else
    {
      //====================================== Sample direction
      double costheta = rng->Rand();
      double theta    = acos(-sqrt(costheta));
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Set position
      new_particle.pos = v1-chi_mesh::Vector(0,0,1.0e-8);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
        ((cell_surface_residualR[lc]<0.0) ? -1.0:1.0);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_R++;
      weights_R += new_particle.w;

      return new_particle;
    }


  }

  return new_particle;
}


//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource::
UniformSampling(chi_montecarlon::RandomNumberGenerator* rng)
{
  int num_local_cells = grid->local_cell_glob_indices.size();
  int lc = 0;
  double sampling_normalization;

  if (true) {
    lc = std::floor( rng->Rand()*(num_local_cells) );
    sampling_normalization = num_local_cells*total_residual;
  }
  else {
    lc = residual_sampler->Sample(rng->Rand());
    sampling_normalization = total_residual;
  }

  int cell_glob_index = grid->local_cell_glob_indices[lc];
  auto cell = grid->cells[cell_glob_index];

  chi_montecarlon::Particle new_particle;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_cell = (chi_mesh::CellSlab*)cell;

    int v0i = slab_cell->v_indices[0];
    int v1i = slab_cell->v_indices[1];

    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    double rn = rng->Rand();

    double center_res = std::fabs(cell_interior_residual[lc]);
    double surfL_res  = std::fabs(cell_surface_residualL[lc]);
    double surfR_res  = std::fabs(cell_surface_residualR[lc]);

    double cell_R = center_res + surfL_res + surfR_res;



    if (cell_R < 1.0e-16) {
      new_particle.alive = false;
      return new_particle;
    }


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTER
    if (rn < 0.3333333)
    {
      //====================================== Sample direction
      double costheta = 2.0*rng->Rand() - 1.0;
      double theta    = acos(costheta);
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Sample position
      double w = rng->Rand();
      new_particle.pos = v0*w + v1*(1.0-w);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
                       ((cell_interior_residual[lc]<0.0) ? -1.0:1.0)*
                       (cell_R/total_residual)*(3.0*center_res/cell_R);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_C++;

      return new_particle;
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT
    else if (rn < 0.666666666)
//    if (rn < (0.5))
    {
      //====================================== Sample direction
      double costheta = rng->Rand();
      double theta    = acos(sqrt(costheta));
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Set position
      new_particle.pos = v0+chi_mesh::Vector(0,0,1.0e-8);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
                       ((cell_surface_residualL[lc]<0.0) ? -1.0:1.0)*
                       (cell_R/total_residual)*(3.0*surfL_res/cell_R);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_L++;
      weights_L += new_particle.w;

      return new_particle;
    }
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT
    else
    {
      //====================================== Sample direction
      double costheta = rng->Rand();
      double theta    = acos(-sqrt(costheta));
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Set position
      new_particle.pos = v1-chi_mesh::Vector(0,0,1.0e-8);

      //====================================== Set quantities
      new_particle.dir = ref_dir;

      new_particle.egrp = 0;
      new_particle.w = sampling_normalization*
                       ((cell_surface_residualR[lc]<0.0) ? -1.0:1.0)*
                       (cell_R/total_residual)*(3.0*surfR_res/cell_R);
//      new_particle.w = std::fabs(sampling_normalization);
      new_particle.cur_cell_ind = cell_glob_index;
//      new_particle.alive = false;
      particles_R++;
      weights_R += new_particle.w;

      return new_particle;
    }


  }

  return new_particle;
}
