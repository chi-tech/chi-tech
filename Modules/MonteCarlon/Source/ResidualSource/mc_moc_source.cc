#include "mc_moc_source.h"

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
#include <ChiMath/Quadratures/quadrature_gausslegendre.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualMOCSource::ResidualMOCSource(
  chi_physics::FieldFunction *in_resid_ff, bool use_uniform_sampling) :
  sample_uniformly(use_uniform_sampling)
{
  resid_ff = in_resid_ff;
  particles_L = 0;
  particles_R = 0;

  weights_L = 0.0;
  weights_R = 0.0;
}

//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void chi_montecarlon::ResidualMOCSource::
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
  cell_dof_phi.resize(num_local_cells);
  cell_sigma_s.resize(num_local_cells);
  cell_sigma_t.resize(num_local_cells);

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

    cell_sigma_s[lc] = xs->sigma_tg[0] - xs->sigma_rg[0];
    cell_sigma_t[lc] = xs->sigma_tg[0];

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      chi_log.Log(LOG_0VERBOSE_1) << "**************** Cell " << cell_glob_index;
      auto slab_cell = (chi_mesh::CellSlabV2*)cell;
      auto cell_fe_view =
        static_cast<SlabFEView*>(
          resid_sdm_pwl->MapFeView(cell_glob_index));

      //==================================== Creating current cell dof-mapping
      std::vector<int> dofs_to_map(cell_fe_view->dofs);
      std::vector<int> cells_to_map;
      std::vector<int> cur_cell_mapping;
      cells_to_map.resize(cell_fe_view->dofs);
      for (int i=0; i<cell_fe_view->dofs; i++)
      {
        dofs_to_map[i] = i;
        cells_to_map[i] = cell_glob_index;
      }

      ff_interp.CreatePWLDMapping(
        resid_ff,dofs_to_map,cells_to_map,&cur_cell_mapping);

      chi_log.Log(LOG_0VERBOSE_1)
        << "dof 0 phi=" << field[cur_cell_mapping[0]] << "\n"
        << "dof 1 phi=" << field[cur_cell_mapping[1]];

      //==================================== Creating adj cell dof-mapping
      std::vector<std::vector<int>> adj_cell_mapping_f;

      int num_faces = 2;
      for (int f=0; f<num_faces; f++)
      {
        int adj_cell_index = slab_cell->faces[f].neighbor;

        std::vector<int> adj_mapping;

        if (adj_cell_index >= 0)
        {
          int adj_num_dofs = 2;
          for (int i=0; i<adj_num_dofs; i++)
          {
            dofs_to_map[i] = i;
            cells_to_map[i] = adj_cell_index;
          }

          ff_interp.CreatePWLDMapping(
            resid_ff,dofs_to_map,cells_to_map,&adj_mapping);

          chi_log.Log(LOG_0VERBOSE_1)
            << "adj_cell " << adj_cell_index << "\n"
            << "dof 0 phi=" << field[adj_mapping[0]] << "\n"
            << "dof 1 phi=" << field[adj_mapping[1]];
        }//if not bndry

        adj_cell_mapping_f.push_back(adj_mapping);
      }//for faces

      double phi     = 0.5*field[cur_cell_mapping[0]] +
                       0.5*field[cur_cell_mapping[1]];

      double q = source->source_value_g[0];

      //==================================== Computing average dof flux
      cell_averages[lc] = phi;
      cell_dof_phi[lc].resize(2,0.0);
      for (int f=0; f<num_faces; f++)
      {
        int num_face_verts = 1;
        for (int fi=0; fi<num_face_verts; fi++)
        {
          if (slab_cell->faces[f].neighbor >= 0)
          {
            cell_dof_phi[lc][f] =
              0.5*field[cur_cell_mapping[f]] +
              0.5*field[adj_cell_mapping_f[f][abs(f-1)]];
//            cell_dof_phi[lc][f] = 0.0; //TODO: remove if RMC
          }
          else
          {
            double phi_adj = 0.0;
//            if (f==0)
//              phi_adj = 1.0;
//            else
              phi_adj = field[cur_cell_mapping[f]];

            cell_dof_phi[lc][f] = 0.5*field[cur_cell_mapping[f]] +
                                  0.5*phi_adj;

//            cell_dof_phi[lc][f] = 0.0; //TODO: remove before flight


          }
        }//for face verts
      }//for face
    }//slab
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "chi_montecarlon: Unsupported cell type encountered in "
        << "call to ResidualSource::Initialize.";
      exit(EXIT_FAILURE);
    }
  }//for local cells


  for (int lc=0; lc<num_local_cells; lc++)
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "Cell nodal values "
      << std::setw(3) << lc
      << std::setw(12) << cell_dof_phi[lc][0]
      << std::setw(12) << cell_dof_phi[lc][1]
      << " phi_avg=" << cell_dof_phi[lc][0]*0.5+cell_dof_phi[lc][1]*0.5;
  }

  //================================================== Raytrace subdivisions
  num_subdivs = 80;
  size_t num_angles = 32;
  quadrature.Initialize(num_angles);

  cell_phi_star.resize(
    num_local_cells,
    std::vector<double>(num_subdivs+1,0.0));
  cell_z_i_star.resize(
    num_local_cells,
    std::vector<double>(num_subdivs+1,0.0));

  for (size_t n=0; n<quadrature.abscissae.size(); n++)
  {
    chi_log.Log(LOG_0VERBOSE_0)
    << "Angle " << n << " " << quadrature.abscissae[n];
    double mu = quadrature.abscissae[n];
    double q_weight = quadrature.weights[n];
    if (mu>0.0)
    {
      double psi_z_i = -0.5*(cell_dof_phi[0][0]-1.0);


      for (int lc=0; lc<num_local_cells; lc++)
      {
        auto cell = grid->cells[lc];

        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
        if (cell->Type() == chi_mesh::CellType::SLAB)
        {
          auto slab_cell = (chi_mesh::CellSlabV2*)cell;

          int v0i = slab_cell->vertex_ids[0];
          int v1i = slab_cell->vertex_ids[1];

          chi_mesh::Vertex& v0 = *grid->nodes[v0i];
          chi_mesh::Vertex& v1 = *grid->nodes[v1i];

          double dz = (v1-v0).Norm();
          double dzstar = dz/num_subdivs;
          double z_i = v0.z;
          double s_t = cell_sigma_t[lc]/mu;

          for (int i=0; i<(num_subdivs+1); i++)
          {
            double zistar = z_i + dzstar*i;

            if (i>0)
            {
              psi_z_i = psi_z_i*exp(-s_t*dzstar);

              double w = (dzstar*i-0.5*dzstar)/dz;
              double phi_zstar =
                  (1.0-w)*cell_dof_phi[lc][0] +
                      (w-0.0)*cell_dof_phi[lc][1];
              double q = 0.5*(cell_sigma_s[lc] - cell_sigma_t[lc])*phi_zstar;

              double nabla_phi = 0.5*(cell_dof_phi[lc][1]-cell_dof_phi[lc][0])/dz;

              q -= mu*nabla_phi;

              psi_z_i += q*(1.0-exp(-s_t*dzstar))/cell_sigma_t[lc];

              if (n==17)
              {
                chi_log.Log(LOG_0VERBOSE_1)
                  << zistar << " " << q;
              }
            }

            cell_phi_star[lc][i] += q_weight*psi_z_i;
            cell_z_i_star[lc][i]  = zistar;
          }
        }//slab
      }//for local cell
    }//if mu>0.0
    else
    {
      double psi_z_i = -0.5*(cell_dof_phi[num_local_cells-1][1]-0.0);

      for (int lc=(num_local_cells-1); lc>=0; lc--)
      {
        auto cell = grid->cells[lc];

        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
        if (cell->Type() == chi_mesh::CellType::SLAB)
        {
          auto slab_cell = (chi_mesh::CellSlabV2*)cell;

          int v0i = slab_cell->vertex_ids[0];
          int v1i = slab_cell->vertex_ids[1];

          chi_mesh::Vertex& v0 = *grid->nodes[v0i];
          chi_mesh::Vertex& v1 = *grid->nodes[v1i];

          double dz = (v1-v0).Norm();
          double dzstar = dz/num_subdivs;
          double z_i = v0.z;
          double s_t = cell_sigma_t[lc]/std::fabs(mu);

          for (int i=0; i<(num_subdivs+1); i++)
          {
            double zistar = z_i + dzstar*i;

            if (i>0)
            {
              psi_z_i = psi_z_i*exp(-s_t*dzstar);

              double w = (dzstar*i-0.5*dzstar)/dz;
              double phi_zstar =
                  (w-0.0)*cell_dof_phi[lc][0] +
                      (1.0-w)*cell_dof_phi[lc][1];
              double q = 0.5*(cell_sigma_s[lc] - cell_sigma_t[lc])*phi_zstar;

              double nabla_phi = 0.5*(cell_dof_phi[lc][1]-cell_dof_phi[lc][0])/dz;

              q -= mu*nabla_phi;

              psi_z_i += q*(1.0-exp(-s_t*dzstar))/cell_sigma_t[lc];
            }


            int irev = num_subdivs - i; //reverse i
            cell_phi_star[lc][irev] += q_weight*psi_z_i;
            cell_z_i_star[lc][irev]  = zistar;
          }
        }//slab
      }//for local cell
    }//if mu<0.0
  }//for mu

  //============================================= Constructing cell totals
  total_abs_source = 0.0;
  cell_abs_total_source.resize(num_local_cells, 0.0);
  cell_total_source.resize(num_local_cells, 0.0);
  cell_subintvl_source.resize(num_local_cells);
  double dz = 5.0/num_local_cells;
  for (int lc=0; lc<num_local_cells; lc++)
  {
    cell_subintvl_source[lc].resize(num_subdivs, 0.0);
    for (int i=0; i<(num_subdivs+1); i++)
    {
//      chi_log.Log(LOG_0VERBOSE_1)
//        << cell_z_i_star[lc][i] << " "
//        << cell_phi_star[lc][i];

      if (i<num_subdivs)
      {
        double phi_avg =
          0.5*cell_phi_star[lc][i] +
          0.5*cell_phi_star[lc][i+1];

        cell_abs_total_source[lc] += std::fabs(phi_avg)*cell_sigma_s[lc]*dz/num_subdivs;
        cell_total_source[lc] += phi_avg;
        cell_subintvl_source[lc][i] = phi_avg*cell_sigma_s[lc]*dz/num_subdivs;
      }
    }//for subdivs
    total_abs_source += cell_abs_total_source[lc];
  }//for local cell

  chi_log.Log(LOG_0VERBOSE_1)
      << "\nTotal absolute source = " << total_abs_source;

  //============================================= Constructing cdf
  cell_cdf.resize(num_local_cells);
  double cumulator = 0.0;
  for (int lc=0; lc<(num_local_cells); lc++)
  {
    cumulator +=cell_abs_total_source[lc];
    cell_cdf[lc] = cumulator/total_abs_source;
  }

  cell_sampler = new chi_math::CDFSampler(cell_cdf);

  //============================================= Print cell averages
  for (int lc=0; lc<num_local_cells; lc++)
  {
    chi_log.Log(LOG_0VERBOSE_1)
        << "Cell avg uncollided flux "
        << lc << " "
        << cell_total_source[lc]/num_subdivs;
  }
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualMOCSource::
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
chi_montecarlon::Particle chi_montecarlon::ResidualMOCSource::
DirectSampling(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  int num_local_cells = grid->local_cell_glob_indices.size();
  int lc = 0;
//  lc = std::floor( rng->Rand()*(num_local_cells) );
  lc = cell_sampler->Sample(rng->Rand());

  int cell_glob_index = grid->local_cell_glob_indices[lc];
  auto cell = grid->cells[cell_glob_index];


  //====================================== Sample direction
  double costheta = 2.0*rng->Rand() - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_cell = (chi_mesh::CellSlabV2*)cell;

    int v0i = slab_cell->vertex_ids[0];
    int v1i = slab_cell->vertex_ids[1];

    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    //====================================== Sample position
    double dz = (v1-v0).Norm();
    double dzstar = (v1-v0).Norm()/num_subdivs;
    double s = rng->Rand();
    double w = rng->Rand();

    double sampling_normalization = total_abs_source;
    double z=0.0;
    double cumulator = 0.0;
    for (int i=0; i<num_subdivs; i++)
    {
      cumulator += std::fabs(cell_subintvl_source[lc][i]);
      if (s < (cumulator/std::fabs(cell_abs_total_source[lc])))
      {
        z = cell_z_i_star[lc][i] + w*dzstar;
        sampling_normalization *=
          ((cell_subintvl_source[lc][i]<0.0) ? -1.0 : 1.0);
        break;
      }
    }

    new_particle.pos = chi_mesh::Vector(0.0,0.0,z);
    new_particle.egrp = 0;
    new_particle.w = sampling_normalization;

    new_particle.cur_cell_ind = cell_glob_index;

    if (w<0.0) new_particle.alive = false;

  }

  return new_particle;
}


//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualMOCSource::
UniformSampling(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  int num_local_cells = grid->local_cell_glob_indices.size();
  int lc = 0;
  lc = std::floor( rng->Rand()*(num_local_cells) );
//  lc = cell_sampler->Sample(rng->Rand());

  int cell_glob_index = grid->local_cell_glob_indices[lc];
  auto cell = grid->cells[cell_glob_index];


  //====================================== Sample direction
  double costheta = 2.0*rng->Rand() - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
  if (cell->Type() == chi_mesh::CellType::SLAB)
  {
    auto slab_cell = (chi_mesh::CellSlabV2*)cell;

    int v0i = slab_cell->vertex_ids[0];
    int v1i = slab_cell->vertex_ids[1];

    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    //====================================== Sample position
    double dz = (v1-v0).Norm();
    double dzstar = (v1-v0).Norm()/num_subdivs;
    double w = rng->Rand();

    double sampling_normalization = num_local_cells;
    double z=0.0;
    int i = std::floor(rng->Rand()*num_subdivs);
    z = cell_z_i_star[lc][i] + w*dzstar;
    sampling_normalization *= cell_subintvl_source[lc][i]*num_subdivs;

    new_particle.pos = chi_mesh::Vector(0.0,0.0,z);
    new_particle.egrp = 0;
    new_particle.w = sampling_normalization;

    new_particle.cur_cell_ind = cell_glob_index;
  }

  return new_particle;
}
