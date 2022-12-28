#include "chi_lua.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Raytrace/raytracing.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include "ChiPhysics/FieldFunction2/fieldfunction2.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_sim_tests
{

int chiSimTest93_RayTracing(lua_State* L)
{
  const std::string fname = "chiSimTest93_RayTracing";
  chi::log.Log() << "chiSimTest93_RayTracing";

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Set parameters
  const size_t num_groups = 1;
  const size_t num_moments = 1;
  const double sigma_t = 0.27;

  //============================================= Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr fem_sdm_ptr = chi_math::SpatialDiscretization_PWLD::New(grid_ptr);
  const auto& fem_sdm = *fem_sdm_ptr;

  chi_math::UnknownManager phi_uk_man;
  for (size_t m=0; m<num_moments; ++m)
    phi_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);

  const size_t num_fem_local_dofs = fem_sdm.GetNumLocalDOFs(phi_uk_man);
  const size_t num_fem_globl_dofs = fem_sdm.GetNumGlobalDOFs(phi_uk_man);

  chi::log.Log() << "Num local FEM DOFs: " << num_fem_local_dofs;
  chi::log.Log() << "Num globl FEM DOFs: " << num_fem_globl_dofs;

  //============================================= Define tallies
  std::vector<double> fem_tally(num_fem_local_dofs,0.0);

  //============================================= Define particle
  //                                              data structure
  typedef chi_mesh::Vector3 Vec3;
  struct Particle
  {
    Vec3 position = {0.0,0.0,0.0};
    Vec3 direction = {0.0,0.0,0.0};
    int  energy_group = 0;
    double weight = 1.0;

    uint64_t cell_id = 0;

    bool alive = true;
  };

  //============================================= Define source position
  //                                              and find cell containing it
  const Vec3 source_pos = {0.0,0.0,0.0};

  chi_mesh::Cell const* source_cell_ptr = nullptr;

  for (auto& cell : grid.local_cells)
    if (grid.CheckPointInsideCell(cell, source_pos))
    {
      source_cell_ptr = &cell;
      break;
    }
  if (source_cell_ptr == nullptr)
    throw std::logic_error(fname + ": Source cell not found.");

  const uint64_t source_cell_id = source_cell_ptr->global_id;

  //============================================= Define lambdas
  chi_math::RandomNumberGenerator rng;
  auto SampleRandomDirection = [&rng]()
  {
    double costheta = 2.0*rng.Rand() - 1.0;
    double theta    = acos(costheta);
    double varphi   = rng.Rand()*2.0*M_PI;

    return chi_mesh::Vector3{sin(theta) * cos(varphi),
                             sin(theta) * sin(varphi),
                             cos(theta)};
  };

  /**PWLD Tally*/
  auto ContributePWLDTally = [&fem_sdm,&grid,&fem_tally,&phi_uk_man,&sigma_t](
    const chi_mesh::Cell& cell,
    const Vec3& positionA,
    const Vec3& positionB,
    const Vec3& omega,
    const int g,
    double weight)
  {
    const auto& cell_mapping = fem_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    std::vector<double> segment_lengths;
    chi_mesh::PopulateRaySegmentLengths(/*input*/grid,
                                        /*input*/cell,
                                        /*input*/positionA,
                                        /*input*/positionB,
                                        /*input*/omega,
                                        /*output*/segment_lengths);

    std::vector<double> shape_values;
    std::vector<double> shape_values_a;
    std::vector<double> shape_values_b;

    cell_mapping.ShapeValues(positionA, shape_values_a);

    double d_run_sum = 0.0;
    for (const auto& segment_length : segment_lengths)
    {
      d_run_sum += segment_length;
      const double& d = d_run_sum;

      cell_mapping.ShapeValues(positionA+omega*d, shape_values_b);

      const auto& N_i = shape_values_a;
      const auto& N_f = shape_values_b;

      for (size_t i=0; i<num_nodes; ++i)
      {
        for (size_t m=0; m < num_moments; ++m)
        {
          const int64_t dof_map = fem_sdm.MapDOFLocal(cell,i,phi_uk_man,m,g);

          const double& ell = segment_length;

          double w_avg  = (N_i[i] / sigma_t) * (1.0 - exp(-sigma_t * ell));
          w_avg += ((N_f[i] - N_i[i]) / (sigma_t * sigma_t * ell)) *
                   (1.0 - (1+sigma_t*ell)*exp(-sigma_t*ell));
          w_avg *= weight/ell;

          fem_tally[dof_map] += ell * w_avg ;
        }//for moment m
      }//for node i

      shape_values_a = shape_values_b;
      weight *= exp(-sigma_t*segment_length);
    }//for d
  };

  /**Lambda to get cell bounding box.*/
  auto GetCellApproximateSize = [&grid](const chi_mesh::Cell& cell)
  {
    const auto& v0 = grid.vertices[cell.vertex_ids[0]];
    double xmin = v0.x, xmax = v0.x;
    double ymin = v0.y, ymax = v0.y;
    double zmin = v0.z, zmax = v0.z;

    for (uint64_t vid : cell.vertex_ids)
    {
      const auto& v = grid.vertices[vid];

      xmin = std::min(xmin, v.x); xmax = std::max(xmax, v.x);
      ymin = std::min(ymin, v.y); ymax = std::max(ymax, v.y);
      zmin = std::min(zmin, v.z); zmax = std::max(zmax, v.z);
    }

    return (chi_mesh::Vector3(xmin, ymin, zmin) -
            chi_mesh::Vector3(xmax, ymax, zmax)).Norm();
  };

  //============================================= Create raytracer
  std::vector<double> cell_sizes(grid.local_cells.size(), 0.0);
  for (const auto& cell : grid.local_cells)
    cell_sizes[cell.local_id] = GetCellApproximateSize(cell);

  chi_mesh::RayTracer ray_tracer(grid, cell_sizes);

  //============================================= Run rays
  const auto PWLD =
    chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;

  const size_t num_particles = 10'000'000;
  for (size_t n=0; n<num_particles; ++n)
  {
    if (n % size_t(num_particles/10.0) == 0)
      std::cout << "#particles = " << n << "\n";
    //====================================== Create the particle
    const auto omega = SampleRandomDirection();
    Particle particle{/*position=*/  source_pos,
                      /*direction=*/ omega,
                      /*e_group=*/   0,
                      /*weight=*/    1.0,
                      /*cell_id=*/   source_cell_id,
                      /*alive=*/     true};

    while (particle.alive)
    {
      //=============================== Get the current cell
      const auto& cell = grid.cells[particle.cell_id];

      //=============================== Perform the trace
      //                                to the next surface
      auto destination_info = ray_tracer.TraceRay(cell,
                                                  particle.position,
                                                  particle.direction);

      const Vec3& end_of_track_position = destination_info.pos_f;

      //=============================== Make tally contribution
      const int g = particle.energy_group;
      if (fem_sdm.type == PWLD)
        ContributePWLDTally(cell,
                            particle.position,     //positionA
                            end_of_track_position, //positionB
                            particle.direction,    //omega
                            g,                     //
                            particle.weight);      //weight at A

      //=============================== Process cell transfer
      //                                or death
      if (not destination_info.particle_lost)
      {
        const auto& f = destination_info.destination_face_index;
        const auto& current_cell_face = cell.faces[f];

        if (current_cell_face.has_neighbor)
          particle.cell_id = current_cell_face.neighbor_id;
        else
          particle.alive = false; //Death at the boundary
      }
      else
      {
        std::cout << "particle" << n << " lost "
                  << particle.position.PrintStr() << " "
                  << particle.direction.PrintStr() << " "
                  << "\n";
        break;
      }

      const auto& pA = particle.position;
      const auto& pB = end_of_track_position;
      particle.weight *= exp(-sigma_t*(pB-pA).Norm()); //Attenuation
      particle.position = end_of_track_position;
    }//while ray alive

  }//for ray n

  //============================================= Post process tallies
  for (const auto& cell : grid.local_cells)
  {
    //====================================== Compute mass matrix
    //                                       and its inverse
    const auto& cell_mapping = fem_sdm.GetCellMapping(cell);
    const auto& qp_data = cell_mapping.MakeVolumeQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const double normalization = num_particles*cell_mapping.CellVolume();

    MatDbl M(num_nodes, VecDbl(num_nodes, 0.0));
    VecDbl Vi(num_nodes, 0.0);
    for (auto qp : qp_data.QuadraturePointIndices())
      for (size_t i=0; i<num_nodes; ++i)
      {
        for (size_t j=0; j<num_nodes; ++j)
          M[i][j] += qp_data.ShapeValue(i,qp) * qp_data.ShapeValue(j, qp) *
                     qp_data.JxW(qp);

        Vi[i] += qp_data.ShapeValue(i,qp) * qp_data.JxW(qp);
      }

    auto M_inv = chi_math::Inverse(M);

    //====================================== Apply projection
    VecDbl b(num_nodes, 0.0);
    for (size_t m=0; m<num_moments; ++m)
      for (size_t g=0; g<num_groups; ++g)
      {
        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t imap = fem_sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          b[i] = fem_tally[imap]/num_particles;
        }

        auto c = chi_math::MatMul(M_inv, b);

        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t imap = fem_sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          fem_tally[imap] = c[i];
        }
      }//for group g

  }//for cell

  //============================================= Localize zeroth moment
  //This routine extracts a single moment vector
  //from the vector that contains multiple moments
  const chi_math::UnknownManager m0_uk_man(
    {chi_math::Unknown(chi_math::UnknownType::VECTOR_N,num_groups)});
  const size_t num_m0_dofs = fem_sdm.GetNumLocalDOFs(m0_uk_man);

  std::vector<double> m0_phi(num_m0_dofs, 0.0);

  fem_sdm.CopyVectorWithUnknownScope(fem_tally,   //from vector
                                     m0_phi,      //to vector
                                     phi_uk_man,  //from dof-structure
                                     0,           //from unknown-id
                                     m0_uk_man,   //to dof-structure
                                     0);          //to unknown-id

  //============================================= Create field function
  auto ff = std::make_shared<chi_physics::FieldFunction2>(
    "Phi",                                           //Text name
    fem_sdm_ptr,                                         //Spatial Discr.
    chi_math::Unknown(chi_math::UnknownType::VECTOR_N,num_groups) //Unknown
  );

  ff->UpdateFieldVector(m0_phi);
  ff->ExportToVTK(fname);

  return 0;
}

}//namespace chi_unit_sim_tests