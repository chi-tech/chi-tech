#ifndef CHITECH_CHIMESH_RAYTRACER_H
#define CHITECH_CHIMESH_RAYTRACER_H

#include <utility>

#include "raytracing_old.h"

namespace chi_mesh
{
  struct RayTracerOutputInformation
  {
    double       distance_to_surface = 0.0;
    Vector3      pos_f;
    unsigned int destination_face_index = 0;
    uint64_t     destination_face_neighbor = 0;
    bool         particle_lost = false;
    std::string  lost_particle_info;
  };

//###################################################################
/**Raytracer object.*/
  class RayTracer
  {
  private:
    const chi_mesh::MeshContinuum& reference_grid;
    std::vector<double> cell_sizes;
    double epsilon_nudge      = 1.0e-8;
    double backward_tolerance = 1.0e-10;
    double extension_distance = 1.0e5;
  public:
    bool   perform_concavity_checks = true;

    explicit
    RayTracer(const chi_mesh::MeshContinuum& grid,
              std::vector<double> in_cell_sizes,
              bool   in_perform_concavity_checks = true) :
      reference_grid(grid),
      cell_sizes(std::move(in_cell_sizes)),
      perform_concavity_checks(in_perform_concavity_checks)
    {}

  private:
    const chi_mesh::MeshContinuum& Grid() const;

    void SetTolerancesFromCellSize(double cell_size)
    {
      epsilon_nudge = cell_size * 1.0e-2;
      backward_tolerance = cell_size * 1.0e-10;
      extension_distance = 3.0 * cell_size;
    }

  public:
    /**Traces a ray with an initial position either within the cell or
     * on the cell surface, and with a direction vector pointing inward
     * toward the cell. If the ray-trace fails the particle will be marked
     * as lost.*/
    RayTracerOutputInformation
    TraceRay(const Cell& cell,
             Vector3& pos_i,
             Vector3& omega_i,
             int function_depth=0);

    /**Traces a ray with an initial position, presumed to be outside the cell,
     * to an incident face.*/
    RayTracerOutputInformation
    TraceIncidentRay(const Cell& cell,
                     const Vector3& pos_i,
                     const Vector3& omega_i);

  private:
    void TraceSlab(const Cell& cell,
                   Vector3& pos_i,
                   Vector3& omega_i,
                   bool& intersection_found,
                   bool& backward_tolerance_hit,
                   RayTracerOutputInformation& oi);
    void TracePolygon(const Cell& cell,
                      Vector3& pos_i,
                      Vector3& omega_i,
                      bool& intersection_found,
                      bool& backward_tolerance_hit,
                      RayTracerOutputInformation& oi);
    void TracePolyhedron(const Cell& cell,
                         Vector3& pos_i,
                         Vector3& omega_i,
                         bool& intersection_found,
                         bool& backward_tolerance_hit,
                         RayTracerOutputInformation& oi);
  };
}

#endif //CHITECH_CHIMESH_RAYTRACER_H
