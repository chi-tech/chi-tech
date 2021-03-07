#include "ChiMesh/chi_mesh.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LineMesh/chi_linemesh.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h"
#include "ChiMesh/SurfaceMesh/chi_surfacemesh.h"
#include "ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a 2D orthogonal mesh from a set of vertices in x and y.
 * The vertices along a dimension merely represents the divisions. They
 * are not the complete vertices defining a cell. For example:
\code
std::vector<chi_mesh::Vertex> vertices_x = {0.0,1.0,2.0};
std::vector<chi_mesh::Vertex> vertices_y = {0.0,1.0,2.0};
chi_mesh::Create2DOrthoMesh(vertices_x,vertices_y);
\endcode

This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.

 */
void chi_mesh::Create2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                 std::vector<double>& vertices_1d_y)
{
  auto surface_mesh =
    chi_mesh::SurfaceMesh::CreateFromDivisions(vertices_1d_x,vertices_1d_y);

  //======================================== Get current mesh handler
  auto handler = chi_mesh::GetCurrentHandler();
  handler->surface_mesh_stack.push_back(surface_mesh);

  //======================================== Add it to boundary
  auto bndry = new chi_mesh::Boundary;
  bndry->initial_mesh_continuum->surface_mesh = surface_mesh;

  //======================================== Add boundary to region
  auto region = new chi_mesh::Region;
  region->boundaries.push_back(bndry);

  handler->region_stack.push_back(region);

  //======================================== Create meshers
  handler->surface_mesher = new chi_mesh::SurfaceMesherPredefined;
  handler->volume_mesher = new chi_mesh::VolumeMesherPredefined2D;

  handler->surface_mesher->Execute();
}

//###################################################################
/**Creates a 3D orthogonal mesh from a set of vertices in x, y and z.
 * The vertices along a dimension merely represents the divisions. They
 * are not the complete vertices defining a cell. For example:
\code
std::vector<double> vertices_x = {0.0,1.0,2.0};
std::vector<double> vertices_y = {0.0,1.0,2.0};
std::vector<double> vertices_z = {0.0,1.0,2.0};
chi_mesh::Create3DOrthoMesh(vertices_x,vertices_y,vertices_z);
\endcode

This code will create a 2x2x2 mesh with \f$ \vec{x} \in [0,2]^3 \f$.

 */
void chi_mesh::Create3DOrthoMesh(std::vector<double>& vertices_x,
                                 std::vector<double>& vertices_y,
                                 std::vector<double>& vertices_z)
{
  auto surface_mesh =
    chi_mesh::SurfaceMesh::CreateFromDivisions(vertices_x,vertices_y);

  //======================================== Get current mesh handler
  auto handler = chi_mesh::GetCurrentHandler();
  handler->surface_mesh_stack.push_back(surface_mesh);

  //======================================== Add it to boundary
  auto bndry = new chi_mesh::Boundary;
  bndry->initial_mesh_continuum->surface_mesh = surface_mesh;

  //======================================== Add boundary to region
  auto region = new chi_mesh::Region;
  region->boundaries.push_back(bndry);

  handler->region_stack.push_back(region);

  //======================================== Create meshers
  handler->surface_mesher = new chi_mesh::SurfaceMesherPredefined;
  auto extruder = new chi_mesh::VolumeMesherExtruder(surface_mesh);
  handler->volume_mesher = extruder;

  //======================================== Populate layers
  std::vector<double> distances;
  distances.reserve(vertices_z.size()-1);

  for (size_t v=0; v<(vertices_z.size()-1); ++v)
    distances.push_back(vertices_z[v+1]-vertices_z[v]);

  std::vector<Vertex> zverts;
  zverts.reserve(vertices_z.size());
  zverts.emplace_back(0.0,0.0,0.0);
  double last_distance = 0.0;
  for (auto& d : distances)
  {
    zverts.emplace_back(0.0,0.0,last_distance+d);
    last_distance += d;
  }

  for (auto distance : distances)
  {
    chi_mesh::VolumeMesherExtruder::MeshLayer new_layer;
    new_layer.height = distance;
    new_layer.sub_divisions = 1;
    extruder->input_layers.push_back(new_layer);
  }

  handler->surface_mesher->Execute();
}