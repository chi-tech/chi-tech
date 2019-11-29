#include "chi_mesh.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LineMesh/chi_linemesh.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h"
#include "ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h"

#include "chi_log.h"
extern ChiLog chi_log;

//###################################################################
/**Creates a 1D slab mesh from a set of vertices.*/
void chi_mesh::Create1DSlabMesh(std::vector<chi_mesh::Vertex> vertices)
{
  //======================================== Checks if vertices are empty
  if (vertices.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_mesh::Create1DSlabMesh. Empty vertex list.";
    exit(EXIT_FAILURE);
  }

  //======================================== Get current mesh handler
  auto handler = chi_mesh::GetCurrentHandler();

  //======================================== Reorient 1D verts along z
  std::vector<double> distances;
  distances.reserve(vertices.size()-1);

  for (size_t v=0; v<(vertices.size()-1); ++v)
    distances.push_back((vertices[v+1]-vertices[v]).Norm());

  std::vector<Vertex> zverts;
  zverts.reserve(vertices.size());
  zverts.emplace_back(0.0,0.0,0.0);
  double last_distance = 0.0;
  for (auto& d : distances)
  {
    zverts.emplace_back(0.0,0.0,last_distance+d);
    last_distance += d;
  }

  //======================================== Create line mesh
  auto line_mesh = new chi_mesh::LineMesh;
  line_mesh->vertices = std::move(zverts);
  handler->linemesh_stack.push_back(line_mesh);

  //======================================== Add it to boundary
  auto bndry = new chi_mesh::Boundary;
  bndry->initial_mesh_continuum.line_mesh = line_mesh;

  //======================================== Add boundary to region
  auto region = new chi_mesh::Region;
  region->boundaries.push_back(bndry);

  handler->region_stack.push_back(region);

  //======================================== Create meshers
  handler->surface_mesher = new chi_mesh::SurfaceMesherPredefined;
  handler->volume_mesher = new chi_mesh::VolumeMesherLinemesh1D;

  handler->surface_mesher->Execute();

}