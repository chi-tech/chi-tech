#include "chi_mesh.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LineMesh/chi_linemesh.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h"
#include "ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h"
#include "ChiMesh/SurfaceMesh/chi_surfacemesh.h"
#include "ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a 1D slab mesh from a set of vertices.*/
void chi_mesh::Create1DSlabMesh(std::vector<double>& vertices)
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
    distances.push_back(vertices[v+1]-vertices[v]);

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
  bndry->initial_mesh_continuum.surface_mesh = surface_mesh;

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
  bndry->initial_mesh_continuum.surface_mesh = surface_mesh;

  //======================================== Add boundary to region
  auto region = new chi_mesh::Region;
  region->boundaries.push_back(bndry);

  handler->region_stack.push_back(region);

  //======================================== Create meshers
  handler->surface_mesher = new chi_mesh::SurfaceMesherPredefined;
  auto extruder = new chi_mesh::VolumeMesherExtruder;
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
    auto new_layer = new MeshLayer;
    new_layer->height = distance;
    new_layer->sub_divisions = 1;
    extruder->input_layers.push_back(new_layer);
  }

  handler->surface_mesher->Execute();
}

//###################################################################
/**Tensor product of two vectors.
 * \f$ \vec{\vec{T}} = \vec{x} \otimes \vec{y} \f$*/
chi_mesh::TensorRank2Dim3 chi_mesh::Vector3::OTimes(const Vector3& that) const
{
  TensorRank2Dim3 new_t;
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      new_t[i](j) = this->operator[](i)*that[j];

  return new_t;
}

//###################################################################
/**Dot product of vector and a rank-2 tensor.
 * \f$ \vec{w} = \vec{x} \bullet \vec{\vec{T}} \f$*/
chi_mesh::Vector3 chi_mesh::Vector3::Dot(const chi_mesh::TensorRank2Dim3& that) const
{
  chi_mesh::Vector3 new_vec;
  for (int i=0; i<3; ++i)
    new_vec(i) = this->Dot(that.t[i]);

  return new_vec;
}

//###################################################################
/**Returns a 3D vector multiplied by the given scalar from the left.
 * \f$ \vec{w} = \alpha \vec{x}\f$*/
chi_mesh::Vector3 operator*(const double value,const chi_mesh::Vector3& that)
{
  chi_mesh::Vector3 newVector;
  newVector.x = that.x*value;
  newVector.y = that.y*value;
  newVector.z = that.z*value;

  return newVector;
}

//###################################################################
/**Dot product of rank-2 tensor with a vector.
 * \f$ \vec{w} = \vec{\vec{T}} \bullet \vec{x} \f$*/
chi_mesh::Vector3 chi_mesh::TensorRank2Dim3::Dot(const chi_mesh::Vector3& v) const
{
  chi_mesh::Vector3 newVector;
  for (int i=0; i<3; ++i)
    newVector(i) = t[i].Dot(v);

  return newVector;
}

//###################################################################
/**Returns the diagonal of a rank-2 dim-3 tensor as a vector3.
 * \f$ \vec{w} = \text{diag} \vec{\vec{T}} \f$*/
chi_mesh::Vector3 chi_mesh::TensorRank2Dim3::Diag() const
{
  chi_mesh::Vector3 newVector;
  for (int i=0; i<3; ++i)
    newVector(i) = t[i][i];

  return newVector;
}

//###################################################################
/**Rank-2 dim-3 tensor multiplied from the left with a scalar.
 * \f$ \vec{\vec{W}} = \alpha \vec{\vec{T}} \f$*/
chi_mesh::TensorRank2Dim3
  operator*(const double value, const chi_mesh::TensorRank2Dim3& that)
{
  chi_mesh::TensorRank2Dim3 new_t = that;
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
      new_t[i](j) *= value;
  }

  return new_t;
}