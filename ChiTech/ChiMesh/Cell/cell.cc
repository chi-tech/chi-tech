#include "cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiDataTypes/byte_array.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Provides the text name associated with a cell type.*/
std::string chi_mesh::CellTypeName(const CellType type)
{
  switch (type)
  {
    case CellType::GHOST:               return "GHOST";
    case CellType::SLAB:                return "SLAB";
//    case CellType::SPHERICAL_SHELL:     return "SPHERICAL_SHELL";
//    case CellType::CYLINDRICAL_ANNULUS: return "CYLINDRICAL_ANNULUS";
    case CellType::TRIANGLE:            return "TRIANGLE";
    case CellType::QUADRILATERAL:       return "QUADRILATERAL";
    case CellType::POLYGON:             return "POLYGON";
    case CellType::TETRAHEDRON:         return "TETRAHEDRON";
    case CellType::HEXAHEDRON:          return "HEXAHEDRON";
    case CellType::POLYHEDRON:          return "POLYHEDRON";
    default: return "NONE";
  }
}

//###################################################################
/**Copy constructor*/
chi_mesh::Cell::Cell(const Cell &other) :
  cell_type(other.cell_type),
  cell_sub_type(other.cell_sub_type),
  global_id(other.global_id),
  local_id(other.local_id),
  partition_id(other.partition_id),
  centroid(other.centroid),
  material_id(other.material_id),
  vertex_ids(other.vertex_ids),
  faces(other.faces)
{
}

//###################################################################
/**Move constructor*/
chi_mesh::Cell::Cell(Cell &&other) noexcept :
  cell_type(other.cell_type),
  cell_sub_type(other.cell_sub_type),
  global_id(other.global_id),
  local_id(other.local_id),
  partition_id(other.partition_id),
  centroid(other.centroid),
  material_id(other.material_id),
  vertex_ids(std::move(other.vertex_ids)),
  faces(std::move(other.faces))
{
}


//###################################################################
/**Copy operator.*/
chi_mesh::Cell& chi_mesh::Cell::operator=(const Cell &other)
{
  global_id    = other.global_id;
  local_id     = other.local_id;
  partition_id = other.partition_id;
  centroid     = other.centroid;
  material_id  = other.material_id;
  vertex_ids   = other.vertex_ids;
  faces        = other.faces;

  return *this;
}


//###################################################################
/**Determines the neighbor's partition and whether its local or not.*/
bool chi_mesh::CellFace::
  IsNeighborLocal(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor) return false;
  if (chi_mpi.process_count == 1) return true;

  auto& adj_cell = grid.cells[neighbor_id];

  return (adj_cell.partition_id == static_cast<uint64_t>(chi_mpi.location_id));
}

//###################################################################
/**Determines the neighbor's partition.*/
int chi_mesh::CellFace::
  GetNeighborPartitionID(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor) return -1;
  if (chi_mpi.process_count == 1) return 0;

  auto& adj_cell = grid.cells[neighbor_id];

  return adj_cell.partition_id;
}

//###################################################################
/**Determines the neighbor's local id.*/
int chi_mesh::CellFace::
  GetNeighborLocalID(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor) return -1;
  if (chi_mpi.process_count == 1) return neighbor_id; //cause global_ids=local_ids

  auto& adj_cell = grid.cells[neighbor_id];

  if (adj_cell.partition_id != chi_mpi.location_id)
    throw std::logic_error("Cell local ID requested from a non-local cell.");

  return adj_cell.local_id;
}

//###################################################################
/**Determines the neighbor's associated face.*/
int chi_mesh::CellFace::
  GetNeighborAssociatedFace(const chi_mesh::MeshContinuum& grid) const
{
  const auto& cur_face = *this; //just for readability
  //======================================== Check index validity
  if ((not cur_face.has_neighbor) || (not cur_face.IsNeighborLocal(grid)))
  {
    std::stringstream outstr;
    outstr
      << "Invalid cell index encountered in call to "
      << "CellFace::GetNeighborAssociatedFace. Index points "
      << "to either a boundary"
      << "or a non-local cell.";
    throw std::logic_error(outstr.str());
  }

  const auto& adj_cell = grid.local_cells[cur_face.GetNeighborLocalID(grid)];

  int associated_face = -1;
  std::set<uint64_t> cfvids(cur_face.vertex_ids.begin(),
                            cur_face.vertex_ids.end()); //cur_face vertex ids

  //======================================== Loop over adj cell faces
  int af=-1;
  for (const auto& adj_face : adj_cell.faces)
  {
    ++af;
    std::set<uint64_t> afvids(adj_face.vertex_ids.begin(),
                              adj_face.vertex_ids.end()); //adj_face vertex ids

    if (afvids == cfvids) {associated_face = af; break;}
  }

  //======================================== Check associated face validity
  if (associated_face<0)
  {
    std::stringstream outstr;
    outstr
      << "Could not find associated face in call to "
      << "CellFace::GetNeighborAssociatedFace.\n"
      << "Reference face with centroid at: "
      << cur_face.centroid.PrintS() << "\n"
      << "Adjacent cell: " << adj_cell.global_id << "\n";
    for (size_t afi=0; afi < adj_cell.faces.size(); afi++)
    {
      outstr
        << "Adjacent cell face " << afi << " centroid "
        << adj_cell.faces[afi].centroid.PrintS();
    }
    throw std::runtime_error(outstr.str());
  }

  return associated_face;
}

//###################################################################
/**Computes the face area.*/
double chi_mesh::CellFace::ComputeFaceArea(chi_mesh::MeshContinuum& grid) const
{
  if (vertex_ids.size() <= 1)
    return 1.0;
  else if (vertex_ids.size() == 2)
  {
    const auto& v0 = grid.vertices[vertex_ids[0]];
    const auto& v1 = grid.vertices[vertex_ids[1]];

    return (v1 - v0).Norm();
  }
  else
  {
    double area = 0.0;
    auto& v2 = centroid;
    const auto num_verts = vertex_ids.size();
    for (uint64_t v=0; v<num_verts; ++v)
    {
      uint64_t vid0 = vertex_ids[v];
      uint64_t vid1 = (v < (num_verts-1))? vertex_ids[v+1] : vertex_ids[0];

      const auto& v0 = grid.vertices[vid0];
      const auto& v1 = grid.vertices[vid1];

      auto v01 = v1-v0;
      auto v02 = v2-v0;

      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0,v01);
      J.SetColJVec(1,v02);
      J.SetColJVec(2,chi_mesh::Vector3(1.0,1.0,1.0));

      area += 0.5*std::fabs(J.Det());
    }

    return area;
  }

}

//###################################################################
/**Serializes a face into a vector of bytes.*/
chi_data_types::ByteArray chi_mesh::CellFace::Serialize() const
{
  chi_data_types::ByteArray raw;

  raw.Write<size_t>(vertex_ids.size());
  for (uint64_t vid : vertex_ids)
    raw.Write<uint64_t>(vid);

  raw.Write<chi_mesh::Vector3>(normal);
  raw.Write<chi_mesh::Vector3>(centroid);
  raw.Write<bool>(has_neighbor);
  raw.Write<uint64_t>(neighbor_id);

  return raw;
}

//###################################################################
/**Deserializes a face from a set of raw data*/
chi_mesh::CellFace chi_mesh::CellFace::
  DeSerialize(const chi_data_types::ByteArray& raw,
              size_t& address)
{

  CellFace face;

  const size_t num_face_verts = raw.Read<size_t>(address, &address);
  face.vertex_ids.reserve(num_face_verts);
  for (size_t fv=0; fv<num_face_verts; ++fv)
    face.vertex_ids.push_back(raw.Read<uint64_t>(address, &address));

  face.normal       = raw.Read<chi_mesh::Vector3>(address, &address);
  face.centroid     = raw.Read<chi_mesh::Vector3>(address, &address);
  face.has_neighbor = raw.Read<bool>             (address, &address);
  face.neighbor_id  = raw.Read<uint64_t>         (address, &address);

  return face;
}


//###################################################################
/**Provides string information of the face.*/
std::string chi_mesh::CellFace::ToString() const
{
  std::stringstream outstr;

  outstr << "num_vertex_ids: " << vertex_ids.size() << "\n";
  {
    size_t counter=0;
    for (uint64_t vid : vertex_ids)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  outstr << "normal: " << normal.PrintS() << "\n";
  outstr << "centroid: " << centroid.PrintS() << "\n";
  outstr << "has_neighbor: " << has_neighbor << "\n";
  outstr << "neighbor_id: " << neighbor_id << "\n";

  return outstr.str();
}


//###################################################################
/**Recomputes the face centroid assuming the mesh vertices
 * have been transformed.*/
void chi_mesh::CellFace::
  RecomputeCentroid(const chi_mesh::MeshContinuum &grid)
{
  centroid = chi_mesh::Vector3(0,0,0);
  for (uint64_t vid : vertex_ids)
    centroid += grid.vertices[vid];
  centroid /= static_cast<double>(vertex_ids.size());
}


//###################################################################
/**Serializes a cell into a vector of bytes.*/
chi_data_types::ByteArray chi_mesh::Cell::Serialize() const
{
  chi_data_types::ByteArray raw;

  raw.Write<uint64_t>(global_id);
  raw.Write<uint64_t>(local_id);
  raw.Write<uint64_t>(partition_id);
  raw.Write<chi_mesh::Vector3>(centroid);
  raw.Write<int>(material_id);

  raw.Write<CellType>(cell_type);
  raw.Write<CellType>(cell_sub_type);

  raw.Write<size_t>(vertex_ids.size());
  for (uint64_t vid : vertex_ids)
    raw.Write<uint64_t>(vid);

  raw.Write<size_t>(faces.size());
  for (const auto& face : faces)
    raw.Append(face.Serialize());

  return raw;
}

//###################################################################
/**Deserializes a cell from a vector of bytes.*/
chi_mesh::Cell chi_mesh::Cell::
  DeSerialize(const chi_data_types::ByteArray& raw,
              size_t& address)
{
  typedef chi_mesh::Vector3 Vec3;
  auto cell_global_id = raw.Read<uint64_t>(address, &address);
  auto cell_local_id  = raw.Read<uint64_t>(address, &address);
  auto cell_prttn_id  = raw.Read<uint64_t>(address, &address);
  auto cell_centroid  = raw.Read<Vec3>    (address, &address);
  auto cell_matrl_id  = raw.Read<int>     (address, &address);

  auto cell_type      = raw.Read<CellType>(address, &address);
  auto cell_sub_type  = raw.Read<CellType>(address, &address);

  Cell cell(cell_type, cell_sub_type);
  cell.global_id    = cell_global_id;
  cell.local_id     = cell_local_id;
  cell.partition_id = cell_prttn_id;
  cell.centroid     = cell_centroid;
  cell.material_id  = cell_matrl_id;

  auto num_vertex_ids = raw.Read<size_t>(address, &address);
  cell.vertex_ids.reserve(num_vertex_ids);
  for (size_t v=0; v<num_vertex_ids; ++v)
    cell.vertex_ids.push_back(raw.Read<uint64_t>(address, &address));

  auto num_faces = raw.Read<size_t>(address, &address);
  cell.faces.reserve(num_faces);
  for (size_t f=0; f<num_faces; ++f)
    cell.faces.push_back(chi_mesh::CellFace::DeSerialize(raw, address));

  return cell;
}


//###################################################################
/**Provides string information of the cell.*/
std::string chi_mesh::Cell::ToString() const
{
  std::stringstream outstr;

  outstr << "cell_type: "     << CellTypeName(cell_type) << "\n";
  outstr << "cell_sub_type: " << CellTypeName(cell_sub_type) << "\n";
  outstr << "global_id: "     << global_id << "\n";
  outstr << "local_id: "      << local_id << "\n";
  outstr << "partition_id: "  << partition_id << "\n";
  outstr << "centroid: "      << centroid.PrintS() << "\n";
  outstr << "material_id: "   << material_id << "\n";

  outstr << "num_vertex_ids: " << vertex_ids.size() << "\n";
  {
    size_t counter=0;
    for (uint64_t vid : vertex_ids)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  {
    outstr << "num_faces: " << faces.size() << "\n";
    size_t f=0;
    for (const auto& face : faces)
      outstr << "Face " << f++ << ":\n" << face.ToString();
  }

  return outstr.str();
}

//###################################################################
/**Recomputes the cell centroid and all face centroids assuming
 * the mesh vertices have been transformed.*/
void chi_mesh::Cell::
  RecomputeCentroidsAndNormals(const chi_mesh::MeshContinuum &grid)
{
  const auto k_hat = Vector3(0,0,1);

  centroid = chi_mesh::Vector3(0,0,0);
  for (uint64_t vid : vertex_ids)
    centroid += grid.vertices[vid];
  centroid /= static_cast<double>(vertex_ids.size());

  for (auto& face : faces)
  {
    face.RecomputeCentroid(grid);

    if (cell_type == CellType::POLYGON)
    {
      const auto v0 = grid.vertices[face.vertex_ids[0]];
      const auto v1 = grid.vertices[face.vertex_ids[1]];

      const auto v01 = v1 - v0;

      face.normal = v01.Cross(k_hat).Normalized();
    }
    else if (cell_type == CellType::POLYHEDRON)
    {
      // A face of a polyhedron can itself be a polygon
      // which can be multifaceted. Here we need the
      // average normal over all the facets computed
      // using an area-weighted average.
      const size_t num_face_verts = face.vertex_ids.size();
      double total_area = 0.0;
      auto weighted_normal = Vector3(0,0,0);
      for (size_t fv=0; fv<num_face_verts; ++fv)
      {
        size_t fvp1 = (fv < (num_face_verts-1))? fv+1 : 0;

        uint64_t fvid_m = face.vertex_ids[fv  ];
        uint64_t fvid_p = face.vertex_ids[fvp1];

        auto leg_m = grid.vertices[fvid_m] - face.centroid;
        auto leg_p = grid.vertices[fvid_p] - face.centroid;

        auto vn = leg_m.Cross(leg_p);

        double area = 0.5*vn.Norm();
        total_area += area;

        weighted_normal = weighted_normal + area*vn.Normalized();
      }
      weighted_normal = weighted_normal/total_area;

      face.normal = weighted_normal.Normalized();
    }
  }
}