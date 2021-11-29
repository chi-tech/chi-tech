#include "chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_log.h"
extern ChiLog&  chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <algorithm>

//###################################################################
/**Populates a face histogram.
 *
 * \param master_tolerance Multiple histograms will only be attempted
 * if the ratio of the maximum dofs-per-face to the average dofs-per-face
 * is greater than this value. Default 1.2.
 *
 * \param slave_tolerance While traversing a sorted list of dofs-per-face,
 * a new bin will only be generated when the ratio of the listed dofs-per-face
 * to a running bin average exceeds this value. Defualt 1.1.
 *
 * The function populates face_categories which is a structure containing
 * pairs. Pair.first is the max dofs-per-face for the category and Pair.second
 * is the number of faces in this category.
 *
 * */
void chi_mesh::MeshContinuum::
  BuildFaceHistogramInfo(double master_tolerance, double slave_tolerance)
{
  if (face_histogram_available) return;

  //================================================== Fill histogram
  std::vector<size_t> face_size_histogram;
  for (const auto& cell : local_cells)
    for (const auto& face : cell.faces)
      face_size_histogram.push_back(face.vertex_ids.size());

  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  //================================================== Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  //================================================== Compute average and ratio
  size_t smallest_face   = face_size_histogram.front();
  size_t largest_face    = face_size_histogram.back();
  size_t total_num_faces = face_size_histogram.size();
  double average_dofs_per_face = (double)total_face_dofs_count/total_num_faces;

  std::stringstream outstr;
  outstr << "\nSmallest face = " << smallest_face;
  outstr << "\nLargest face = " << largest_face;
  outstr << "\nTotal face dofs = " << total_face_dofs_count;
  outstr << "\nTotal faces = " << face_size_histogram.size();
  outstr << "\nAverage dofs/face = " << average_dofs_per_face;
  outstr << "\nMax to avg ratio = " << largest_face/average_dofs_per_face;
  chi_log.Log(LOG_ALLVERBOSE_1) << outstr.str();

  //================================================== Determine number of bins
  size_t last_bin_num_faces = total_num_faces;
  if ((largest_face/average_dofs_per_face) > master_tolerance)
  {
    chi_log.Log(LOG_ALLVERBOSE_1)
    << "The ratio of max face dofs to average face dofs "
    << "is larger than " << master_tolerance
    << ", therefore a binned histogram "
    << "will be constructed.";

    //====================================== Build categories
    size_t running_total_face_dofs = 0;
    size_t running_face_count = 0;
    size_t running_face_size = face_size_histogram[0];

    double running_average = face_size_histogram[0];

    for (size_t f=0; f<total_num_faces; ++f)
    {
      if ((face_size_histogram[f]/running_average) > slave_tolerance)
      {
        face_categories.emplace_back(running_face_size,running_face_count);
        running_total_face_dofs = 0;
        running_face_count = 0;
      }

      running_face_size = face_size_histogram[f];
      running_total_face_dofs += face_size_histogram[f];
      running_face_count++;
      running_average = (double)running_total_face_dofs/double(running_face_count);
      last_bin_num_faces = running_face_count;
    }
  }
  face_categories.emplace_back(largest_face,last_bin_num_faces);

  //================================================== Verbose print bins
  outstr.str(std::string());
  outstr
  << "A total of " << face_categories.size()
  << " bins were created:\n";

  size_t bin_counter = -1;
  for (auto bins : face_categories)
  {
    outstr
    << "Bin " << ++bin_counter << ": "
    << bins.second << " faces with max face dofs " << bins.first << "\n";
  }

  chi_log.Log(LOG_ALLVERBOSE_1) << outstr.str();

  face_histogram_available = true;
}

//###################################################################
/**Gets the number of face-histogram categories.*/
size_t chi_mesh::MeshContinuum::NumberOfFaceHistogramBins()
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  return face_categories.size();
}

//###################################################################
/**Maps the face-histogram category number for a given face size.*/
size_t chi_mesh::MeshContinuum::MapFaceHistogramBins(size_t num_face_dofs)
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  size_t category_counter = -1;
  for (auto category : face_categories)
  {
    category_counter++;
    if (num_face_dofs <= category.first)
      return category_counter;
  }

  return 0;
}

//###################################################################
/**Maps the face-histogram category number for a given face size.*/
size_t chi_mesh::MeshContinuum::GetFaceHistogramBinDOFSize(size_t category)
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  size_t face_dof_size;

  try {
    face_dof_size = face_categories.at(category).first;
  }
  catch (std::out_of_range& o){
    chi_log.Log(LOG_ALLWARNING)
    << "Fault detected in chi_mesh::MeshContinuum::"
    << "GetFaceHistogramBinDOFSize.";
    return 0;
  }

  return face_dof_size;
}

//###################################################################
/**Check whether a cell is local by attempting to find the key in
 * the native index map.*/
bool chi_mesh::MeshContinuum::IsCellLocal(uint64_t cell_global_index) const
{
  auto native_index = global_cell_id_to_native_id_map.find(cell_global_index);

  if (native_index != global_cell_id_to_native_id_map.end())
    return true;

  return false;
}


//###################################################################
/**Check whether a cell is a boundary by checking if the key is
 * found in the native or foreign cell maps.*/
bool chi_mesh::MeshContinuum::IsCellBndry(uint64_t cell_global_index) const
{
  auto native_index = global_cell_id_to_native_id_map.find(cell_global_index);
  auto foreign_index = global_cell_id_to_foreign_id_map.find(cell_global_index);

  auto no_native = global_cell_id_to_native_id_map.end();
  auto no_foreign = global_cell_id_to_foreign_id_map.end();

  if ( (native_index == no_native) and (foreign_index == no_foreign))
    return true;

  return false;
}



//###################################################################
/**General map vertices*/
void chi_mesh::MeshContinuum::
FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                       std::vector<short>& dof_mapping) const
{
  int associated_face = cur_face.GetNeighborAssociatedFace(*this);
  //======================================== Check index validity
  if ((not cur_face.has_neighbor) || (not cur_face.IsNeighborLocal(*this)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedVertices. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  auto& adj_cell = local_cells[cur_face.GetNeighborLocalID(*this)];

  dof_mapping.reserve(cur_face.vertex_ids.size());

  const auto& adj_face = adj_cell.faces[associated_face];

  for (auto cfvid : cur_face.vertex_ids)
  {
    bool found = false;
    short afv = 0;
    for (auto afvid : adj_face.vertex_ids)
    {
      if (cfvid == afvid)
      {
        dof_mapping.push_back((short)afv);
        found = true;
        break;
      }
      afv++;
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Face DOF mapping failed in call to "
        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
           "node."
        << cur_face.neighbor_id << " " << cur_face.centroid.PrintS();
      exit(EXIT_FAILURE);
    }
  }

}


//###################################################################
/**Computes the centroid from nodes specified by the given list.*/
chi_mesh::Vector3 chi_mesh::MeshContinuum::
  ComputeCentroidFromListOfNodes(const std::vector<uint64_t> &list) const
{
  if (list.empty())
  {
    chi_log.Log(LOG_ALLERROR) << "ComputeCentroidFromListOfNodes, empty list";
    exit(EXIT_FAILURE);
  }
  chi_mesh::Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + vertices[node_id];

  return centroid/double(list.size());
}

//###################################################################
/**Counts the number of cells within a logical volume across all
 * partitions.*/
size_t chi_mesh::MeshContinuum::
  CountCellsInLogicalVolume(chi_mesh::LogicalVolume &log_vol) const
{
  size_t local_count=0;
  for (const auto& cell : local_cells)
    if (log_vol.Inside(cell.centroid))
      ++local_count;

  size_t global_count=0;

  MPI_Allreduce(&local_count,           //sendbuf
                &global_count,          //recvbuf
                1,                      //count
                MPI_UNSIGNED_LONG_LONG, //datatype
                MPI_SUM,                //op
                MPI_COMM_WORLD);        //communicator

  return global_count;
}