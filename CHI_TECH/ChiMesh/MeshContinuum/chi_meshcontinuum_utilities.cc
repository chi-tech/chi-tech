#include "chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell_slab.h"

#include <boost/graph/bandwidth.hpp>

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

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
  for (auto c : local_cell_glob_indices)
  {
    auto cell = cells[c];

    for (auto face : cell->faces)
      face_size_histogram.push_back(face.vertex_ids.size());
  }
  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  //================================================== Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  //================================================== Compute average and ratio
  size_t smallest_face = face_size_histogram.front();
  size_t largest_face = face_size_histogram.back();
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
      running_average = (double)running_total_face_dofs/running_face_count;
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
/**Check whether a cell is local*/
bool chi_mesh::MeshContinuum::IsCellLocal(int cell_global_index)
{
  if (cell_global_index<0)
  {
    return false;
  }
  else
  {
    auto cell = cells[cell_global_index];
    if (cell->partition_id == chi_mpi.location_id)
    {
      return true;
    } else
    {
      return false;
    }

  }
  return false;
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

  size_t face_dof_size = 0;

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
/**Check whether a cell is a boundary*/
bool chi_mesh::MeshContinuum::IsCellBndry(int cell_global_index)
{
  if (cell_global_index<0)
    return true;

  return false;
}



//###################################################################
/**Generalized find associated cell face*/
int chi_mesh::MeshContinuum::FindAssociatedFace(chi_mesh::CellFace& cur_face,
                                                int adj_cell_g_index,bool verbose)
{
  //======================================== Check index validity
  if (IsCellBndry(adj_cell_g_index) || (!IsCellLocal(adj_cell_g_index)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedFace. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check cell validity by index
//  if (adj_cell_g_index >= cells.size())
//  {
//    chi_log.Log(LOG_ALLERROR)
//      << "Invalid cell index encountered in call to "
//      << "MeshContinuum::FindAssociatedFace. Index is out of cell index bounds.";
//    exit(EXIT_FAILURE);
//  }

  chi_mesh::Cell* adj_cell = cells[adj_cell_g_index];

  int associated_face = -1;

  //======================================== Loop over adj cell faces
  for (int af=0; af < adj_cell->faces.size(); af++)
  {
    //Assume face matches
    bool face_matches = true; //Now disprove it
    //================================= Loop over adj cell face verts
    for (int afv=0; afv < adj_cell->faces[af].vertex_ids.size(); afv++)
    {
      //========================== Try and find them in the reference face
      bool found = false;
      for (int cfv=0; cfv<cur_face.vertex_ids.size(); cfv++)
      {
        if (cur_face.vertex_ids[cfv] == adj_cell->faces[af].vertex_ids[afv])
        {
          found = true;
          break;
        }
      }//for cfv

      if (!found) {face_matches = false; break;}
    }//for afv

    if (face_matches) {associated_face = af; break;}
  }


  //======================================== Check associated face validity
  if (associated_face<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Could not find associated face in call to "
      << "MeshContinuum::FindAssociatedFace. Reference face with centroid at \n"
      << cur_face.centroid.PrintS();
    for (int af=0; af < adj_cell->faces.size(); af++)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Adjacent cell face " << af << " centroid "
        << adj_cell->faces[af].centroid.PrintS();
    }
    exit(EXIT_FAILURE);
  }

  //======================================== Verbose output
  if (verbose)
  {
    std::stringstream out_string;

    out_string
      << "Adj cell " << adj_cell_g_index << ":\n"
      << "face " << associated_face << "\n";
    for (int v=0; v < adj_cell->faces[associated_face].vertex_ids.size(); v++)
    {
      out_string
        << "vertex " << v << " "
        << adj_cell->faces[associated_face].vertex_ids[v]
        << "\n";
    }
    out_string
      << "Cur cell "
      << adj_cell->faces[associated_face].neighbor << ":\n";
    for (int v=0; v<cur_face.vertex_ids.size(); v++)
    {
      out_string
        << "vertex " << v << " "
        << cur_face.vertex_ids[v]
        << "\n";
    }
    chi_log.Log(LOG_ALL) << out_string.str();
  }

  return associated_face;
}



//###################################################################
/**General map vertices*/
void chi_mesh::MeshContinuum::
FindAssociatedVertices(chi_mesh::CellFace& cur_face,
                       int adj_cell_g_index, int associated_face,
                       std::vector<int>& dof_mapping)
{
  //======================================== Check index validity
  if (IsCellBndry(adj_cell_g_index) || (!IsCellLocal(adj_cell_g_index)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedVertices. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check cell validity by index
//  if (adj_cell_g_index >= cells.size())
//  {
//    chi_log.Log(LOG_ALLERROR)
//      << "Invalid cell index encountered in call to "
//      << "MeshContinuum::FindAssociatedVertices. Index is out of cell index bounds.";
//    exit(EXIT_FAILURE);
//  }

  chi_mesh::Cell* adj_cell = cells[adj_cell_g_index];


  for (int cfv=0; cfv<cur_face.vertex_ids.size(); cfv++)
  {
    bool found = false;
    for (int afv=0;
         afv < adj_cell->faces[associated_face].vertex_ids.size(); afv++)
    {
      if (cur_face.vertex_ids[cfv] ==
        adj_cell->faces[associated_face].vertex_ids[afv])
      {
        dof_mapping.push_back(afv);
        found = true;
        break;
      }
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Face DOF mapping failed in call to "
        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
           "node."
        << adj_cell_g_index << " " << cur_face.centroid.PrintS();
      exit(EXIT_FAILURE);
    }

  }//for cfv

}


//###################################################################
/**Computes the centroid from nodes specified by the given list.*/
chi_mesh::Vector3 chi_mesh::MeshContinuum::
  ComputeCentroidFromListOfNodes(const std::vector<int> &list)
{
  if (list.empty())
  {
    chi_log.Log(LOG_ALLERROR) << "ComputeCentroidFromListOfNodes, empty list";
    exit(EXIT_FAILURE);
  }
  chi_mesh::Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + *vertices[node_id];

  return centroid/list.size();
}