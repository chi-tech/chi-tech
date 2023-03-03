#include "chi_meshcontinuum.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_grid_vtk_utils.h"

#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkInformation.h>
#include <vtkModelMetadata.h>

#include <vtkExodusIIWriter.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

#include "chi_runtime.h"
#include "chi_log.h"

#define scvtkid static_cast<vtkIdType>

//###################################################################
/**Exports just the mesh to VTK format.*/
void chi_mesh::MeshContinuum::
  ExportCellsToExodus(const std::string& file_base_name) const
{
  const std::string fname = "chi_mesh::MeshContinuum::ExportCellsToExodus";
  chi::log.Log() << "Exporting mesh to Exodus file with base " << file_base_name;

  if (chi::mpi.process_count != 1)
    throw std::logic_error(fname + ": Currently this routine is only allowed "
                                   "in serial.");

  const auto& grid = *this;

  //============================================= Separate cells by material id
  //                                              and cell-type
  typedef std::pair<int, chi_mesh::CellType> MatIDCellType;
  typedef std::map<MatIDCellType, std::vector<uint64_t>> CType2IDListMap;
  CType2IDListMap mat_id_cell_id_list_map;
  std::map<int, unsigned int> mat_id_usage;
  for (const auto& cell : grid.local_cells)
  {
    mat_id_cell_id_list_map
    [{cell.material_id_, cell.SubType()}].push_back(cell.local_id_);
    mat_id_usage[cell.material_id_] = 0;
  }

  //============================================= Create unstructured meshes
  //                                              for each material-type pair
  vtkNew<vtkMultiBlockDataSet> grid_blocks;
  {
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkPoints> points;

    vtkNew<vtkIdTypeArray> global_node_id_list;
    global_node_id_list->SetName("GlobalNodeId");

    vtkNew<vtkIdTypeArray> global_elem_id_list;
    global_elem_id_list->SetName("GlobalElementId");

    vtkNew<vtkIntArray> block_id_list;
    block_id_list->SetName("BlockID");

    //========================== Build vertex set
    std::set<uint64_t> vid_set;
    for (const auto& cell : grid.local_cells)
      for (uint64_t vid : cell.vertex_ids_)
        vid_set.insert(vid);

    //========================== Build vertex map
    std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
    {
      uint64_t mapped_id = 0;
      for (uint64_t vid : vid_set)
        vertex_map[vid] = mapped_id++;
    }

    //========================== Load vertices
    for (uint64_t vid : vid_set)
    {
      const auto& vertex = grid.vertices[vid];
      points->InsertNextPoint(vertex.x,vertex.y,vertex.z);
      global_node_id_list->InsertNextValue(scvtkid(vid+1));
    }

    //========================== Load cells

    for (const auto& cell : grid.local_cells)
    {
      chi_mesh::UploadCellGeometryContinuous(cell, vertex_map, ugrid);
      global_elem_id_list->InsertNextValue(scvtkid(cell.global_id_+1));
      block_id_list->InsertNextValue(cell.material_id_);
    }//for local cells

    //========================== Set arrays
    ugrid->SetPoints(points);
    ugrid->GetPointData()->AddArray(global_node_id_list);
    ugrid->GetCellData()->AddArray(global_elem_id_list);
    ugrid->GetCellData()->AddArray(block_id_list);

    ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");
    ugrid->GetCellData()->SetActiveGlobalIds("GlobalElementId");

    //========================== Set block
    grid_blocks->SetBlock(0, ugrid);

    chi::log.Log()
      << "Writing grid block "
      << " Number of cells: " << ugrid->GetNumberOfCells()
      << " Number of points: " << ugrid->GetNumberOfPoints();
  }//end of grid_blocks assignment

  //============================================= Separate faces by boundary
  //                                              id
  struct FaceInfo
  {
    chi_mesh::CellFace const* face_ptr;
    uint64_t source_cell_id;
    int source_face_id;
  };
  typedef std::vector<FaceInfo> ListOfFaces;
  std::map<uint64_t, ListOfFaces> boundary_id_faces_map;
  for (const auto& cell : grid.local_cells)
  {
    int f=0;
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
        boundary_id_faces_map[face.neighbor_id_].push_back(
          {&face,cell.global_id_,f});
      ++f;
    }
  }

  //============================================= Make
  vtkNew<vtkMultiBlockDataSet> nodesets_blocks;
  vtkNew<vtkMultiBlockDataSet> sidesets_blocks;
  for (const auto& [bndry_id, face_list] : boundary_id_faces_map)
  {
    const std::string block_name = grid.GetBoundaryIDMap().at(bndry_id);

    //====================================== NodeSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      vtkNew<vtkIdTypeArray> node_global_ids;
      node_global_ids->SetName("GlobalNodeId");

      //========================== Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      //========================== Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      //========================== Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x,vertex.y,vertex.z);
        node_global_ids->InsertNextValue(scvtkid(vid+1));
      }

      //========================== Load cells
      for (uint64_t vid : vid_set)
      {
        std::vector<vtkIdType> cell_vids = {scvtkid(vertex_map[vid])};
        ugrid->InsertNextCell(VTK_VERTEX,
                              scvtkid(1),
                              cell_vids.data());
      }

      ugrid->SetPoints(points);
      ugrid->GetPointData()->AddArray(node_global_ids);

      ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");

      nodesets_blocks->SetBlock(bndry_id, ugrid);
      nodesets_blocks->GetMetaData(bndry_id)->
        Set(vtkCompositeDataSet::NAME(), block_name);

      chi::log.Log()
      << "Writing nodeset block " << block_name
      << " Number of cells: " << ugrid->GetNumberOfCells()
      << " Number of points: " << ugrid->GetNumberOfPoints();
    }

    //====================================== SideSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      vtkNew<vtkIdTypeArray> src_cell_global_ids;
      vtkNew<vtkIntArray>    src_cell_face_id;
      src_cell_global_ids->SetName("SourceElementId");
      src_cell_face_id->SetName("SourceElementSide");

      //========================== Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      //========================== Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      //========================== Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x,vertex.y,vertex.z);
      }

      //========================== Load faces
      for (const auto& face_info : face_list)
      {
        UploadFaceGeometry(*face_info.face_ptr, vertex_map, ugrid);
        src_cell_global_ids->InsertNextValue(scvtkid(face_info.source_cell_id));
        src_cell_face_id->InsertNextValue(face_info.source_face_id);
      }

      ugrid->SetPoints(points);
      ugrid->GetCellData()->AddArray(src_cell_global_ids);
      ugrid->GetCellData()->AddArray(src_cell_face_id);

      sidesets_blocks->SetBlock(bndry_id, ugrid);
      sidesets_blocks->GetMetaData(bndry_id)->
        Set(vtkCompositeDataSet::NAME(), block_name);

      chi::log.Log()
        << "Writing sideset block " << block_name
        << " Number of cells: " << ugrid->GetNumberOfCells()
        << " Number of points: " << ugrid->GetNumberOfPoints();
    }//End of side-set
  }

  //============================================= Write the file
  vtkNew<vtkMultiBlockDataSet> main_block;
  main_block->SetBlock(0, grid_blocks);
  main_block->SetBlock(1, nodesets_blocks);
  main_block->GetMetaData(1)->Set(vtkCompositeDataSet::NAME(), "Node Sets");
  main_block->SetBlock(2, sidesets_blocks);
  main_block->GetMetaData(2)->Set(vtkCompositeDataSet::NAME(), "Side Sets");

  vtkNew<vtkExodusIIWriter> writer;
  writer->SetBlockIdArrayName("BlockID");

  writer->SetFileName((file_base_name + ".e").c_str());

  writer->SetInputData(main_block);

  writer->WriteOutGlobalNodeIdArrayOff();
  writer->WriteOutGlobalElementIdArrayOff();
  writer->WriteOutBlockIdArrayOff();

  writer->Write();

  auto em = writer->GetModelMetadata();

  chi::log.Log() << "Num Blocks   :  " << em->GetNumberOfBlocks();
  chi::log.Log() << "Num Node Sets:  " << em->GetNumberOfNodeSets();
  chi::log.Log() << "Num Side Sets:  " << em->GetNumberOfSideSets();

  chi::log.Log() << "Done exporting mesh to VTK.";
  MPI_Barrier(MPI_COMM_WORLD);
}