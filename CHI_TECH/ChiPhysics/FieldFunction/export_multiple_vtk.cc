#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include <vtkInformation.h>

//###################################################################
/**Exports multiple field functions to a VTK file collection.*/
void chi_physics::FieldFunction::
  ExportMultipleFFToVTK(const std::string& file_base_name,
                        std::vector<FieldFunction*> ff_list)
{
  chi_log.Log(LOG_0) << "Exporting field functions to VTK.";

  //============================================= Check ff_list populated
  if (ff_list.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "ExportMultipleFFToVTK: Empty field-function list.";
    exit(EXIT_FAILURE);
  }

  //============================================= Check spatial discretization
  auto ff_type = ff_list.front()->type;
  for (auto ff : ff_list)
    if (ff->type != ff_type)
    {
      chi_log.Log(LOG_ALLERROR)
        << "ExportMultipleFFToVTK: Dissimilar field-function type encountered "
           "in the supplied field-function list. "
           "Currently this function requires "
           "all field-functions used in this call to have the same "
           "spatial discretization.";
      exit(EXIT_FAILURE);
    }
  auto ff_spatial_discretization = ff_list.front()->spatial_discretization;

  //============================================= Check grid
  auto grid = ff_list.front()->spatial_discretization->ref_grid;
  for (auto ff : ff_list)
    if (ff->spatial_discretization->ref_grid != grid)
    {
      chi_log.Log(LOG_ALLERROR)
        << "ExportMultipleFFToVTK: Differing grids encountered "
           "in the supplied field-function list. "
           "Currently this function requires "
           "all field-functions used in this call to refer to the same"
           "grid/mesh.";
      exit(EXIT_FAILURE);
    }

  //============================================= Instantiate grid
  auto ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  //============================================= Populate points
  auto points = vtkSmartPointer<vtkPoints>::New();
  int vc=0;
  for (const auto& cell : grid->local_cells)
    for (int vid : cell.vertex_ids)
    {
      auto vertex = grid->vertices[vid];

      points->InsertPoint(vc++,vertex->x, vertex->y, vertex->z);
    }
  ugrid->SetPoints(points);

  //============================================= Populate cells
  auto material_array         = vtkSmartPointer<vtkIntArray>::New();
  auto partition_number_array = vtkSmartPointer<vtkIntArray>::New();

  material_array->SetName("Material");
  partition_number_array->SetName("Partition");
  vc=0;
  for (const auto& cell : grid->local_cells)
  {
    material_array->InsertNextValue(cell.material_id);
    partition_number_array->InsertNextValue(cell.partition_id);

    //================================= Build cell vertices
    int num_verts = cell.vertex_ids.size();
    std::vector<vtkIdType> vertex_ids;
    vertex_ids.reserve(num_verts);
    for (int vid : cell.vertex_ids)
      vertex_ids.push_back(vc++);

    //================================= Handle cell specific items
    if (cell.Type() == chi_mesh::CellType::SLAB)
      ugrid->
        InsertNextCell(VTK_LINE,
                       num_verts,
                       vertex_ids.data());
    else if (cell.Type() == chi_mesh::CellType::POLYGON)
      ugrid->
        InsertNextCell(VTK_POLYGON,
                       num_verts,
                       vertex_ids.data());
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto faces = vtkSmartPointer<vtkCellArray>::New();

      int num_faces = cell.faces.size();
      for (const auto& face : cell.faces)
      {
        int num_fverts = face.vertex_ids.size();
        std::vector<vtkIdType> fvertex_ids;
        fvertex_ids.reserve(num_fverts);

        for (int fv=0; fv<num_fverts; fv++)
          for (int v=0; v<cell.vertex_ids.size(); ++v)
            if (face.vertex_ids[fv] == cell.vertex_ids[v])
              fvertex_ids.push_back(vertex_ids[v]);
//        for (int vid : face.vertex_ids)
//          fvertex_ids.push_back(vid);

        faces->InsertNextCell(num_fverts,fvertex_ids.data());
      }//for faces

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,
                       num_verts,
                       vertex_ids.data(),
                       num_faces,
                       faces->GetPointer());
    }//polyhedron
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "ExportMultipleFFToVTK: Unsupported cell type encountered.";
      exit(EXIT_FAILURE);
    }
  }//for cell

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_number_array);



  //=============================================
  if (ff_type == FieldFunctionType::FV)
  {
    auto fv = (SpatialDiscretization_FV*)ff_spatial_discretization;
    for (auto ff : ff_list)
    {
      if (ff->unknown_manager== nullptr) continue;

      int ref_unknown = ff->ref_set;
      auto unknown = ff->unknown_manager->unknowns[ref_unknown];

      if (unknown->type == chi_math::UnknownType::SCALAR)
      {
        auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
        unk_arr->SetName(unknown->text_name.c_str());

        for (auto& cell : grid->local_cells)
        {
          int local_mapping =
            fv->MapDOFLocal(&cell,ff->unknown_manager,ref_unknown);

          double value = (*ff->field_vector_local)[local_mapping];

          unk_arr->InsertNextValue(value);
        }

        ugrid->GetCellData()->AddArray(unk_arr);
      }//scalar
    }
  }

  //=============================================
  if (ff_type == FieldFunctionType::CFEM_PWL or
      ff_type == FieldFunctionType::DFEM_PWL)
  {
    for (auto ff : ff_list)
    {
      if (ff->unknown_manager== nullptr) continue;

      int ref_unknown = ff->ref_set;
      auto unknown = ff->unknown_manager->unknowns[ref_unknown];

      int N = ff->unknown_manager->GetTotalUnknownSize();

      if (unknown->type == chi_math::UnknownType::SCALAR)
      {
        int component = ff->unknown_manager->MapUnknown(ref_unknown,0);

        auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
        unk_arr->SetName(unknown->text_name.c_str());

        int c=-1;
        for (auto& cell : grid->local_cells)
        {
          for (int vid : cell.vertex_ids)
          {
            ++c;
            int local_mapping = c*N + component;

            double value = (*ff->field_vector_local)[local_mapping];

            unk_arr->InsertNextValue(value);
          }//for vid
        }//for cell

        ugrid->GetPointData()->AddArray(unk_arr);
      }//scalar
      if (unknown->type == chi_math::UnknownType::VECTOR_N)
      {
        for (int comp=0; comp<unknown->num_components; ++comp)
        {
          int component = ff->unknown_manager->MapUnknown(ref_unknown,comp);

          auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
          unk_arr->SetName(unknown->component_text_names[comp].c_str());

          int c=-1;
          for (auto& cell : grid->local_cells)
          {
            for (int vid : cell.vertex_ids)
            {
              ++c;
              int local_mapping = c*N + component;

              double value = (*ff->field_vector_local)[local_mapping];

              unk_arr->InsertNextValue(value);
            }//for vid
          }//for cell

          ugrid->GetPointData()->AddArray(unk_arr);
        }//for c
      }//scalar
    }
  }

  //============================================= Construct file name
  std::string base_filename     = std::string(file_base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Write master file
  if (chi_mpi.location_id == 0)
  {
    std::string pvtu_file_name = base_filename + std::string(".pvtu");

    auto pgrid_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    pgrid_writer->EncodeAppendedDataOff();
    pgrid_writer->SetFileName(pvtu_file_name.c_str());
    pgrid_writer->SetNumberOfPieces(chi_mpi.process_count);
    pgrid_writer->SetStartPiece(chi_mpi.location_id);
    pgrid_writer->SetEndPiece(chi_mpi.process_count-1);
    pgrid_writer->SetInputData(ugrid);

    pgrid_writer->Write();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Serial output each piece
  auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  chi_log.Log(LOG_0) << "Done exporting field functions to VTK.";
}