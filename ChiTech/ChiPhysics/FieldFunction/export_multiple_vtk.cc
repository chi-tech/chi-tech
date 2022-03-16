#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include <vtkInformation.h>

//###################################################################
/**Exports multiple field functions to a VTK file collection.*/
void chi_physics::FieldFunction::
  ExportMultipleFFToVTK(const std::string& file_base_name,
                        const std::vector<std::shared_ptr<chi_physics::FieldFunction>>& ff_list)
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
  auto ff_type = ff_list.front()->spatial_discretization->type;
  for (auto& ff : ff_list)
    if (ff->spatial_discretization->type != ff_type)
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
  for (auto& ff : ff_list)
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

  //============================================= Instantiate VTK grid
  auto ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  //============================================= Populate VTK points
  auto points = vtkSmartPointer<vtkPoints>::New();
  int64_t vertex_counter=0;
  for (const auto& cell : grid->local_cells)
    for (uint64_t vid : cell.vertex_ids)
    {
      auto vertex = grid->vertices[vid];

      points->InsertPoint(vertex_counter++, vertex.x, vertex.y, vertex.z);
    }
  ugrid->SetPoints(points);

  //============================================= Populate VTK cells
  auto material_array         = vtkSmartPointer<vtkIntArray>::New();
  auto partition_number_array = vtkSmartPointer<vtkIntArray>::New();

  material_array->SetName("Material");
  partition_number_array->SetName("Partition");
  vertex_counter=0;
  for (const auto& cell : grid->local_cells)
  {
    material_array->InsertNextValue(cell.material_id);
    partition_number_array->InsertNextValue(static_cast<int>(cell.partition_id));

    //================================= Build cell vertices
    size_t num_verts = cell.vertex_ids.size();
    std::vector<vtkIdType> vertex_ids;
    vertex_ids.reserve(num_verts);
    for (auto& vid : cell.vertex_ids)
      vertex_ids.push_back(vertex_counter++);

    //================================= Handle cell specific items
    if (cell.Type() == chi_mesh::CellType::SLAB)
      ugrid->
        InsertNextCell(VTK_LINE,
                       static_cast<vtkIdType>(num_verts),
                       vertex_ids.data());
    else if (cell.Type() == chi_mesh::CellType::POLYGON)
      ugrid->
        InsertNextCell(VTK_POLYGON,
                       static_cast<vtkIdType>(num_verts),
                       vertex_ids.data());
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      vtkNew<vtkIdList> faces;

      size_t num_faces = cell.faces.size();
      for (const auto& face : cell.faces)
      {
        size_t num_fverts = face.vertex_ids.size();
        std::vector<vtkIdType> fvertex_ids;
        fvertex_ids.reserve(num_fverts);

        for (int fv=0; fv<num_fverts; fv++)
          for (int v=0; v<cell.vertex_ids.size(); ++v)
            if (face.vertex_ids[fv] == cell.vertex_ids[v])
              fvertex_ids.push_back(vertex_ids[v]);

        faces->InsertNextId(static_cast<vtkIdType>(num_fverts));
        for (auto vid : fvertex_ids)
          faces->InsertNextId(vid);
      }//for faces

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,
                       static_cast<vtkIdType>(num_verts),
                       vertex_ids.data(),
                       static_cast<vtkIdType>(num_faces),
                       faces->GetPointer(0));
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
  typedef chi_math::SpatialDiscretizationType SDMType;
  typedef SpatialDiscretization_FV SDMFV;

  if (ff_type == SDMType::FINITE_VOLUME)
  {
    auto fv_ptr = std::dynamic_pointer_cast<SDMFV>(ff_spatial_discretization);
    if (not fv_ptr) throw std::logic_error(std::string(__FUNCTION__) +
                                       ": Failed to obtain fv-sdm.");
    auto& fv = *fv_ptr;
    for (auto& ff : ff_list)
    {
      unsigned int ref_unknown = ff->ref_variable;
      const auto& unknown = ff->unknown_manager.unknowns[ref_unknown];

      if (unknown.type == chi_math::UnknownType::SCALAR)
      {
        auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
        unk_arr->SetName(unknown.text_name.c_str());

        for (auto& cell : grid->local_cells)
        {
          int64_t local_mapping =
            fv.MapDOFLocal(cell,0,ff->unknown_manager,ref_unknown,0);

          double value = (*ff->field_vector_local)[local_mapping];

          unk_arr->InsertNextValue(value);
        }

        ugrid->GetCellData()->AddArray(unk_arr);
      }//scalar
    }
  }

  //=============================================
  if (ff_type == SDMType::PIECEWISE_LINEAR_CONTINUOUS or
      ff_type == SDMType::PIECEWISE_LINEAR_DISCONTINUOUS)
  {
    int ff_number = -1;
    for (auto& ff : ff_list)
    {
      unsigned int ref_unknown = ff->ref_variable;
      const auto& unknown = ff->unknown_manager.unknowns[ref_unknown];
      ff_number++;

      unsigned int N = ff->unknown_manager.GetTotalUnknownStructureSize();

      if (unknown.type == chi_math::UnknownType::SCALAR)
      {
        unsigned int component = ff->unknown_manager.MapUnknown(ref_unknown, 0);

        auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
        if (unknown.text_name.empty())
          unk_arr->SetName((std::string("FF_")+
                            std::to_string(ff_number)).c_str());
        else
          unk_arr->SetName((std::string("FF_") +
                            std::to_string(ff_number) +
                            unknown.text_name).c_str());

        uint64_t c=0;
        for (auto& cell : grid->local_cells)
        {
          for (int v=0; v < cell.vertex_ids.size(); ++v)
          {
            uint64_t local_mapping = c*N + component;
            ++c;

            double value = (*ff->field_vector_local)[local_mapping];

            unk_arr->InsertNextValue(value);
          }//for vid
        }//for cell

        ugrid->GetPointData()->AddArray(unk_arr);
      }//scalar
      if (unknown.type == chi_math::UnknownType::VECTOR_2 or
          unknown.type == chi_math::UnknownType::VECTOR_3 or
          unknown.type == chi_math::UnknownType::VECTOR_N)
      {
        for (int comp=0; comp<unknown.num_components; ++comp)
        {
          unsigned int component = ff->unknown_manager.MapUnknown(ref_unknown, comp);

          auto unk_arr = vtkSmartPointer<vtkDoubleArray>::New();
          if (unknown.component_text_names[comp].empty())
            unk_arr->SetName((std::string("FF_") +
                              std::to_string(ff_number) +
                              std::string("Component_") +
                              std::to_string(comp)).c_str());
          else
            unk_arr->SetName((std::string("FF_") +
                              std::to_string(ff_number) +
                              unknown.component_text_names[comp]).c_str());

          uint64_t c=0;
          for (auto& cell : grid->local_cells)
          {
            for (size_t v=0; v < cell.vertex_ids.size(); ++v)
            {
              uint64_t local_mapping = c*N + component;
              ++c;

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