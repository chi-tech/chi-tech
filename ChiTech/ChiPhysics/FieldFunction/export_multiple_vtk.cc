#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>

#include <vtkInformation.h>

//###################################################################
/**Exports multiple field functions to a VTK file collection.*/
void chi_physics::FieldFunction::
  ExportMultipleFFToVTK(
    const std::string& file_base_name,
    const std::vector<std::shared_ptr<chi_physics::FieldFunction>>& ff_list)
{
  const std::string fname = __FUNCTION__;
  chi::log.Log() << "Exporting field functions to VTK with file base \""
                     << file_base_name << "\"";

  //============================================= Check ff_list populated
  if (ff_list.empty())
  {
    chi::log.LogAllError()
      << "ExportMultipleFFToVTK: Empty field-function list.";
    chi::Exit(EXIT_FAILURE);
  }

  //============================================= Check spatial discretizations
  // All the spatial discretizations need to be
  // the same
  auto ff_type = ff_list.front()->spatial_discretization->type;
  for (auto& ff : ff_list)
    if (ff->spatial_discretization->type != ff_type)
    {
      chi::log.LogAllError()
        << "ExportMultipleFFToVTK: Dissimilar field-function type encountered "
           "in the supplied field-function list. "
           "Currently this function requires "
           "all field-functions used in this call to have the same "
           "spatial discretization.";
      chi::Exit(EXIT_FAILURE);
    }
  auto ff_sdm = ff_list.front()->spatial_discretization;
  const auto ff_sdm_type = ff_sdm->type;

  //============================================= Check grid
  // All of the field-functions need to refer
  // to the same grid
  auto grid = ff_list.front()->spatial_discretization->ref_grid;
  for (auto& ff : ff_list)
    if (ff->spatial_discretization->ref_grid != grid)
    {
      chi::log.LogAllError()
        << "ExportMultipleFFToVTK: Differing grids encountered "
           "in the supplied field-function list. "
           "Currently this function requires "
           "all field-functions used in this call to refer to the same"
           "grid/mesh.";
      chi::Exit(EXIT_FAILURE);
    }

  //============================================= Determine cell/point data
  typedef chi_math::SpatialDiscretizationType SDMType;

  enum class CellDataScope
  {
    CELLDATA = 0,
    POINTDATA = 1
  };

  CellDataScope ff_data_scope;
  switch (ff_sdm_type)
  {
    case SDMType::FINITE_VOLUME:
      ff_data_scope = CellDataScope::CELLDATA;
      break;
    case SDMType::PIECEWISE_LINEAR_CONTINUOUS:
    case SDMType::PIECEWISE_LINEAR_DISCONTINUOUS:
      ff_data_scope = CellDataScope::POINTDATA;
      break;
    default:
      throw std::logic_error(fname + ": Unsupported spatial discretization "
                                     "method encountered.");
  }

  //============================================= Instantiate VTK items
  vtkNew<vtkUnstructuredGrid>         ugrid;
  vtkNew<vtkPoints>                   points;
  vtkNew<vtkIntArray>                 material_array;
  vtkNew<vtkUnsignedIntArray>         partition_id_array;

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  //############################################# Populate cell information
  int64_t node_count=0;
  for (const auto& cell : grid->local_cells)
  {
    UploadCellGeometry(*grid, cell, node_count, points, ugrid);

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);
  }//for local cells
  ugrid->SetPoints(points);

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);

  //============================================= Upload cell/point data
  auto cell_data = ugrid->GetCellData();
  auto point_data = ugrid->GetPointData();
  for (const auto &ff: ff_list)
  {
    unsigned int ref_unknown_id = ff->ref_variable;
    const auto &unknown = ff->unknown_manager.unknowns[ref_unknown_id];

    for (unsigned int comp=0; comp<unknown.num_components; ++comp)
    {
      vtkNew<vtkDoubleArray> unk_arr;
      unk_arr->SetName((ff->text_name +
                        unknown.text_name+
                        unknown.component_text_names[comp]).c_str());

      //Populate the array here
      for (const auto& cell : grid->local_cells)
      {
        size_t num_nodes = ff_sdm->GetCellNumNodes(cell);
        for (int n=0; n<num_nodes; ++n)
        {
          const int64_t nmap = ff_sdm->MapDOFLocal(cell,n,ff->unknown_manager,
                                                   ref_unknown_id,
                                                   comp);
          unk_arr->InsertNextValue((*ff->field_vector_local)[nmap]);
        }//for node
      }//for cell

      switch (ff_data_scope)
      {
        case CellDataScope::CELLDATA:  cell_data->AddArray(unk_arr); break;
        case CellDataScope::POINTDATA: point_data->AddArray(unk_arr);break;
      }
    }//for component
  }//for ff

  //============================================= Construct file name
  std::string base_filename     = std::string(file_base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi::mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Write master file
  if (chi::mpi.location_id == 0)
  {
    std::string pvtu_file_name = base_filename + std::string(".pvtu");

    auto pgrid_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    pgrid_writer->EncodeAppendedDataOff();
    pgrid_writer->SetFileName(pvtu_file_name.c_str());
    pgrid_writer->SetNumberOfPieces(chi::mpi.process_count);
    pgrid_writer->SetStartPiece(chi::mpi.location_id);
    pgrid_writer->SetEndPiece(chi::mpi.process_count-1);
    pgrid_writer->SetInputData(ugrid);

    pgrid_writer->Write();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Serial output each piece
  auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  chi::log.Log() << "Done exporting field functions to VTK.";
}