#include "pwlc.h"

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polyhedron.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_cylindrical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab_spherical.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_polygon_cylindrical.h"

#include "chi_log.h"

#include "ChiTimer/chi_timer.h"
#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

#define InvalidCoordinateSystem(fname) \
std::invalid_argument((fname) + \
": Unsupported coordinate system type encountered.");

#define UnsupportedCellType(fname) \
std::invalid_argument((fname) + \
": Unsupported cell type encountered.");

//###################################################################
/**Makes a shared_ptr CellPWLView for a cell based on its type.*/
std::shared_ptr<chi_math::CellMappingFE_PWL>
  chi_math::SpatialDiscretization_PWLC::
  MakeCellMappingFE(const chi_mesh::Cell& cell) const
{
  const std::string fname = "SpatialDiscretization_PWLC::MakeCellMapping";

  typedef SlabMappingFE_PWL                SlabSlab;
  typedef SlabMappingFE_PWL_Cylindrical    SlabCyli;
  typedef SlabMappingFE_PWL_Spherical      SlabSphr;
  typedef PolygonMappingFE_PWL             Polygon;
  typedef PolygonMappingFE_PWL_Cylindrical PolygonCyli;
  typedef PolyhedronMappingFE_PWL          Polyhedron;

  using namespace chi_math;
  using namespace std;

  switch (cell.Type())
  {
    case chi_mesh::CellType::SLAB:
    {
      const auto& vol_quad = line_quad_order_arbitrary;

      switch (cs_type)
      {
        case CoordinateSystemType::CARTESIAN:
          return make_shared<SlabSlab>(cell, ref_grid, vol_quad);
        case CoordinateSystemType::CYLINDRICAL:
          return make_shared<SlabCyli>(cell, ref_grid, vol_quad);
        case CoordinateSystemType::SPHERICAL:
          return make_shared<SlabSphr>(cell, ref_grid, vol_quad);
        default:
          throw InvalidCoordinateSystem(fname)
      }
    }
    case chi_mesh::CellType::POLYGON:
    {
      const auto& vol_quad = tri_quad_order_arbitrary;
      const auto& area_quad = line_quad_order_arbitrary;

      switch (cs_type)
      {
        case CoordinateSystemType::CARTESIAN:
          return make_shared<Polygon>(cell, ref_grid, vol_quad, area_quad);
        case CoordinateSystemType::CYLINDRICAL:
          return make_shared<PolygonCyli>(cell, ref_grid, vol_quad,area_quad);
        default:
          throw InvalidCoordinateSystem(fname)
      }
    }
    case chi_mesh::CellType::POLYHEDRON:
    {
      const auto& vol_quad = tet_quad_order_arbitrary;
      const auto& area_quad = tri_quad_order_arbitrary;

      switch (cs_type)
      {
        case CoordinateSystemType::CARTESIAN:
          return make_shared<Polyhedron>(cell,ref_grid,vol_quad,area_quad);
        default:
          throw InvalidCoordinateSystem(fname)
      }
    }
    default:
      throw UnsupportedCellType(fname)
  }
}

//###################################################################
/**Adds a PWL Finite Element for each cell of the local problem.*/
void chi_math::SpatialDiscretization_PWLC::PreComputeCellSDValues()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log0Verbose1() << chi::program_timer.GetTimeString()
                              << " Add cell SD-values.";

  size_t num_local_cells = ref_grid->local_cells.size();

  //================================================== Create empty view
  //                                                 for each cell
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_CELL_MAPPINGS)
    {
      chi::log.Log0Verbose1() << chi::program_timer.GetTimeString()
                                  << " Computing unit integrals.";
      if (!mapping_initialized)
      {
        for (const auto& cell : ref_grid->local_cells)
          cell_mappings.push_back(MakeCellMappingFE(cell));

        mapping_initialized = true;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Unit integrals
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_UNIT_INTEGRALS)
    {
    chi::log.Log0Verbose1() << chi::program_timer.GetTimeString()
                                << " Computing unit integrals.";
      if (not integral_data_initialized)
      {
        fe_unit_integrals.reserve(num_local_cells);
        for (size_t lc=0; lc<num_local_cells; ++lc)
        {
          UIData ui_data;

          auto cell_fe_view = GetCellMappingFE(lc);
          cell_fe_view->ComputeUnitIntegrals(ui_data);

          fe_unit_integrals.push_back(std::move(ui_data));
        }

        integral_data_initialized = true;
      }
    }//if compute unit intgrls
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Quadrature data
  {
    using namespace chi_math::finite_element;
    if (setup_flags & SetupFlags::COMPUTE_QP_DATA)
    {
      chi::log.Log0Verbose1() << chi::program_timer.GetTimeString()
                                  << " Computing quadrature data.";
      if (not qp_data_initialized)
      {
        fe_vol_qp_data.reserve(num_local_cells);
        fe_srf_qp_data.reserve(num_local_cells);
        for (size_t lc=0; lc<num_local_cells; ++lc)
        {
          fe_vol_qp_data.emplace_back();
          fe_srf_qp_data.emplace_back();

          auto cell_fe_view = GetCellMappingFE(lc);
          cell_fe_view->InitializeAllQuadraturePointData(fe_vol_qp_data.back(),
                                                         fe_srf_qp_data.back());
        }

        qp_data_initialized = true;
      }
    }//if init qp data
  }
  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log0Verbose1() << chi::program_timer.GetTimeString()
                              << " Done adding cell SD-values.";

}//AddViewOfLocalContinuum

//###################################################################
/**Returns a locally stored finite element view.*/
std::shared_ptr<chi_math::CellMappingFE_PWL>
  chi_math::SpatialDiscretization_PWLC::GetCellMappingFE(uint64_t cell_local_index)
{

  if (mapping_initialized)
  {
    try
    {
      return cell_mappings.at(cell_local_index);
    }
    catch (const std::out_of_range& o)
    {
      chi::log.LogAllError()
        << "SpatialDiscretization_PWLC::MapFeView "
           "Failure to map Finite Element View. The view is either not"
           "available or the supplied local index is invalid.";
     chi::Exit(EXIT_FAILURE);
    }
  }
  else
  {
    return MakeCellMappingFE(ref_grid->local_cells[cell_local_index]);
  }

  return nullptr;
}
