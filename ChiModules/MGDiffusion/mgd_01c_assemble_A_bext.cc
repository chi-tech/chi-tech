#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

//#include "ChiMesh/MeshHandler/chi_meshhandler.h"
//
//#include "mg_diffusion_bndry.h"
//
//#include "ChiPhysics/FieldFunction/fieldfunction.h"
//

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"
//#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
//#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Assemble_A_bext()
{
  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  //============================================= Assemble the system
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();

    const auto& xs   = matid_to_xs_map.at(cell.material_id);
    const auto& qext = matid_to_src_map.at(cell.material_id);

    const size_t num_nodes = cell_mapping.NumNodes();

    std::vector<MatDbl> Acell;
    Acell.resize(mg_diffusion::Solver::num_groups);
    std::vector<VecDbl> rhs_cell;
    rhs_cell.resize(mg_diffusion::Solver::num_groups);

    for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
    {
      Acell[g].resize(num_nodes, VecDbl(num_nodes, 0.0));
      rhs_cell[g].resize(num_nodes, 0.0);
    }

 
    for (size_t i=0; i<num_nodes; ++i)
    {
      for (size_t j=0; j<num_nodes; ++j)
      {
        double entry_mij = 0.0;
        double entry_kij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_mij +=  qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp)
                        * qp_data.JxW(qp);
          entry_kij +=  qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp))
                        * qp_data.JxW(qp);
        }//for qp
        for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
        {
          const double Dg     = xs->diffusion_coeff[g];
          const double sigr_g = xs->sigma_removal[g];
          Acell[g][i][j] = entry_mij * sigr_g + entry_kij * Dg;
        }
      }//for j
      double entry_rhsi = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        entry_rhsi += qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
        rhs_cell[g][i] = entry_rhsi * (qext->source_value_g[g]);
    }//for i
 
    //======================= Flag nodes for being on a boundary
    std::vector<int> dirichlet_count(num_nodes, 0);
    std::vector<double> dirichlet_value(num_nodes, 0.0);

    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      // not a boundary face
	    if (face.has_neighbor) continue; 
	  
      const auto& bndry = boundaries[face.neighbor_id];

      // Robin boundary
      if (bndry.type == BoundaryType::Robin)
      { 
        const auto  qp_face_data = cell_mapping.MakeFaceQuadraturePointData( f );
        const size_t num_face_nodes = face.vertex_ids.size();

        const auto& aval = bndry.mg_values[0];
        const auto& bval = bndry.mg_values[1];
        const auto& fval = bndry.mg_values[2];

//        chi::log.Log() << "Boundary  set as Robin with a,b,f = ("
//                    << aval << ","
//                    << bval << ","
//                    << fval << ") ";
        // true Robin when a!=0, otherwise, it is a Neumann:
        // Assert if b=0
        if (std::fabs(bval) < 1e-8)
          throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");
        
        // loop over nodes of that face
        for (size_t fi=0; fi<num_face_nodes; ++fi)
        {
          const uint i = cell_mapping.MapFaceNode(f,fi);
            
          double entry_rhsi = 0.0;
          for (size_t qp : qp_face_data.QuadraturePointIndices() )
            entry_rhsi +=  qp_face_data.ShapeValue(i, qp) * qp_face_data.JxW(qp);
          for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
            rhs_cell[g][i] +=  fval / bval * entry_rhsi;
            
          // only do this part if true Robin (i.e., a!=0)
          if (std::fabs(aval) > 1.0e-8)
          {
            for (size_t fj=0; fj<num_face_nodes; ++fj)
            {
              const uint j = cell_mapping.MapFaceNode(f,fj);
          
              double entry_aij = 0.0;
              for (size_t qp : qp_face_data.QuadraturePointIndices())
                entry_aij +=  qp_face_data.ShapeValue(i, qp) *qp_face_data.ShapeValue(j, qp)
                                * qp_face_data.JxW(qp);
              for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
                Acell[g][i][j] += aval / bval * entry_aij;
            }//for fj
          }//end true Robin
        }//for fi
      }//if Robin
	  
	    // Dirichlet boundary
	    if (bndry.type == BoundaryType::Dirichlet)
      {
 		    const size_t num_face_nodes = face.vertex_ids.size();

        const auto& boundary_value = bndry.mg_values[0];

        // loop over nodes of that face
        for (size_t fi=0; fi<num_face_nodes; ++fi)
        {
          const uint i = cell_mapping.MapFaceNode(f,fi);
          dirichlet_count[i] += 1;
          dirichlet_value[i] += boundary_value;
        }//for fi
      }//if Dirichlet
	  
    }//for face f
 
    //======================= Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); //node-mapping
    for (size_t i=0; i<num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);
 
    //======================= Assembly into system
    for (size_t i=0; i<num_nodes; ++i)
    {
      if (dirichlet_count[i]>0) //if Dirichlet boundary node
      {
        for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
          MatSetValue(A[g], imap[i], imap[i], 1.0, ADD_VALUES);
		    // because we use CFEM, a given node is common to several faces
		    const double aux = dirichlet_value[i]/dirichlet_count[i];
        for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
          VecSetValue(bext[g], imap[i], aux, ADD_VALUES);
      }
      else
      {
        for (size_t j=0; j<num_nodes; ++j)
        {
          if (dirichlet_count[j]==0) // not related to a dirichlet node
            for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
              MatSetValue(A[g], imap[i], imap[j], Acell[g][i][j], ADD_VALUES);
		      else // related to a dirichlet node
          {
		        const double aux = dirichlet_value[j]/dirichlet_count[j];
            for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
              rhs_cell[g][i] -= Acell[g][i][j]*aux;
          }
        }//for j
        for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g)
          VecSetValue(bext[g], imap[i], rhs_cell[g][i], ADD_VALUES);
      }
    }//for i
  }//for cell
 
  chi::log.Log() << "Global assembly";

  for (uint g=0; g<mg_diffusion::Solver::num_groups; ++g) {
    MatAssemblyBegin(A[g], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A[g], MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(bext[g]);
    VecAssemblyEnd(bext[g]);
  }
 
  chi::log.Log() << "Done global assembly";

}