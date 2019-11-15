#include "diffusion_solver.h"

#include <ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h>
#include <ChiPhysics/PhysicsMaterial/property00_thermconductivity.h>
#include <ChiPhysics/PhysicsMaterial/property01_scalarvalue.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>
#include <ChiPhysics/chi_physics.h>
extern ChiPhysics chi_physics_handler;

#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/***/
void chi_diffusion::Solver::GetMaterialProperties(int mat_id,
                                                  double& diffCoeff,
                                                  double& sourceQ,
                                                  double& sigmaa)
{
  if (mat_id<0)
  {
    chi_log.Log(LOG_0ERROR)
      << "Cell encountered with no material id. ";
    exit(EXIT_FAILURE);
  }

  if (mat_id>=chi_physics_handler.material_stack.size())
  {
    chi_log.Log(LOG_0ERROR)
      << "Cell encountered with material id pointing to "
         "non-existing material.";
    exit(EXIT_FAILURE);
  }

  chi_physics::Material* material =
    chi_physics_handler.material_stack[mat_id];

  //====================================== Process material properties
  diffCoeff = 1.0;
  sourceQ   = 1.0;
  sigmaa    = 0.0;

  //We absolutely need the diffusion coefficient so process error
  if ((property_map_D < 0) || (property_map_D >= material->properties.size()))
  {
    chi_log.Log(LOG_0ERROR)
      << "Solver diffusion coefficient mapped to property index "
      << property_map_D << " is not a valid index for material \""
      << material->name <<"\" id " << mat_id;
    exit(EXIT_FAILURE);
  }

  //For now we can only support scalar values so lets check that
  if (dynamic_cast<chi_physics::ScalarValue*>
      (material->properties[property_map_D]) ||
      dynamic_cast<chi_physics::ThermalConductivity*>
      (material->properties[property_map_D]))
  {
    diffCoeff = material->properties[property_map_D]->GetScalarValue();
  }
  else
  {
    chi_log.Log(LOG_0ERROR)
      << "Solver diffusion coefficient mapped to property index "
      << property_map_D << " is not a valid property type"
      << " for material \""
      << material->name <<"\" id " << mat_id
      << ". Currently SCALAR_VALUE and THERMAL_CONDUCTIVITY are the "
      << "only supported types.";
    exit(EXIT_FAILURE);
  }


  if ((property_map_q < material->properties.size()) &&
      (property_map_q >= 0))
  {
    if (dynamic_cast<chi_physics::ScalarValue*>
        (material->properties[property_map_q]))
    {
      sourceQ = material->properties[property_map_q]->GetScalarValue();
    }
    else
    {
      chi_log.Log(LOG_0ERROR)
        << "Source value mapped to property index "
        << property_map_q << " is not a valid property type"
        << " for material \""
        << material->name <<"\" id " << mat_id
        << ". Currently SCALAR_VALUE is the "
        << "only supported type.";
      exit(EXIT_FAILURE);
    }
  }

  if ((property_map_sigma < 0) ||
      (property_map_sigma >= material->properties.size()))
  {
//    chi_log.Log(LOG_0WARNING)
//      << "Solver absorbtion coefficient mapped to property index "
//      << property_map_sigma << " is not a valid index for material \""
//      << material->name <<"\" id " << mat_id
//      << ". Property will be set to zero.";

  }
  else
  {
    sigmaa = material->properties[property_map_sigma]->GetScalarValue();
  }


}

//###################################################################
/**Gets material properties various sources.*/
void chi_diffusion::Solver::GetMaterialProperties(int mat_id,
                                                  int cell_glob_index,
                                                  int cell_dofs,
                                                  std::vector<double>& diffCoeff,
                                                  std::vector<double>& sourceQ,
                                                  std::vector<double>& sigmaa,
                                                  int group,
                                                  int moment)
{
  if (mat_id<0)
  {
    chi_log.Log(LOG_0ERROR)
      << "Cell encountered with no material id. ";
    exit(EXIT_FAILURE);
  }

  if (mat_id>=chi_physics_handler.material_stack.size())
  {
    chi_log.Log(LOG_0ERROR)
      << "Cell encountered with material id pointing to "
         "non-existing material.";
    exit(EXIT_FAILURE);
  }

  chi_physics::Material* material =
    chi_physics_handler.material_stack[mat_id];

  //====================================== Process material properties
  diffCoeff.resize(cell_dofs,1.0);
  sourceQ.resize(cell_dofs,0.0);
  sigmaa.resize(cell_dofs,0.0);

  //####################################################### REGULAR MATERIAL
  if (material_mode == DIFFUSION_MATERIALS_REGULAR)
  {
    //We absolutely need the diffusion coefficient so process error
    if ((property_map_D < 0) || (property_map_D >= material->properties.size()))
    {
      chi_log.Log(LOG_0ERROR)
        << "Solver diffusion coefficient mapped to property index "
        << property_map_D << " is not a valid index for material \""
        << material->name <<"\" id " << mat_id;
      exit(EXIT_FAILURE);
    }

    //For now we can only support scalar values so lets check that
    if (dynamic_cast<chi_physics::ScalarValue*>
        (material->properties[property_map_D]) ||
        dynamic_cast<chi_physics::ThermalConductivity*>
        (material->properties[property_map_D]))
    {
      diffCoeff.assign(cell_dofs,
                       material->properties[property_map_D]->GetScalarValue());
    }
    else
    {
      chi_log.Log(LOG_0ERROR)
        << "Solver diffusion coefficient mapped to property index "
        << property_map_D << " is not a valid property type"
        << " for material \""
        << material->name <<"\" id " << mat_id
        << ". Currently SCALAR_VALUE and THERMAL_CONDUCTIVITY are the "
        << "only supported types.";
      exit(EXIT_FAILURE);
    }


    if ((property_map_q < material->properties.size()) &&
        (property_map_q >= 0))
    {
      if (dynamic_cast<chi_physics::ScalarValue*>
          (material->properties[property_map_q]))
      {
        sourceQ.assign(cell_dofs,
                       material->properties[property_map_q]->GetScalarValue());
      }
      else
      {
        chi_log.Log(LOG_0ERROR)
          << "Source value mapped to property index "
          << property_map_q << " is not a valid property type"
          << " for material \""
          << material->name <<"\" id " << mat_id
          << ". Currently SCALAR_VALUE is the "
          << "only supported type.";
        exit(EXIT_FAILURE);
      }
    }

    if (not ((property_map_sigma < 0) ||
             (property_map_sigma >= material->properties.size())))
    {
      sigmaa.assign(cell_dofs,
                    material->properties[property_map_sigma]->GetScalarValue());
    }
  }//regular

  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        SCALAR       Q
  else if (material_mode == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p=0; p<material->properties.size(); p++)
    {
      if (dynamic_cast<chi_physics::TransportCrossSections*>
          (material->properties[p]))
      {
        auto xs = (chi_physics::TransportCrossSections*)material->properties[p];

        if (!xs->diffusion_initialized)
          xs->ComputeDiffusionParameters();

        diffCoeff.assign(cell_dofs,xs->diffg[group]);
        sigmaa.assign(cell_dofs,xs->sigma_rg[group]);
        transportxs_found = true;
      }
    }//for properties

    if (!transportxs_found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material encountered with no tranport xs"
           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((property_map_q < material->properties.size()) &&
        (property_map_q >= 0))
    {
      if (dynamic_cast<chi_physics::ScalarValue*>
          (material->properties[property_map_q]))
      {
        sourceQ.assign(cell_dofs,
                       material->properties[property_map_q]->GetScalarValue());
      }
      else
      {
        chi_log.Log(LOG_0ERROR)
          << "Source value mapped to property index "
          << property_map_q << " is not a valid property type"
          << " for material \""
          << material->name <<"\" id " << mat_id
          << ". Currently SCALAR_VALUE is the "
          << "only supported type.";
        exit(EXIT_FAILURE);
      }
    }
  }//transport xs TTR
  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        FIELDFUNC    Q
  else if (material_mode == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p=0; p<material->properties.size(); p++)
    {
      if (dynamic_cast<chi_physics::TransportCrossSections*>
          (material->properties[p]))
      {
        auto xs = (chi_physics::TransportCrossSections*)material->properties[p];

        if (!xs->diffusion_initialized)
          xs->ComputeDiffusionParameters();

        diffCoeff.assign(cell_dofs,xs->diffg[group]);
        sigmaa.assign(cell_dofs,xs->sigma_rg[group]);
        transportxs_found = true;
      }
    }//for properties

    if (!transportxs_found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material encountered with no tranport xs"
           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((q_field != nullptr) and (grid->IsCellLocal(cell_glob_index)))
    {
      std::vector<int> mapping;
      std::vector<int> pwld_nodes;
      std::vector<int> pwld_cells;

      for (int i=0; i<cell_dofs; i++)
      {
        pwld_nodes.push_back(i);
        pwld_cells.push_back(cell_glob_index);
      }

      chi_mesh::FieldFunctionInterpolation ffinterp;
      ffinterp.grid_view = grid;
      ffinterp.CreatePWLDMapping(q_field->num_grps,
                                 q_field->num_moms,
                                 group-gi,moment,
                                 pwld_nodes,pwld_cells,
                                 *q_field->local_cell_dof_array_address,
                                 &mapping);

      for (int i=0; i<cell_dofs; i++)
      {
//        sourceQ[i] = q_field->field_vector_local->operator[](mapping[i]);
        try {
          sourceQ[i] = q_field->field_vector_local->at(mapping[i]);
        }
        catch (const std::out_of_range& o)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Mapping error i=" << i
            << " mapping[i]=" << mapping[i]
            << " g=" << group << "(" << G << ")"
            << " ffsize=" << q_field->field_vector_local->size()
            << " dof_count=" << pwld_local_dof_count
            << " cell_loc=" << grid->cells[cell_glob_index]->partition_id;
          exit(EXIT_FAILURE);
        }

      }


    }
    else if (q_field == nullptr)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material source set to field function however"
           " the field is empty or not set.";
      exit(EXIT_FAILURE);
    }
  }//transport xs TTF
  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        FIELDFUNC    Q
  else if (material_mode == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JFULL)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p=0; p<material->properties.size(); p++)
    {
      if (dynamic_cast<chi_physics::TransportCrossSections*>
          (material->properties[p]))
      {
        auto xs = (chi_physics::TransportCrossSections*)material->properties[p];

        if (!xs->diffusion_initialized)
          xs->ComputeDiffusionParameters();

        diffCoeff.assign(cell_dofs,xs->D_jfull);
        sigmaa.assign(cell_dofs,xs->sigma_a_jfull);
        transportxs_found = true;
      }
    }//for properties

    if (!transportxs_found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material encountered with no tranport xs"
           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((q_field != nullptr) and (grid->IsCellLocal(cell_glob_index)))
    {
      std::vector<int> mapping;
      std::vector<int> pwld_nodes;
      std::vector<int> pwld_cells;

      for (int i=0; i<cell_dofs; i++)
      {
        pwld_nodes.push_back(i);
        pwld_cells.push_back(cell_glob_index);
      }

      chi_mesh::FieldFunctionInterpolation ffinterp;
      ffinterp.grid_view = grid;
      ffinterp.CreatePWLDMapping(q_field->num_grps,
                                 q_field->num_moms,
                                 0,0,
                                 pwld_nodes,pwld_cells,
                                 *q_field->local_cell_dof_array_address,
                                 &mapping);

      for (int i=0; i<cell_dofs; i++)
      {
//        sourceQ[i] = q_field->field_vector_local->operator[](mapping[i]);
        try {
          sourceQ[i] = q_field->field_vector_local->at(mapping[i]);
        }
        catch (const std::out_of_range& o)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Mapping error i=" << i
            << " mapping[i]=" << mapping[i]
            << " g=" << group << "(" << G << ")"
            << " ffsize=" << q_field->field_vector_local->size()
            << " dof_count=" << pwld_local_dof_count
            << " cell_loc=" << grid->cells[cell_glob_index]->partition_id;
          exit(EXIT_FAILURE);
        }

      }


    }
    else if (q_field == nullptr)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material source set to field function however"
           " the field is empty or not set.";
      exit(EXIT_FAILURE);
    }
  }//transport xs TTF
  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        FIELDFUNC    Q
  else if (material_mode == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTF_JPART)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p=0; p<material->properties.size(); p++)
    {
      if (dynamic_cast<chi_physics::TransportCrossSections*>
          (material->properties[p]))
      {
        auto xs = (chi_physics::TransportCrossSections*)material->properties[p];

        if (!xs->diffusion_initialized)
          xs->ComputeDiffusionParameters();

        diffCoeff.assign(cell_dofs,xs->D_jpart);
        sigmaa.assign(cell_dofs,xs->sigma_a_jpart);
        transportxs_found = true;
      }
    }//for properties

    if (!transportxs_found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material encountered with no tranport xs"
           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((q_field != nullptr) and (grid->IsCellLocal(cell_glob_index)))
    {
      std::vector<int> mapping;
      std::vector<int> pwld_nodes;
      std::vector<int> pwld_cells;

      for (int i=0; i<cell_dofs; i++)
      {
        pwld_nodes.push_back(i);
        pwld_cells.push_back(cell_glob_index);
      }

      chi_mesh::FieldFunctionInterpolation ffinterp;
      ffinterp.grid_view = grid;
      ffinterp.CreatePWLDMapping(q_field->num_grps,
                                 q_field->num_moms,
                                 0,0,
                                 pwld_nodes,pwld_cells,
                                 *q_field->local_cell_dof_array_address,
                                 &mapping);

      for (int i=0; i<cell_dofs; i++)
      {
//        sourceQ[i] = q_field->field_vector_local->operator[](mapping[i]);
        try {
          sourceQ[i] = q_field->field_vector_local->at(mapping[i]);
        }
        catch (const std::out_of_range& o)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Mapping error i=" << i
            << " mapping[i]=" << mapping[i]
            << " g=" << group << "(" << G << ")"
            << " ffsize=" << q_field->field_vector_local->size()
            << " dof_count=" << pwld_local_dof_count
            << " cell_loc=" << grid->cells[cell_glob_index]->partition_id;
          exit(EXIT_FAILURE);
        }

      }


    }
    else if (q_field == nullptr)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Diffusion Solver: Material source set to field function however"
           " the field is empty or not set.";
      exit(EXIT_FAILURE);
    }
  }//transport xs TTF
  else
  {
    chi_log.Log(LOG_0ERROR)
      << "Diffusion Solver: Invalid material mode.";
    exit(EXIT_FAILURE);
  }


}



//###################################################################
/**Attempts to apply Dirichlet boundary conditions on the ith row
 * of the matrix. If it did apply the dirichlet boundary condition
 * then it will return true which should logically be used to prevent
 * other matrix columns to be populated.
 *
 * \param ir The row which is to receive the possible application of
 *           the dirichlet boundary condition.
 * \param ir_boundary_type Pointer to a index that will receive the
 *                         boundary type.
 * \param iref The reference vertex index to be used for boundary type
 *             information.
 *
 * Note on ir. In a multigroup simulation ir would be the dof number plus
 * offset calculations (i.e. idof*G+g).
 *
 *
 * */
bool chi_diffusion::Solver::ApplyDirichletI(int ir,
                                            int *ir_boundary_type,
                                            int iref,
                                            bool suppress_assembly)
{
  int irefr = iref;
  if (irefr<0)
    irefr = ir;
  //================= Check if i is on boundary
  *ir_boundary_type = NO_BOUNDARY; //NO_BOUNDARY
  int ir_boundary_index= 0;
  if (nodal_boundary_numbers[irefr]<0)
  {
    ir_boundary_index = abs(nodal_boundary_numbers[irefr])-1;
    *ir_boundary_type = boundaries[ir_boundary_index]->type;
  }

  if (*ir_boundary_type == DIFFUSION_DIRICHLET)
  {
    double ir_mat_entry =1.0;
    if (!suppress_assembly)
      MatSetValue(Aref,ir,ir,ir_mat_entry,ADD_VALUES);
    auto dirich_bound =
      (chi_diffusion::BoundaryDirichlet*)boundaries[ir_boundary_index];
    double bvalue = dirich_bound->boundary_value;
    VecSetValue(bref,ir,bvalue,ADD_VALUES);
    VecSetValue(xref,ir,bvalue,INSERT_VALUES);


    return true; //Applied
  }
  return false; //Not applied
}

//###################################################################
/**Applies Dirichlet boundary conditions on the jth DOF. This
 * prevent non-symmetric matrix entries for dirichlet boundary
 * conditions.
 *
 * \param jr The column which is to receive an entry.
 * \param ir The row which is to receive an entry.
 * \param ir_boundary_type Pointer to a index that will receive the
 *                         boundary type.
 * \param jref The reference vertex index to be used for boundary
 *             type information.
 *
 *             */
bool chi_diffusion::Solver::ApplyDirichletJ(int jr,int ir,
                                            double jr_mat_entry,
                                            int *jr_boundary_type,int jref)
{
  int jrefr = jref;
  if (jrefr<0)
    jrefr = jr;
  //================= Check if j is on boundary
  *jr_boundary_type = NO_BOUNDARY; //NO_BOUNDARY
  int jr_boundary_index= 1;
  if (nodal_boundary_numbers[jrefr]<0)
  {
    jr_boundary_index = abs(nodal_boundary_numbers[jrefr])-1;
    *jr_boundary_type  = boundaries[jr_boundary_index]->type;
  }

  //================= If not dirichlet then
  if (*jr_boundary_type == DIFFUSION_DIRICHLET)
  {
    auto dirich_bound =
      (chi_diffusion::BoundaryDirichlet*)boundaries[jr_boundary_index];
    double bvalue = -1*jr_mat_entry*dirich_bound->boundary_value;
    VecSetValues(bref,1,&ir,&bvalue,ADD_VALUES);

    return true; //Applied
  }//if dirichlet
  return false; //Not applied
}
