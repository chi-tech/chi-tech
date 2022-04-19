#ifndef CHITECH_LBS_ADJ_RESPONSE_FUNCTION_H
#define CHITECH_LBS_ADJ_RESPONSE_FUNCTION_H

#include "ChiMesh/Cell/cell.h"

namespace lbs_adjoint
{

struct ResponseFunctionDesignation
{
  const std::string              name;
  const chi_mesh::LogicalVolume& logical_volume;
  const std::string              lua_functional;

  explicit ResponseFunctionDesignation(std::string in_name,
                                       const chi_mesh::LogicalVolume& in_logical_volume,
                                       std::string  in_lua_function_name) :
    name(std::move(in_name)),
    logical_volume(in_logical_volume),
    lua_functional(std::move(in_lua_function_name))
  {}

  std::vector<double> GetMGResponse(const chi_mesh::Cell& cell,size_t num_groups) const;
};

}//namespace lbs_adjoint

#endif //CHITECH_LBS_ADJ_RESPONSE_FUNCTION_H
