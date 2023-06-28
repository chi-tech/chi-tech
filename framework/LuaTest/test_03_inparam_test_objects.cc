#include "test_03_inparam_test_objects.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_tests
{

// ##################################################################
RegisterChiObject(chi_unit_tests, TestObject);

chi::InputParameters TestObject::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("General test object");
  params.SetDocGroup("DocUnitTests");

  params.AddOptionalParameter("solver_type", "A", "The solver type.");
  params.AddRequiredParameter<std::string>(
    "coupled_field", "The text name of the coupled field.");
  params.AddRequiredParameterBlock(
    "sub_obj1", "A block of parameters for chi_unit_tests::TestSubObject");

  chi::ParameterBlock sub_obj2_param_block("sub_obj2");
  sub_obj2_param_block.AddParameter("num_groups", 99);
  params.AddOptionalParameterBlock(
    "sub_obj2",
    sub_obj2_param_block,
    "A block of parameters for chi_unit_tests::TestSubObject");

  params.AddOptionalParameter(
    "limiter_type", 1, "Type of limiter to use in the solver");
  params.MarkParamaterDeprecatedWarning("limiter_type");

  params.AddOptionalParameter("scheme", "Zorba", "What scheme to use");
  params.MarkParamaterDeprecatedError("scheme");

  params.AddRequiredParameter<bool>("format", "What output format to use");
  params.MarkParamaterDeprecatedError("format");

  params.AddOptionalParameter("use_my_stuff", false, "Yeah please do");
  params.MarkParamaterRenamed("use_my_stuff", "Renamed to \"use_zaks_stuff\".");

  params.AddRequiredParameter<bool>("use_ragusas_stuff", "If you want");
  params.MarkParamaterRenamed("use_ragusas_stuff",
                              "Renamed to \"use_complicated_stuff\".");

  params.AddOptionalParameter<int>(
    "groupset_num_subsets",
    1,
    "The number of subsets to apply to the set of groups in this set. This is "
    "useful for increasing pipeline size for parallel simulations");

  using namespace chi_data_types;
  params.ConstrainParameterRange("groupset_num_subsets",
                                 AllowableRangeLowLimit::New<int>(1));
  return params;
}

TestObject::TestObject(const chi::InputParameters& params)
  : solver_type_(params.GetParamValue<std::string>("solver_type")),
    sub_obj1_(MakeInpParamsForObj(TestSubObject, params.GetParam("sub_obj1"))),
    sub_obj2_(MakeInpParamsForObj(TestSubObject, params.GetParam("sub_obj2")))
{
  Chi::log.Log() << "TestObject created "
                 << "solver_type=" << solver_type_;
}

// ##################################################################
RegisterChiObject(chi_unit_tests, TestSubObject);

chi::InputParameters TestSubObject::GetInputParameters()
{
  chi::InputParameters params;
  params.SetGeneralDescription("General test sub-object");
  params.SetDocGroup("DocUnitTests");

  params.AddRequiredParameter<size_t>(
    "num_groups", "Number of groups to use in the simulation");

  return params;
}

TestSubObject::TestSubObject(const chi::InputParameters& params)
  : num_groups_(params.GetParamValue<size_t>("num_groups"))
{
  Chi::log.Log() << "TestSubObject created "
                 << "num_groups=" << num_groups_;
}

// ##################################################################
RegisterChiObject(chi_unit_tests, ChildTestObject);

chi::InputParameters ChildTestObject::GetInputParameters()
{
  chi::InputParameters params = TestObject::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "General test child-object inheriting option from parent");
  params.SetDocGroup("DocUnitTests");
  // clang-format on

  params.ChangeExistingParamToOptional("coupled_field", "Q");
  params.ChangeExistingParamToRequired<std::string>("solver_type");

  params.AddOptionalParameter(
    "num_sub_groups", 1, "Number of sub-groups to use in the simultion");

  return params;
}

ChildTestObject::ChildTestObject(const chi::InputParameters& params)
  : TestObject(params),
    num_sub_groups_(params.GetParamValue<int>("num_sub_groups"))
{
  Chi::log.Log() << "ChildTestObject created "
                 << "num_sub_groups=" << num_sub_groups_;
}

} // namespace chi_unit_tests