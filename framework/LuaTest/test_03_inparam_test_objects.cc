#include "test_03_inparam_test_objects.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_tests
{

// ##################################################################
RegisterChiObject(chi_unit_tests, TestObject);

chi_objects::InputParameters TestObject::GetInputParameters()
{
  chi_objects::InputParameters params;

  params.AddOptionalParameter("solver_type", "A");
  params.AddRequiredParameter<std::string>("coupled_field");
  params.AddRequiredParameterBlock("sub_obj1");

  chi_objects::ParameterBlock sub_obj2_param_block("sub_obj2");
  sub_obj2_param_block.AddParameter("num_groups", 99);
  params.AddOptionalParameterBlock("sub_obj2", sub_obj2_param_block);

  params.AddOptionalParameter("limiter_type", 1);
  params.MarkParamaterDeprecatedWarning("limiter_type");

  params.AddOptionalParameter("scheme", "Zorba");
  params.MarkParamaterDeprecatedError("scheme");

  params.AddRequiredParameter<bool>("format");
  params.MarkParamaterDeprecatedError("format");

  params.AddOptionalParameter("use_my_stuff", false);
  params.MarkParamaterRenamed("use_my_stuff", "Renamed to \"use_zaks_stuff\".");

  params.AddRequiredParameter<bool>("use_ragusas_stuff");
  params.MarkParamaterRenamed("use_ragusas_stuff",
                              "Renamed to \"use_complicated_stuff\".");

  return params;
}

TestObject::TestObject(const chi_objects::InputParameters& params)
  : solver_type_(params.GetParamValue<std::string>("solver_type")),
    sub_obj1_(MakeInpParamsForObj(TestSubObject, params.GetParam("sub_obj1"))),
    sub_obj2_(MakeInpParamsForObj(TestSubObject, params.GetParam("sub_obj2")))
{
  chi::log.Log() << "TestObject created "
                 << "solver_type=" << solver_type_;
}

// ##################################################################
RegisterChiObject(chi_unit_tests, TestSubObject);

chi_objects::InputParameters TestSubObject::GetInputParameters()
{
  chi_objects::InputParameters params;

  params.AddRequiredParameter<size_t>("num_groups");

  return params;
}

TestSubObject::TestSubObject(const chi_objects::InputParameters& params)
  : num_groups_(params.GetParamValue<size_t>("num_groups"))
{
  chi::log.Log() << "TestSubObject created "
                 << "num_groups=" << num_groups_;
}

// ##################################################################
RegisterChiObject(chi_unit_tests, ChildTestObject);

chi_objects::InputParameters ChildTestObject::GetInputParameters()
{
  chi_objects::InputParameters params = TestObject::GetInputParameters();

  params.ChangeExistingParamToOptional("coupled_field", "Q");
  params.ChangeExistingParamToRequired<std::string>("solver_type");

  params.AddOptionalParameter("num_sub_groups", 1);

  return params;
}

ChildTestObject::ChildTestObject(const chi_objects::InputParameters& params)
  : TestObject(params),
    num_sub_groups_(params.GetParamValue<int>("num_sub_groups"))
{
  chi::log.Log() << "ChildTestObject created "
                 << "num_sub_groups=" << num_sub_groups_;
}

} // namespace chi_unit_tests