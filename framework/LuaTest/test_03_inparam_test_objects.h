#ifndef CHITECH_TEST_03_INPARAM_TEST_OBJECTS_H
#define CHITECH_TEST_03_INPARAM_TEST_OBJECTS_H

#include "ChiObject/chi_object.h"
#include "ChiParameters/input_parameters.h"

namespace chi_unit_tests
{

class TestSubObject : public chi_objects::ChiObject
{
private:
  const size_t num_groups_;
public:
  static chi_objects::InputParameters GetInputParameters();
  explicit TestSubObject(const chi_objects::InputParameters& params);
};

class TestObject : public chi_objects::ChiObject
{
private:
  const std::string solver_type_;
  TestSubObject sub_obj1_;
  TestSubObject sub_obj2_;
public:
  static chi_objects::InputParameters GetInputParameters();
  explicit TestObject(const chi_objects::InputParameters& params);
};

class ChildTestObject : public TestObject
{
private:
  const int num_sub_groups_;
public:
  static chi_objects::InputParameters GetInputParameters();
  explicit ChildTestObject(const chi_objects::InputParameters& params);
};



}

#endif // CHITECH_TEST_03_INPARAM_TEST_OBJECTS_H
