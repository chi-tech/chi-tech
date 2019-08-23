/** \defgroup LuaFieldFunc Field Functions
 * \ingroup LuaPhysics
 *
 * The concept of a field function is a mature concept in Hydrodynamic
 * simulation codes such as StarCCM+ and OpenFOAM. So what is it?
 *
 * \image html "FieldFunction.png" width=500px
 *
 * In the most basic sense a field function is the values of a field as a
 * function of space and time. This can be as simple as cell average values
 * located at cell centroids. Such a simple case would be a field function
 * that is defined by a set of delta functions, therefore the field is only
 * defined at those specific values via its associated grid (mesh).
 * In practice however the cell average values
 * will have do be constant on an entire cell and any position within a given
 * cell should give this value. Hence the association with a discretization
 * method. Given a grid and a discretization method the value of a field
 * function can be obtained without the need to know any information of where
 * that value is stored or from which physics solver this value was populated.
 *
 * This concept will be developed further.
 * */