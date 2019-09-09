/**\page DevManAddingLua Adding Lua-input to the system

Lua is the main scripting system used by ChiTech.

\section devman03_sec0 Adding Lua-input

## Step 1: Create a source file
First create a `.cc` file in the location where you want the lua-input to
reside. Example: `CHI_TECH/lua_test.cc`
\n\n
If the file is not covered by a `CMakeLists.txt` recursion then follow
the instructions for \ref DevManAddingCode "adding a code to the system".

## Step 2: Define a standard lua template with your input function name
In our case the test function will be `LuaTest`. We add `chi` to the front of
all input language to distinguish it from other lua functions.
\n\n
Lua-wrapper functions have no real arguments other than the lua-state and
always return an integer representing the number of items returned. For more
information on how to handle parameters in lua C-api consult the
<a href="https://www.lua.org/pil/24.html">lua-documentation</a>.

\code
#include <ChiLua/chi_lua.h>

//###################################################################
///This is a lua test function.
int chiLuaTest(lua_State* L)
{
  return 0;
}
\endcode

## Step 3: Add your code to it
In this simple example we will simply take the first argument and print it
to `stdout`. We also include argument protection.

\code
#include <ChiLua/chi_lua.h>
#include <iostream>

//###################################################################
///This is a lua test function.
///\param argument1 Any Argument of any type.
int chiLuaTest(lua_State* L)
{
  //============================== Optional argument protection
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError("chiLuaTest",1,num_args);

  //============================== Obtain first argument from stack
  const char* argument_1 = lua_tostring(L,1);

  //============================== Print to screen
  std::cout << "LuaTest: " << argument_1 << std::endl;

  return 0;
}
\endcode

## Step 4: Add it to the lua-register
In the project directory go to the file `CHI_TECH/ChiLua/chi_lua_register.h`
and add a `RegisterFunction` call to the register, registering the lua function
call to be connected with the lua console.
\n\n
If you would like for the function call to be listed under a category then
proceed the `RegisterFunction` call with a comment starting exactly with
`//module:` which will create a module with text matching the text following
this piece.
\n
\n
The module tabulation in the manual and the documentation for the parameters
will be processed with a script therefore it is important to comment
parameters and return values. Some specific notes in the lua-register file

 - `//module:` provides a categorization of the function on the documenation
 main page.
 - `//\ref` will provide a link anchor on the documenation main page.

\code
 .
 .
 .
RegisterFunction(chiLBSGroupsetSetGroupSubsets)
RegisterFunction(chiLBSGroupsetSetIterativeMethod)
RegisterFunction(chiLBSGroupsetSetResidualTolerance)
RegisterFunction(chiLBSGroupsetSetMaxIterations)
RegisterFunction(chiLBSGroupsetSetGMRESRestartIntvl)
RegisterFunction(chiLBSGroupsetSetWGDSA)
RegisterFunction(chiLBSGroupsetSetTGDSA)

//module:Test scripts
RegisterFunction(chiLuaTest)


\endcode


## Final step: Test it!!

Create a script input file and call your function. If it is linked incorrectly
the lua scripting system will throw an error like

`Test.lua:2: attempt to call a nil value (global 'chiLuaTest')`
 */