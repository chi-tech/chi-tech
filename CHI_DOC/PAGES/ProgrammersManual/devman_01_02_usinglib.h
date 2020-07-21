/**\page DevManUsingLib Using Chi-Tech as a library
 *
 * \section devman01_02 So you want to use Chi-Tech as a library
 *
 * \subsection devman01_02_01 Step 1 - Make a application directory and empty source file
 *
 * Wherever you like, create a folder that will contain your application.
 *
 * \code
 * mkdir TestApp
 * cd TestApp
 * \endcode
 *
 * Now create an empty source file
 *
 * \code
 * >test.cc
 * \endcode
 *
 * ## _
 * \subsection devman01_02_02 Step 2 - Create CMakeLists.txt file
 *
 * In the folder you just created create a text file called `CMakeLists.txt`
 *
 * \code
 * >CMakeLists.txt
 * \endcode
 *
 * ## _
 * \subsection devman01_02_03 Step 3 - Edit the cmake file to connect to Chi-Tech
 *
 * There are 3 custom values to be included in the `CMakeLists.txt` file.
 *  - The desired application/executable name, for now we will call that <B>test</B>.
 *  - The Chi-Tech "downstream" cmake include, `Downstream.cmake`.
 *  - The source file names. For now we will only have `test.cc`
 *
\code
cmake_minimum_required(VERSION 3.2)

set(TARGET test)
project(${TARGET} C CXX)
set(CHI_TECH_DIR "~/Desktop/ChiTech/chi-tech")
include("${CHI_TECH_DIR}/CHI_RESOURCES/Macros/Downstream.cmake")

file (GLOB_RECURSE SOURCES "*.cc")
add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})
\endcode

The `set(TARGET test)` line sets the name of the eventual executable
to be used. You can use any name other than `test`.
The `set(CHI_TECH_DIR "~/Desktop/ChiTech/chi-tech")` line sets the directory
of an installed chi-tech. Change this line to point to the appropriate
directory.
The `include` statement allows you to connect to all the resources connected
to Chi-Tech, including lua, PETSc, etc.
You will have to find the location where you compiled Chi-Tech in order to
properly specify the location of `Downstream.cmake`. Once this is properly
specified, the cmake-variable `CHI_LIBS` will be defined and the necessary
include-files will be usable.
 The line `file (GLOB_RECURSE SOURCES "*.cc")` adds all `.cc` files,
 in the current directory, to
 the list of sources to compile. To specify specific
 source-file names use `set(SOURCES "test.cc")`. Alternatively cmake
 functionality can be googled to determine how to add sub-directories.

## _

\subsection devman01_02_04 Step 4 - Add the basic code to `test.cc`

\code
#include "chi_runtime.h"

int main(int argc, char* argv[])
{
    ChiTech::Initialize(argc,argv);

    //We will add code here

    ChiTech::Finalize();
    return 0;
}
\endcode

This code is the minimum needed to have everything available in Chi-Tech.
The basic namespace access is provided via `#include "chi_runtime.h"`

MPI initialization and PETSc initialization is handled via the call to
`ChiTech::Initialize()`.

Finally all MPI and PETSc related items are destroyed via the call to
`ChiTech::Finalize()`.

## _

\subsection devman01_02_05 Step 5 - Compile the code
Well, we should probably state here that you should make sure ChiTech compiled.
This can be done by running the test cases. If all is well then proceed to make
build directory (```TestApp/build```) for this app.

\code
mkdir build
cd build
\endcode

Next execute cmake from this directory, giving it the directory-location of the
CMakeLists.txt created in \ref devman01_02_02 "step 2".

\code
cmake ../
\endcode

If all is well here then simply run make:

\code
make -j4
\endcode

##_

 * \section devman01_03 Adding lua-input to library-using apps
 * For this section we refer to \ref devman01_02 "the section above" where we
 * created a test application.
 *
 * \subsection devman01_03_01 Step 1 - Create folder to contain lua sources.
 * Purely for categorization, create a ```lua``` folder in the ```TestApp```
 * folder
 *
 * \code
 * cd ../TestApp
 * mkdir lua
 * cd lua
 * \endcode
 *
 * \subsection devman01_03_02 Step 2 - Create a source file
 * Let us pretend that we will have a status update message that we want to be
 * printed from a lua call. Create a file called ```testapp_statusmess.cc```
 * with the following contents:
 *
 * \code
 * #include"ChiLua/chi_lua.h"
 * #include "chi_log.h"
 *
 * int chiPrintStatus(lua_State *L)
 * {
 *     ChiLog& log = ChiLog::GetInstance();
 *
 *     log.Log(LOG_0) << "Hello from here";
 *
 *     return 0;
 * }
 * \endcode
 *
 * ## _
 *
 * \subsection devman01_03_03a Step 3 - Modify the CMakeLists.txt file to compile this code
 * In preparation for adding many more sources to the ```lua``` folder, edit the
 * CMakeLists.txt file as follows:
 *
\code
cmake_minimum_required(VERSION 3.2)

SET(TARGET test)
project(${TARGET} C CXX)

include(~/Desktop/ChiTech/chi-tech/CHI_RESOURCES/Macros/Downstream.cmake)
include(~/Desktop/ChiTech/chi-tech/CHI_RESOURCES/Macros/Filter.cmake)

set(SOURCES "test.cc" )
file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lua/*.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES})

add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})
\endcode

This will allow any ```.cc``` file in the ```lua``` folder to be compiled with
the project. The line

 \code
 file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lua/*.cc")
 \endcode

 recursively searches through the specified folder and pumps any filenames
 ending with .cc into the variable ```MORE_SOURCES```. This variable gets added
 to the manually specified ```SOURCES``` variable using the line

 \code
 set(SOURCES ${SOURCES} ${MORE_SOURCES})
 \endcode

 ## _
 * \subsection devman01_03_04a Step 4 - Modify test.cc to register the lua function
 * Modify your ```test.cc``` source file to be
\code
#include "chi_runtime.h"
#include "ChiConsole/chi_console.h"

int main(int argc, char* argv[])
{
    ChiTech::Initialize(argc,argv);

    ChiConsole& console = ChiConsole::GetInstance();

    auto L = console.consoleState;
    #include "ChiMacros/lua_register_macro.h"

    RegisterFunction(chiPrintStatus);

    ChiTech::Finalize();
    return 0;
}
\endcode

Notice here that we firstly included the header file for the ```ChiConsole```
object. This is the object that allows us to connect to the lua state.
 \code
 #include "ChiConsole/chi_console.h"
 \endcode

 Next we obtained a reference to the console through

 \code
 ChiConsole& console = ChiConsole::GetInstance();
 \endcode

 We then setup the variable ```L```, which is needed by the macro file
 ```ChiMacros/lua_register_macro.h```.

 \code
 auto L = console.consoleState;
 #include "ChiMacros/lua_register_macro.h"
 \endcode

 This macro file exposes numerous lua macros, one of which is the
 ```RegisterFunction(x)``` macro where we will supply the ```x``` as the
 function call we defined earlier:

 \code
 RegisterFunction(chiPrintStatus);
 \endcode
 * */