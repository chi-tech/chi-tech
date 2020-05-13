/**\page DevManUsingLib Using Chi-Tech as a library
 *
 * \section devman01_02 So you want to use Chi-Tech as a library
 *
 * \subsection devman01_02_01 Step 1 - Make a application directory and
 *              and empty source file
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

SET(TARGET test)
project(${TARGET} C CXX)

include(~/Desktop/ChiTech/chi-tech/CHI_RESOURCES/Macros/Downstream.cmake)

set(SOURCES "test.cc" )

add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})
 \endcode

The `SET(TARGET test)` line sets the name of the eventual executable
to be used. You can use any name other than `test`.
The `include` statement allows you to connect to all the resources connected
to Chi-Tech, including lua, PETSc, etc.
You will have to find the location where you compiled Chi-Tech in order to
properly specify the location of `Downstream.cmake`. Once this is properly
specified, the cmake-variable `CHI_LIBS` will be defined and the necessary
include-files will be usable.
The source-file names are specified via the line `set(SOURCES "test.cc")`.
Additional source files can be specified using:

\code
 set(SOURCES ${SOURCES} "another_file.cc")
\endcode

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


 * */