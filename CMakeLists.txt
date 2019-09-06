cmake_minimum_required(VERSION 3.2)

#set(CMAKE_VERBOSE_MAKEFILE ON)
find_program(GNU_CXX_COMPILER NAMES gcc)
message(STATUS ${GNU_COMPILER})

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#------------------------------------------------ DEPENDENCIES
if (NOT (DEFINED ENV{BOOST_ROOT}))
    message(FATAL_ERROR "***** BOOST_ROOT is not set *****")
else()
    set(BOOST_ROOT    "$ENV{BOOST_ROOT}")
    message(STATUS "BOOST_ROOT set to ${BOOST_ROOT}")
endif()

if (NOT (DEFINED ENV{PETSC_ROOT}))
    message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
else()
    set(PETSC_ROOT    "$ENV{PETSC_ROOT}")
    message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")
endif()


#================================================ Set project targets
SET( CMAKE_BUILD_TYPE Release ... FORCE )
SET(TARGET ChiTech)
project(${TARGET} C CXX)

#================================================ Set variables
set(CHI_TECH_DIR "${PROJECT_SOURCE_DIR}/CHI_TECH")
set(CHI_TECH_MOD "${PROJECT_SOURCE_DIR}/CHI_TECH/CHI_MODULES")
set(CHI_TECH_DEP "${PROJECT_SOURCE_DIR}/CHI_RESOURCES/Dependencies")

set(TRIANGLE_ROOT "${CHI_TECH_DEP}/triangle")
set(LUA_ROOT      "${CHI_TECH_DEP}/lua-5.3.5/install")


#================================================ Set cmake variables
find_package(MPI)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CHI_RESOURCES/Macros")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_CXX_COMPILER ${GNU_COMPILER})



#================================================ Include macros
include(GNUInstallDirs)
include(Filter)

#================================================ Precompiler flags
#LUA CONSOLE
add_definitions(-DNO_FONTGENESIS)


#================================================ Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable -std=c++11 -O3 -DNDEBUG")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-format -Wno-unused-but-set-variable")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-write-strings")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-strict-aliasing")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -D_GLIBCXX_USE_CXX11_ABI=0")
#elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

#================================================ Include directories
include_directories("${LUA_ROOT}/include")
include_directories("${CHI_TECH_DEP}/Random123/include")
include_directories("${TRIANGLE_ROOT}")
include_directories("${PROJECT_SOURCE_DIR}/CHI_TECH")
include_directories("${PROJECT_SOURCE_DIR}/CHI_TECH/CHI_MATH/CHI_DISCRETIZATION")
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories("${PETSC_ROOT}/include")
include_directories(SYSTEM "${BOOST_ROOT}/include")


#================================================ Library directories
link_directories("${LUA_ROOT}/lib")
link_directories("${TRIANGLE_ROOT}")
link_directories("${PETSC_ROOT}/lib")
#link_directories("${BOOST_ROOT}/lib")

set(TRIANGLE "${TRIANGLE_ROOT}/triangle.o")

# --------------------------- VTK
find_package(VTK COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
	vtkIOLegacy vtkIOCore 
	vtkIOXML vtkParallelCore
	vtkIOParallelXML REQUIRED)
if (NOT VTK_FOUND)
	message("VTK not found: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")

if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
else ()
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
		TARGETS ${TARGET}
    MODULES ${VTK_LIBRARIES}
    )
endif()

set(LIBS lua m dl ${MPI_CXX_LIBRARIES} petsc ${VTK_LIBRARIES})


#================================================ Default include directories
include_directories("${CHI_TECH_DIR}/CHI_LUA")
include_directories("${CHI_TECH_DIR}/CHI_MPI")
include_directories("${CHI_TECH_DIR}/CHI_LOG")

#================================================ Define source directories
set(SOURCES "${CHI_TECH_DIR}/chi_tech_main.cc")
add_subdirectory("${CHI_TECH_DIR}/CHI_CONSOLE")
add_subdirectory("${CHI_TECH_DIR}/CHI_LIB")
add_subdirectory("${CHI_TECH_DIR}/CHI_LUA")
add_subdirectory("${CHI_TECH_DIR}/CHI_MATH")
add_subdirectory("${CHI_TECH_DIR}/CHI_PHYSICS")
add_subdirectory("${CHI_TECH_DIR}/CHI_GRAPH")

add_subdirectory("${CHI_TECH_DIR}/CHI_TIMER")
add_subdirectory("${CHI_TECH_DIR}/CHI_TOOLS")
add_subdirectory("${CHI_TECH_DIR}/CHI_MESH")
add_subdirectory("${CHI_TECH_DIR}/CHI_MPI")
add_subdirectory("${CHI_TECH_DIR}/CHI_LOG")

add_subdirectory("${CHI_TECH_MOD}/CHI_MONTECARLON")
add_subdirectory("${CHI_TECH_MOD}/CHI_DIFFUSION")
add_subdirectory("${CHI_TECH_MOD}/LinearBoltzmanSolver")



add_executable(${TARGET} ${SOURCES})
target_link_libraries(${TARGET} ${LIBS} ${TRIANGLE})

# |------------ Write Makefile to root directory
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build -s\n\n"
        "clean:\n\t$(MAKE) -C chi_build -s clean\n")
