cmake_minimum_required(VERSION 3.2)

message(STATUS "Configuring ${CMAKE_BUILD_TYPE} build of ChiTech.")

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED PETSC_ROOT)
    if (NOT (DEFINED ENV{PETSC_ROOT}))
        message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
    else()
        set(PETSC_ROOT    "$ENV{PETSC_ROOT}")
    endif()
endif()
message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")

if (NOT DEFINED LUA_ROOT)
    if (NOT (DEFINED ENV{LUA_ROOT}))
        message(FATAL_ERROR "***** LUA_ROOT is not set *****")
    else()
        set(LUA_ROOT    "$ENV{LUA_ROOT}")
    endif()
endif()
message(STATUS "LUA_ROOT set to ${LUA_ROOT}")

#================================================ Set project targets
SET(TARGET ChiTech)
project(${TARGET} C CXX)

#================================================ Set variables
set(CHI_TECH_DIR "${PROJECT_SOURCE_DIR}/CHI_TECH")
set(CHI_TECH_MOD "Modules")

#================================================ Set cmake variables
find_package(MPI)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CHI_RESOURCES/Macros")

if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()


#================================================ Include macros
include(GNUInstallDirs)
include(Filter)

#================================================ Precompiler flags
#LUA CONSOLE
add_definitions(-DNO_FONTGENESIS)


#================================================ Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable -Wno-sign-compare")
  set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force -O3 in release builds (some OS's might downgrade it)
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})

#================================================ Include directories
include_directories("${LUA_ROOT}/include")
include_directories("${PETSC_ROOT}/include")
include_directories("${PROJECT_SOURCE_DIR}/")
include_directories("${PROJECT_SOURCE_DIR}/CHI_TECH")
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


#================================================ Library directories
link_directories("${LUA_ROOT}/lib")
link_directories("${PETSC_ROOT}/lib")
# --------------------------- VTK
find_package(VTK COMPONENTS
             vtkCommonCore vtkCommonDataModel
	         vtkIOLegacy vtkIOCore
	         vtkIOXML vtkParallelCore vtkIOParallelXML
             vtkFiltersCore
             vtkIOEnSight
             REQUIRED PATHS $ENV{VTK_DIR})
if (NOT VTK_FOUND)
  message(FATAL_ERROR "VTK not found: ${VTK_NOT_FOUND_MESSAGE}")
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

set(CHI_LIBS lua m dl ${MPI_CXX_LIBRARIES} petsc ${VTK_LIBRARIES})

#================================================ Default include directories
include_directories("${CHI_TECH_DIR}/ChiLua")
include_directories("${CHI_TECH_DIR}/ChiMPI")
include_directories("${CHI_TECH_DIR}/ChiLog")
include_directories("${PROJECT_SOURCE_DIR}/Modules")

#================================================ Define source directories
set(SOURCES "${CHI_TECH_DIR}/chi_runtime.cc" "${CHI_TECH_DIR}/LuaTest/lua_test.cc")
add_subdirectory("${CHI_TECH_DIR}/ChiConsole")
add_subdirectory("${CHI_TECH_DIR}/ChiLua")
add_subdirectory("${CHI_TECH_DIR}/ChiMath")
add_subdirectory("${CHI_TECH_DIR}/ChiPhysics")
add_subdirectory("${CHI_TECH_DIR}/ChiGraph")

add_subdirectory("${CHI_TECH_DIR}/ChiTimer")
add_subdirectory("${CHI_TECH_DIR}/ChiMesh")
add_subdirectory("${CHI_TECH_DIR}/ChiMPI")
add_subdirectory("${CHI_TECH_DIR}/ChiLog")

add_subdirectory("${CHI_TECH_MOD}")

add_library(ChiLib STATIC ${SOURCES})
add_executable(${TARGET} "${CHI_TECH_DIR}/chi_tech_main.cc")
target_link_libraries(${TARGET} ChiLib)
target_link_libraries(ChiLib ${CHI_LIBS})

# |------------ Write Makefile to root directory
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build \n\n"
        "clean:\n\t$(MAKE) -C chi_build clean\n")

