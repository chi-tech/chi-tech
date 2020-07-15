message (STATUS "Loading Downstream.cmake")

set(CHI_TECH_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#------------------------------------------------ DEPENDENCIES
if (NOT (DEFINED ENV{PETSC_ROOT}))
    message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
else()
    set(PETSC_ROOT    "$ENV{PETSC_ROOT}")
    message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")
endif()

if (NOT (DEFINED ENV{LUA_ROOT}))
    message(FATAL_ERROR "***** LUA_ROOT is not set *****")
else()
    set(LUA_ROOT    "$ENV{LUA_ROOT}")
    message(STATUS "LUA_ROOT set to ${LUA_ROOT}")
endif()

find_package(MPI)

include(GNUInstallDirs)
include(${CHI_TECH_DIR}/CHI_RESOURCES/Macros/Filter.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

#================================================ Include directories
include_directories("${LUA_ROOT}/include")
include_directories("${PETSC_ROOT}/include")

include_directories("${CHI_TECH_DIR}/CHI_TECH")
include_directories("${CHI_TECH_DIR}/CHI_TECH/ChiLua")
include_directories("${CHI_TECH_DIR}/CHI_TECH/ChiMPI")
include_directories("${CHI_TECH_DIR}/CHI_TECH/ChiLog")
include_directories("${CHI_TECH_DIR}/CHI_RESOURCES")
include_directories("${CHI_TECH_DIR}/Modules")
include_directories("${CHI_TECH_DIR}/CHI_TECH/ChiMath/SpatialDiscretization")


include_directories(SYSTEM ${MPI_INCLUDE_PATH})


#================================================ Library directories
link_directories("${LUA_ROOT}/lib")
link_directories("${PETSC_ROOT}/lib")
link_directories("${CHI_TECH_DIR}/chi_build")


# --------------------------- VTK
find_package(VTK COMPONENTS
        vtkCommonCore vtkCommonDataModel
        vtkIOLegacy vtkIOCore
        vtkIOXML vtkParallelCore vtkIOParallelXML
        vtkFiltersCore
        vtkIOEnSight
        REQUIRED)
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

set(CHI_LIBS lua m dl ${MPI_CXX_LIBRARIES} petsc ${VTK_LIBRARIES} ChiLib)
