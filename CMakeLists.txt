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

#------------------------------------------------ DEPENDECIES
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

if (NOT (DEFINED ENV{VTK_ROOT}))
    message(FATAL_ERROR "***** VTK_ROOT is not set *****")
else()
    set(VTK_ROOT    "$ENV{VTK_ROOT}")
    message(STATUS "VTK_ROOT set to ${VTK_ROOT}")
endif()

if (NOT (DEFINED ENV{VTK_VERSION}))
    message(FATAL_ERROR "***** VTK_VERSION is not set *****")
else()
    set(VTK_VERSION    "$ENV{VTK_VERSION}")
    message(STATUS "VTK_VERSION set to ${VTK_VERSION}")
endif()


#set(VTK_ROOT "/usr/local/Cellar/vtk/8.2.0")
#set(VTK_VERSION "-8.2")


option(USE_GRAPHICS "USE_GRAPHICS" OFF)
#set(CHI_USEGRAPHICS ON)
set(CHI_USESURFACEREMESHER OFF)  #pending fix for CHIVECTOR



if (APPLE)
    if (${USE_GRAPHICS})
        set(CHI_USEGRAPHICS ON)
    endif()
else()
  if (${USE_GRAPHICS})
   set(CHI_USEGRAPHICS ON)
  endif()
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
#add_definitions(-DLUA_STATIC)
add_definitions(-Wall -Wno-unused-value)
if (${CHI_USEGRAPHICS})
    add_definitions(-DCHI_USEGRAPHICS)
endif()
if (${CHI_USESURFACEREMESHER})
    add_definitions(-DCHI_USESURFACEREMESHER)
endif()


#================================================ Compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-format -Wno-unused-but-set-variable")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-write-strings")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-strict-aliasing")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -D_GLIBCXX_USE_CXX11_ABI=0")
add_definitions(-fpermissive)
if (${CHI_USEGRAPHICS})
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -fopenmp")
endif()

#================================================ Include directories
include_directories("${LUA_ROOT}/include")
include_directories("${CHI_TECH_DEP}/Random123/include")
include_directories("${TRIANGLE_ROOT}")
include_directories("${PROJECT_SOURCE_DIR}/CHI_TECH")
include_directories("${PROJECT_SOURCE_DIR}/CHI_TECH/CHI_MATH/CHI_DISCRETIZATION")
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories("${PETSC_ROOT}/include")
include_directories("${BOOST_ROOT}/include")

include_directories("${VTK_ROOT}/include/vtk${VTK_VERSION}")

#================================================ Library directories
link_directories("${LUA_ROOT}/lib")
link_directories("${TRIANGLE_ROOT}")
link_directories("${PETSC_ROOT}/lib")
link_directories("${BOOST_ROOT}/lib")

link_directories("${VTK_ROOT}/lib")
link_directories("${VTK_ROOT}/lib64")

set(TRIANGLE "${TRIANGLE_ROOT}/triangle.o")

set(LIBS lua m dl ${MPI_CXX_LIBRARIES} petsc)
set(LIBS ${LIBS} vtkCommonCore${VTK_VERSION} vtkCommonDataModel${VTK_VERSION})
set(LIBS ${LIBS} vtkIOLegacy${VTK_VERSION} vtkIOCore${VTK_VERSION})
set(LIBS ${LIBS} vtkIOXML${VTK_VERSION} vtkParallelCore${VTK_VERSION} )
set(LIBS ${LIBS} vtkIOParallelXML${VTK_VERSION})

if (${CHI_USEGRAPHICS})
    # OpenGL Lib
    find_package(OpenGL REQUIRED)
    include_directories(${OPENGL_INCLUDE_DIR})



    if (APPLE)
        include_directories("/opt/X11/include")
        link_directories("/opt/X11/lib")
        #set(GLEW_DIR "/usr/local")
        set(GLEW_DIR "${CHI_TECH_DEP}/glew-2.1.0")
        include_directories(${GLEW_DIR}/include)
        link_directories(${GLEW_DIR}/lib)
        add_definitions(-DGLEW_STATIC)
        set(LIBS ${LIBS} X11 GLEW GL GLU gomp )
        set (CMAKE_EXE_LINKER_FLAGS "-Wl,-undefined,dynamic_lookup,-no_pie -O3 -DNDEBUG")
    else()
        set(GLEW_DIR "${CHI_TECH_DEP}/GLEW")
        include_directories(${GLEW_DIR}/include)
        add_definitions(-DGLEW_STATIC)
        set(LIBS ${LIBS} X11 GLEW GL GLU gomp)
    endif()


endif()


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
add_subdirectory("${CHI_TECH_MOD}/CHI_NPTRANSPORT")



if (${CHI_USEGRAPHICS})
    add_subdirectory("${CHI_TECH_DIR}/CHI_GRAPHICS")
    add_subdirectory("${CHI_TECH_DIR}/CHI_WINDOWMANAGER")
    #add_subdirectory("${CHI_TECH_DIR}/CHI_SURFACE")
endif()

if(${CHI_USESURFACEREMESHER})
    add_subdirectory("${CHI_TECH_MOD}/CHI_SURFACEREMESHER")
endif()

add_executable(${TARGET} ${SOURCES})
target_link_libraries(${TARGET} ${LIBS} ${TRIANGLE})

# |------------ Write Makefile to root directory
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build -s\n\n"
        "clean:\n\t$(MAKE) -C chi_build -s clean\n")
