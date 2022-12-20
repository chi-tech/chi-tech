add_subdirectory("DiffusionSolver")
add_subdirectory("LinearBoltzmannSolver")
add_subdirectory("LBSCurvilinear")
add_subdirectory("LBKEigenvalueSolver")
add_subdirectory("LBSAdjointSolver")
add_subdirectory("CFEMDiffusion")
add_subdirectory("DFEMDiffusion")
add_subdirectory("LBSTransientSolver")

set(MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/lua/chi_modules_lua.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
