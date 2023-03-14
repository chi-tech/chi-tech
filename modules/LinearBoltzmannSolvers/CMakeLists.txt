add_subdirectory("A_LBSSolver")
add_subdirectory("B_DO_Solver")
add_subdirectory("Ca_DO_SteadyState")
add_subdirectory("Cc_DO_KEigenvalue")
add_subdirectory("D_DO_Transient")
add_subdirectory("Cd_DO_SteadyState_Adjoint")
add_subdirectory("D_DO_RZ_SteadyState")
add_subdirectory("B_MIP_SteadyState")


set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)