add_subdirectory("A_LBSSolver")
add_subdirectory("B_LBSSteadyState")
add_subdirectory("LBSKEigenvalue")
add_subdirectory("LBSTransient")
add_subdirectory("LBSSteadyStateAdjoint")
add_subdirectory("LBSCurvilinearSteadyState")
add_subdirectory("LBSMIPSteadyState")


set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)