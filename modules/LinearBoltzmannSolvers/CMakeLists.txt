add_subdirectory("LBSSteadyState")
add_subdirectory("LBSKEigenvalue")
add_subdirectory("LBSTransient")
add_subdirectory("LBSSteadyStateAdjoint")
add_subdirectory("LBSCurvilinearSteadyState")


set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)