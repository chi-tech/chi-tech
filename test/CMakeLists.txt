file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

set(TEST_SOURCES ${TEST_SOURCES} ${MORE_SOURCES} PARENT_SCOPE)
