if (${CHI_USEGRAPHICS})
    file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
else()
endif()


file (GLOB_RECURSE MORE_SOURCES2 "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
set(SOURCES ${SOURCES} ${MORE_SOURCES} ${MORE_SOURCES2} PARENT_SCOPE)

#set_source_files_properties(chi_console_00_constrdestr.cc PROPERTIES COMPILE_FLAGS -Wno-effc++)