file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file (GLOB_RECURSE MORE_SOURCES2 "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

FOREACH(item ${MORE_SOURCES2})
    IF(${item} MATCHES "FieldFunction/fieldfunction_maindoc.cc")
        LIST(REMOVE_ITEM MORE_SOURCES2 ${item})
    ENDIF(${item} MATCHES "FieldFunction/fieldfunction_maindoc.cc")
ENDFOREACH(item)

FOREACH(item ${MORE_SOURCES2})
    IF(${item} MATCHES "lua/physics_maindoc_page.cc")
        LIST(REMOVE_ITEM MORE_SOURCES2 ${item})
    ENDIF(${item} MATCHES "lua/physics_maindoc_page.cc")
ENDFOREACH(item)

set(SOURCES ${SOURCES} ${MORE_SOURCES} ${MORE_SOURCES2} PARENT_SCOPE)
#set_source_files_properties(chi_console_00_constrdestr.cc PROPERTIES COMPILE_FLAGS -Wno-effc++)

