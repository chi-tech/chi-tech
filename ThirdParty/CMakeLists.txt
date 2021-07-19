add_subdirectory("paraview")

set(TARGET ThirdParty)
add_library(${TARGET} STATIC
    $<TARGET_OBJECTS:ThirdPartyParaView>)
