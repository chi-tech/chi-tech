# |---------------------------------------------------------
# | Module Settings
# |---------------------------------------------------------
cmake_minimum_required(VERSION 3.2)
project(CHI_AUDIO)

file(GLOB_RECURSE CHI_AUDIO_SOURCE "${CHI_TECH_MOD}/CHI_AUDIO/*.cpp")
add_library(CHI_AUDIO STATIC ${CHI_AUDIO_SOURCE})
include_directories("${CHI_TECH_MOD}/CHI_AUDIO")
list(APPEND CHI_EXMOD_LIBS "${CHI_EXMOD_LIBS}" ${CHI_AUDIO_LIB})
MESSAGE( STATUS "Test")

message (STATUS "CHI_AUDIO loaded.")