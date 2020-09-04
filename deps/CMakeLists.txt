##################
#= project info =#
##################
project(JERI)

######################
#= debug vs release =#
######################
cmake_minimum_required(VERSION 2.6)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_VERBOSE_MAKEFILE ON)

####################################
#= fortran compiler and flag info =#
####################################
enable_language(C)

###################
#= general flags =#
###################
set(OPT "-O2") 
if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
  set(OPT " -g -O0")
endif(CMAKE_BUILD_TYPE STREQUAL DEBUG)

set(CMAKE_C_FLAGS "${OPT} -std=c99 -fPIC")

#################################
#= set up actual build process =#
#################################
include_directories(${SIMINT_ROOT}/include)

if(EXISTS ${SIMINT_ROOT}/lib64/libsimint.a)
  set(SIMINT_LIB_PATH ${SIMINT_ROOT}/lib64)
else()
  set(SIMINT_LIB_PATH ${SIMINT_ROOT}/lib)
endif()

link_directories(${SIMINT_LIB_PATH})

add_library(simint SHARED src/simint.c)
target_link_libraries(simint ${SIMINT_LIB_PATH}/libsimint.a)

install(TARGETS simint DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}) 
