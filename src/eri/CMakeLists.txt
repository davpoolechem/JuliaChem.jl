##################
#= project info =#
##################
project(SIMINT.JL)

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
set(OPT "-Ofast -funroll-loops") 
if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
  set(OPT " -g -O0")
endif(CMAKE_BUILD_TYPE STREQUAL DEBUG)

set(CMAKE_C_FLAGS "${OPT} -std=gnu99")

#set(SIMINT_PATH /home/david/Local/Documents/Programs/simint-install)
set(SIMINT_PATH /export/home/david/programs/install/simint)

#################################
#= set up actual build process =#
#################################
include_directories(${SIMINT_PATH}/include)
link_directories(${SIMINT_PATH}/lib64)
add_library(jeri SHARED simint.c)
target_link_libraries(jeri simint)

#install(TARGETS simint DESTINATION /home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri)
install(TARGETS jeri DESTINATION /export/home/david/projects/Julia/JuliaChem.jl/src/eri)
