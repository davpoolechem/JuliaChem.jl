####################
#   initial info  ##
####################
cmake_minimum_required(VERSION 2.8.12)
project(JuliaChem) 

######################
#   set build type  ##
######################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message("-- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

##################
#   enable C++  ##
##################
enable_language(CXX)

#add_library(tutorial_compiler_flags INTERFACE)
#target_compile_features(tutorial_compiler_flags INTERFACE cxx_std_1z)

###################################
#   setup c++ compilation flags   #
###################################
#target_compile_options(cxx_compiler_flags INTERFACE
#  "$<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-Werror=pedantic;-pedantic-errors>>"
#  "$<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-Werror=pedantic;-pedantic-errors>>"
#)

if (CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_CXX_FLAGS "-std=c++1z -O3 -ffast-math")
else()
  set(CMAKE_CXX_FLAGS "-std=c++1z -O0 -g")
endif() 

#################################
#   setup library directories   # 
#################################
set(EIGEN_DIR ${EIGEN_INPUT})
set(BOOST_DIR ${BOOST_INPUT})
set(JLCXX_DIR ${JLCXX_INPUT})
set(LIBINT_DIR ${LIBINT_INPUT})

INCLUDE_DIRECTORIES( SYSTEM
  ${EIGEN_DIR} ${BOOST_DIR}/include 
  ${JLCXX_DIR}/include ${LIBINT_DIR}/include
)

LINK_DIRECTORIES( SYSTEM
  ${JLCXX_DIR}/lib ${LIBINT_DIR}/lib
)

#################################
#   build and install library   # 
#################################
add_library(jeri SHARED src/jeri.cpp)
target_link_libraries(jeri PUBLIC -lcxxwrap_julia -lint2)

install (TARGETS jeri DESTINATION ${PROJECT_SOURCE_DIR}) 
