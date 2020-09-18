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
INCLUDE_DIRECTORIES( SYSTEM
  ${EIGEN_PATH} ${BOOST_PATH}/include 
  ${JLCXX_PATH}/include ${LIBINT_PATH}/include
)

LINK_DIRECTORIES( SYSTEM
  ${JLCXX_PATH}/lib ${LIBINT_PATH}/lib
)

message("-- Found Eigen at ${EIGEN_PATH}")
message("-- Found Boost at ${BOOST_PATH}")
message("-- Found JLCxx at ${JLCXX_PATH}")
message("-- Found Libint at ${LIBINT_PATH}")

#################################
#   build and install library   # 
#################################
add_library(jeri SHARED src/jeri.cpp)
target_link_libraries(jeri PUBLIC -lcxxwrap_julia -lint2)

install (TARGETS jeri DESTINATION ${PROJECT_SOURCE_DIR}) 
