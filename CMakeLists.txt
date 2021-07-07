# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( mesh-2-pc.out )

set (CMAKE_CXX_STANDARD 11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# OpenMesh and its components
find_package( OpenMesh REQUIRED )

if ( NOT OpenMesh_FOUND )

  message(STATUS "This project requires the OpenMesh library, and will not be compiled.")

  return()  

endif()

# Multithreading
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package( Threads REQUIRED )

if ( NOT Threads_FOUND )

  message(STATUS "This project requires the Multithreading capacity, and will not be compiled.")

  return()  

endif()

# Creating entries for target: mesh-2-pc.out
# ############################

add_executable( mesh-2-pc.out main.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh-2-pc.out )

# Link the executable to CGAL and third-party libraries
target_link_libraries(mesh-2-pc.out PRIVATE CGAL::CGAL OpenMeshCore OpenMeshTools Threads::Threads)

