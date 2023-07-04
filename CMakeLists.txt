cmake_minimum_required(VERSION 2.8.12)

project(TrioBinning)

SET(MAX_KMER_SIZE "64" CACHE STRING "MAX_KMER_SIZE")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
        message("Debug mode")
        add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
        if(CMAKE_BUILD_TYPE MATCHES Profile)
                message("Profile mode")
                add_compile_options(-pg)
                set(CMAKE_SHARED_LINKER_FLAGS "-pg")
                set(CMAKE_EXE_LINKER_FLAGS "-pg")
        else(CMAKE_BUILD_TYPE MATCHES Profile)
                message("Release mode")
                add_compile_options(-O3)
        endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp *.hh *.tcc)

# Bifrost library
include(ExternalProject)

ExternalProject_Add(    
  bifrost_lib
  INSTALL_DIR ${PROJECT_SOURCE_DIR}/Bifrost/build
  BINARY_DIR ${PROJECT_SOURCE_DIR}/Bifrost/build
  TMP_DIR ${PROJECT_SOURCE_DIR}/Bifrost/build
  STAMP_DIR ${PROJECT_SOURCE_DIR}/Bifrost/build
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/Bifrost
  PREFIX ${PROJECT_SOURCE_DIR}/Bifrost
  CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DMAX_KMER_SIZE=${MAX_KMER_SIZE} -DCOMPILATION_ARCH=OFF DENABLE_AVX2=OFF
  UPDATE_COMMAND "" 
  INSTALL_COMMAND "")   

ExternalProject_Get_property(bifrost_lib INSTALL_DIR)

add_library(bifrost STATIC IMPORTED)    

include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/Bifrost/src")

set_property(TARGET bifrost PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/Bifrost/build/src/libbifrost.a")

add_dependencies(bifrost bifrost_lib)

# Pthread library
find_package(Threads REQUIRED)

# Zlib library
find_package(ZLIB REQUIRED)

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Zlib not found. Required for to output files")
endif(ZLIB_FOUND)

# Bin
add_executable(TrioBinning ${sources})

# Compilation
add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})
target_compile_options(TrioBinning PUBLIC "-DMAX_KMER_SIZE=${MAX_KMER_SIZE}")

# Linking
target_link_libraries(TrioBinning ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(TrioBinning ${ZLIB_LIBRARIES})
target_link_libraries(TrioBinning bifrost)

install(TARGETS TrioBinning DESTINATION bin)
