# CMake file for building bethe-faster
#------------------------------------------------------------------------------
# Setup
cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
cmake_policy(VERSION 3.5.1)

# - project name
project(bethe-faster LANGUAGES CXX)

# - build type
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Set build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")

option(BUILD_TESTS "Whether to build tests" OFF)
#------------------------------------------------------------------------------
# Compiler flags

# - set C++11 flag
if (NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "-std=c++11")
endif()

if("${CMAKE_BUILD_TYPE} " STREQUAL "Debug ")
    set(CMAKE_CXX_FLAGS "-O0 -ggdb ${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_BUILD_TYPE} " STREQUAL "Release ")
    set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS}")
endif()

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Build products

# - build the sources into an object library
set(BETHE_SRCS src/bethe.cc)

set(BETHE_OBJECTS bethe-obj)
add_library(${BETHE_OBJECTS} OBJECT ${BETHE_SRCS})
set_property(TARGET ${BETHE_OBJECTS} PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(${BETHE_OBJECTS} PRIVATE include ${CMAKE_CURRENT_BINARY_DIR})

# - write the static library
set(BETHE_STATIC_LIBRARY_NAME bethe-faster-static)
add_library(${BETHE_STATIC_LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${BETHE_OBJECTS}>)
set_target_properties(${BETHE_STATIC_LIBRARY_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${BETHE_STATIC_LIBRARY_NAME} PRIVATE include ${CMAKE_CURRENT_BINARY_DIR})

# - write the shared library
set(BETHE_SHARED_LIBRARY_NAME bethe-faster-shared)
add_library(${BETHE_SHARED_LIBRARY_NAME} SHARED $<TARGET_OBJECTS:${BETHE_OBJECTS}>)
set_target_properties(${BETHE_SHARED_LIBRARY_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${BETHE_SHARED_LIBRARY_NAME} PRIVATE include ${CMAKE_CURRENT_BINARY_DIR})
include_directories(include)

if(BUILD_TESTS)
    add_subdirectory(test)
endif()

#------------------------------------------------------------------------------
# Install products

# - directories
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX} DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY DESTINATION include DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY DESTINATION lib DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# - libraries
install(TARGETS ${BETHE_STATIC_LIBRARY_NAME} ${BETHE_SHARED_LIBRARY_NAME} DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(DIRECTORY include/bethe-faster/ DESTINATION include/bethe-faster FILES_MATCHING PATTERN "*.h" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
