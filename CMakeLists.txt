# CMake file for building bethe-faster
#------------------------------------------------------------------------------
# Setup
cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
cmake_policy(VERSION 3.5.1)

# - project name
project(bethe-faster LANGUAGES CXX)

# - build type
set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Set build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")

option(BUILD_TESTS "Whether to build tests" OFF)
option(USE_OPENMP "Whether to build with OpenMP support" ON)

#------------------------------------------------------------------------------
# Compiler flags

# - set C++14 flag
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if("${CMAKE_BUILD_TYPE} " STREQUAL "Debug ")
    if ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "Clang ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -Wno-unknown-pragmas -O0 -ggdb")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "AppleClang ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -Wno-unknown-pragmas -O0 -ggdb")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "GNU ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wno-maybe-uninitialized -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=4 -Wswitch-default -Wundef -Werror -Wconversion -Wno-unknown-pragmas -O0 -ggdb -rdynamic")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "Intel ")
        message(FATAL_ERROR "Unsupported compiler: Intel")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "MSVC ")
        message(FATAL_ERROR "Unsupported compiler: MSVC ")
    endif()

elseif("${CMAKE_BUILD_TYPE} " STREQUAL "Release ")
    if ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wno-unknown-pragmas -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -O3")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "AppleClang ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wno-unknown-pragmas -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -O3")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "GNU ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wno-unknown-pragmas -Wno-maybe-uninitialized -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=4 -Wswitch-default -Wundef -Werror -Wconversion -O3")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel ")
        message(FATAL_ERROR "Unsupported compiler: Intel")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "MSVC ")
        message(FATAL_ERROR "Unsupported compiler: MSVC")
    endif()

elseif("${CMAKE_BUILD_TYPE} " STREQUAL "RelWithDebInfo ")
    if ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -Wno-unknown-pragmas -O3 -ggdb")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "AppleClang ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=4 -Wno-unknown-pragmas -O3 -ggdb")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "GNU ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wno-maybe-uninitialized -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=4 -Wswitch-default -Wundef -Werror -Wconversion -Wno-unknown-pragmas -O3 -ggdb -rdynamic")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel ")
        message(FATAL_ERROR "Unsupported compiler: Intel")

    elseif ("${CMAKE_CXX_COMPILER_ID} " STREQUAL "MSVC ")
        message(FATAL_ERROR "Unsupported compiler: MSVC")
    endif()
endif()

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Build products

find_package(ROOT 5.26.00 REQUIRED COMPONENTS MathMore)
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES})
include(${ROOT_USE_FILE})

if (USE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

# - build the sources into an object library
set(BETHE_SRCS libbethe-faster/Particle.cc libbethe-faster/DetectorHelper.cc libbethe-faster/Propagator.cc libbethe-faster/ParticleFilter.cc libbethe-faster/FilterHelper.cc libbethe-faster/PlotHelper.cc libbethe-faster/ParticleHelper.cc libbethe-faster/QuickPidAlgorithm.cc)

set(BETHE_OBJECTS bethe-obj)
add_library(${BETHE_OBJECTS} OBJECT ${BETHE_SRCS})
set_property(TARGET ${BETHE_OBJECTS} PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(${BETHE_OBJECTS} PRIVATE libbethe-faster ${CMAKE_CURRENT_BINARY_DIR})

# - write the static library
set(BETHE_STATIC_LIBRARY_NAME bethe-faster-static)
add_library(${BETHE_STATIC_LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${BETHE_OBJECTS}>)
set_target_properties(${BETHE_STATIC_LIBRARY_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${BETHE_STATIC_LIBRARY_NAME} PRIVATE libbethe-faster ${CMAKE_CURRENT_BINARY_DIR})

# - write the shared library
set(BETHE_SHARED_LIBRARY_NAME bethe-faster-shared)
add_library(${BETHE_SHARED_LIBRARY_NAME} SHARED $<TARGET_OBJECTS:${BETHE_OBJECTS}>)
set_target_properties(${BETHE_SHARED_LIBRARY_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${BETHE_SHARED_LIBRARY_NAME} PRIVATE libbethe-faster ${CMAKE_CURRENT_BINARY_DIR})

include_directories(libbethe-faster)

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
install(DIRECTORY libbethe-faster/ DESTINATION include/bethe-faster FILES_MATCHING PATTERN "*.h" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
