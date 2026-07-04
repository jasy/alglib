cmake_minimum_required(VERSION 3.14)
project(alglib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_COVERAGE "Enable coverage instrumentation" OFF)

if(ENABLE_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-g --coverage)
    add_link_options(--coverage)
  endif()
endif()

# Prefer building minimal GoogleTest (disable googlemock to save build time)
find_package(Threads REQUIRED)

# Disable building googlemock and samples/tests inside googletest to speed up builds.
# Projects embedding googletest can control these via -D options when configuring.
set(BUILD_GMOCK OFF CACHE BOOL "Do not build googlemock subproject" FORCE)
set(gtest_build_samples OFF CACHE BOOL "Do not build gtest samples" FORCE)
set(gtest_build_tests OFF CACHE BOOL "Do not build gtest internal tests" FORCE)

# Build the vendored googletest subproject (minimal targets only).
add_subdirectory(googletest)

add_executable(alglib sample.cpp)

target_include_directories(alglib PRIVATE ${CMAKE_SOURCE_DIR})

target_link_libraries(alglib PRIVATE gtest_main Threads::Threads)

enable_testing()
add_test(NAME alglib COMMAND alglib)
