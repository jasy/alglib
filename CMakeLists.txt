cmake_minimum_required(VERSION 3.14)
project(alglib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_COVERAGE "Enable coverage instrumentation" OFF)

include(ExternalProject)

set(GTEST_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/googletest")
set(GTEST_BUILD_DIR "${CMAKE_BINARY_DIR}/googletest/build")
set(GTEST_INSTALL_DIR "${CMAKE_BINARY_DIR}/googletest/install")

ExternalProject_Add(googletest
  SOURCE_DIR "${GTEST_SOURCE_DIR}"
  BINARY_DIR "${GTEST_BUILD_DIR}"
  INSTALL_DIR "${GTEST_INSTALL_DIR}"
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_GMOCK=OFF
    -DINSTALL_GTEST=ON
    -Dgtest_build_samples=OFF
    -Dgtest_build_tests=OFF
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_CXX_STANDARD=11
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    -DCMAKE_CXX_EXTENSIONS=OFF
  BUILD_BYPRODUCTS
    "${GTEST_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
    "${GTEST_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
  BUILD_COMMAND "${CMAKE_COMMAND}" --build "${GTEST_BUILD_DIR}" --config Release
  INSTALL_COMMAND "${CMAKE_COMMAND}" --build "${GTEST_BUILD_DIR}" --config Release --target install
  TEST_COMMAND ""
)

add_library(gtest STATIC IMPORTED)
set_target_properties(gtest PROPERTIES
  IMPORTED_LOCATION "${GTEST_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
add_dependencies(gtest googletest)

add_library(gtest_main STATIC IMPORTED)
set_target_properties(gtest_main PROPERTIES
  IMPORTED_LOCATION "${GTEST_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_LINK_LIBRARIES gtest
)
add_dependencies(gtest_main googletest)

# Build sample code for this algorithm library.
add_executable(alglib sample.cpp)
set_target_properties(alglib PROPERTIES COMPILE_WARNING_AS_ERROR ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(alglib PRIVATE -Wall -Wextra -Wpedantic -Wnon-virtual-dtor)
  if(ENABLE_COVERAGE)
    target_compile_options(alglib PRIVATE -g --coverage)
    target_link_libraries(alglib PRIVATE --coverage)
  endif()
endif()
target_include_directories(alglib PRIVATE "${GTEST_INSTALL_DIR}/include")
target_link_libraries(alglib PRIVATE gtest_main)
enable_testing()
include(GoogleTest)
gtest_discover_tests(alglib)

