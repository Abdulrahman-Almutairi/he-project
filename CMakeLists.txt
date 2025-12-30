cmake_minimum_required(VERSION 3.16)
project(PaillierHE LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(paillierhe
  src/bigint/bigint.cpp
  src/bigint/mul_schoolbook.cpp
  src/bigint/mod_arith.cpp
  src/bigint/montgomery.cpp
  src/bigint/prime.cpp
  src/bigint/inv.cpp
  src/bigint/div_exact.cpp
  src/bigint/gcd.cpp
  src/paillier/paillier.cpp
  src/bfv/poly.cpp
)
target_include_directories(paillierhe PUBLIC include)

add_executable(bench_all bench/bench_all.cpp)
target_link_libraries(bench_all PRIVATE paillierhe)
set_target_properties(bench_all PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bench)

add_executable(he_cli tools/he_cli.cpp)
target_link_libraries(he_cli PRIVATE paillierhe)
set_target_properties(he_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(test_all tests/test_all.cpp)
target_link_libraries(test_all PRIVATE paillierhe)
set_target_properties(test_all PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

add_executable(test_bfv_poly tests/test_bfv_poly.cpp)
target_link_libraries(test_bfv_poly PRIVATE paillierhe)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
  target_link_libraries(paillierhe PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(test_all PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(he_cli PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(bench_all PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found (build will be single-threaded)")
endif()