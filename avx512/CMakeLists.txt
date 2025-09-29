cmake_minimum_required(VERSION 3.10)

project(mkhss)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Libraries and packages
find_library(GMP_LIB gmp REQUIRED)
find_library(GMPXX_LIB gmpxx REQUIRED)
find_library(FLINT_LIB flint REQUIRED)
find_package(benchmark REQUIRED)
find_package(OpenSSL REQUIRED)

# Headers
include_directories(include)

# Sources
file(GLOB SOURCES "src/*.cc")
file(GLOB BENCHMARK_SRC "benchmark/*.cc")

# Compiler flags
add_compile_options(-march=native -mavx512f -mavx512ifma -g -O3)

# Executables
add_executable(location_anike app/location_anike.cc ${SOURCES})
add_executable(rns_location_anike app/location_anike.cc ${SOURCES})

add_executable(fuzzy_pake app/fuzzy_pake.cc ${SOURCES})
add_executable(rns_fuzzy_pake app/fuzzy_pake.cc ${SOURCES})

add_executable(test app/test.cc ${SOURCES})
add_executable(benchmark ${BENCHMARK_SRC} ${SOURCES})

# Target-specific compile options
target_compile_options(rns_fuzzy_pake PRIVATE -DUSE_RNS_MULTIPLICATION)
target_compile_options(rns_location_anike PRIVATE -DUSE_RNS_MULTIPLICATION)

# Link libraries
set(MKHSS_LIBS ${GMP_LIB} ${GMPXX_LIB} ${FLINT_LIB} OpenSSL::SSL OpenSSL::Crypto)
target_link_libraries(location_anike PRIVATE ${MKHSS_LIBS})
target_link_libraries(rns_location_anike PRIVATE ${MKHSS_LIBS})
target_link_libraries(fuzzy_pake PRIVATE ${MKHSS_LIBS})
target_link_libraries(rns_fuzzy_pake PRIVATE ${MKHSS_LIBS})
target_link_libraries(benchmark PRIVATE ${MKHSS_LIBS} benchmark::benchmark)
target_link_libraries(test PRIVATE ${MKHSS_LIBS})
