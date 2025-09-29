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
add_compile_options(-g -O3)

# Executables
add_executable(sample app/sample.cc ${SOURCES})
add_executable(baseline_sample app/sample.cc ${SOURCES})
add_executable(location_anike app/location_anike.cc ${SOURCES})
add_executable(baseline_location_anike app/location_anike.cc ${SOURCES})
add_executable(fuzzy_pake app/fuzzy_pake.cc ${SOURCES})
add_executable(baseline_fuzzy_pake app/fuzzy_pake.cc ${SOURCES})
add_executable(test app/test.cc ${SOURCES})
add_executable(wheel app/wheel.cc)
add_executable(benchmark ${BENCHMARK_SRC} ${SOURCES})

# Target-specific compile options
target_compile_options(baseline_sample PRIVATE -DBASELINE_MKHSS)
target_compile_options(baseline_location_anike PRIVATE -DBASELINE_MKHSS)
target_compile_options(baseline_fuzzy_pake PRIVATE -DBASELINE_MKHSS)

# Link libraries
set(MKHSS_LIBS ${GMP_LIB} ${GMPXX_LIB} ${FLINT_LIB} OpenSSL::SSL OpenSSL::Crypto)
target_link_libraries(sample PRIVATE ${MKHSS_LIBS})
target_link_libraries(baseline_sample PRIVATE ${MKHSS_LIBS})
target_link_libraries(location_anike PRIVATE ${MKHSS_LIBS})
target_link_libraries(fuzzy_pake PRIVATE ${MKHSS_LIBS})
target_link_libraries(baseline_fuzzy_pake PRIVATE ${MKHSS_LIBS})
target_link_libraries(baseline_location_anike PRIVATE ${MKHSS_LIBS})
target_link_libraries(benchmark PRIVATE ${MKHSS_LIBS} benchmark::benchmark)
target_link_libraries(test PRIVATE ${MKHSS_LIBS})
