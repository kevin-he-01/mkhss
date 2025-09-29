#include "benchmark/params.hh"

// Main function for Google Benchmark
int main(int argc, char** argv) {
    print_benchmark_parameters();
    
    // Initialize Google Benchmark
    benchmark::Initialize(&argc, argv);

    // Run the benchmarks
    // initialize_moduli_table();
    benchmark::RunSpecifiedBenchmarks();
    return 0;
}
