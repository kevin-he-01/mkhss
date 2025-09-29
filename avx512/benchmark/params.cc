#include "benchmark/params.hh"

mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");
mpz_class test_n_squared = test_n * test_n;

// Pregenerate moduli for optimal MKHSS performance
std::unordered_map<int, mpz_class> moduli_table;

void print_benchmark_parameters() {
    mkhss::parameters test_params;
    // Initialize parameters
    mkhss::create_parameters(test_params, SHORT_EXPONENT_ASSUMPTION, STAT_SEC_PARAM, SEC_PARAM, N_WIDTH, LOG2_B_EVALUATE_MKHSS);
    // Print parameters
    mkhss::print_parameters(test_params);
    
    // for (int i = 1; i <= LOG2_B_MAX; i += STRIDE) {
    //     mkhss::parameters params;
    //     mkhss::create_parameters(params, SHORT_EXPONENT_ASSUMPTION, STAT_SEC_PARAM, SEC_PARAM, N_WIDTH, i);
    //     std::cout << "if log2_b = " << i << ", exponent_bits: " << params.exponent_bits << std::endl;
    // }
    std::cout << std::endl;
}

void initialize_moduli_table() {
    if (!moduli_table.empty()) {
        // Moduli table already initialized
        return;
    }
    std::cerr << "Initializing moduli table..." << std::endl;
    std::vector<int> moduli_lengths = {3072, 3296, 3968, 3520, 4224};
    for (int n_width : moduli_lengths) {
        mpz_class p, q;
        std::cerr << "Generating moduli for n_width = " << n_width << std::endl;
        generate_safe_prime(p.get_mpz_t(), n_width / 2 + (n_width % 2));
        generate_safe_prime(q.get_mpz_t(), n_width / 2);
        moduli_table[n_width] = p * q;
    }
    std::cerr << "Moduli table initialized with " << moduli_table.size() << " entries." << std::endl;
}
