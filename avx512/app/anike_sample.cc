#include "biometric_anike.hh"
#include "baseline_biometric_anike.hh"
#include "naive_biometric_anike.hh"

#ifdef BASELINE_MKHSS
using namespace baseline_anike;
#define mkhss baseline_mkhss
#else
using namespace anike;
#endif

// using namespace naive_anike;

// p = 1846676733488631055471255842241675918524634969943310544305024546510335111754942467855175334391304886738089012161892806684010679369213995855689802335628273984374651322457100666803755783382146601342956837891522725255623204154321951439681667011710360653069057516378690718722316644878184359035522184553662331270365968446822008203609426500960578829207809935403657564047100423746644192595412017524710813488670533965603293464156021031202360140589859051220120924452326967
// q = 1725876987217986047624500396971859681137265214517892022764955133522145650537095430314045088828089319223467973204866689257287084472431639057265828454653668217896505768149095595788843000748135974986150770180016555436749641844269146519701313369089352948184435258565377439750312707561252924698656404987527139949369182619563623089969086836486319123231640554043654596504290380661637096995077501815992908334347140447815056438846057257331379512924564768728797003234542679
mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");

int hamming_distance(const std::vector<bool> &x, const std::vector<bool> &y) {
    assert(x.size() == y.size());
    int distance = 0;
    for (size_t i = 0; i < x.size(); i++) {
        if (x[i] != y[i]) {
            distance++;
        }
    }
    return distance;
}

void output_benchmark_result(const crs &crs, double keygenA_time, double keygenB_time, double keyderA_time, double keyderB_time) {
    // Output in machine-readable format
    print_parameters(std::cout, crs.params);
    mkhss::print_parameters(crs.mkhss_crs.params, std::cout);
    #ifdef BASELINE_MKHSS
    std::cout << "mkhss_scheme: baseline" << std::endl;
    #else
    std::cout << "mkhss_scheme: optimized" << std::endl;
    #endif

    std::cout << "length: " << crs.params.length << std::endl;
    std::cout << "keygenA_time: " << keygenA_time << std::endl;
    std::cout << "keygenB_time: " << keygenB_time << std::endl;
    std::cout << "keyderA_time: " << keyderA_time << std::endl;
    std::cout << "keyderB_time: " << keyderB_time << std::endl;
}

int main(int argc, char **argv) {
    // Test the ANIKE parameters
    int length = argc == 1 ? 2048 : atoi(argv[1]);
    int threshold = length / 3; // DEBUG
    crs crs;
    parameters params;
    public_encoding pe_A, pe_B;
    private_state st_A, st_B;
    create_parameters(params, 128, 128, 3072, length);
    print_parameters(std::cerr, params);
    std::cerr << "setup..." << std::flush;
    setup(crs, params, test_n.get_mpz_t());
    std::cerr << " done" << std::endl;
    mkhss::print_parameters(crs.mkhss_crs.params, std::cerr);

    std::vector<bool> x(params.length), y(params.length);
    // Fill with random bits
    srand(time(nullptr));
    for (int i = 0; i < params.length; i++) {
        x[i] = rand() % 2;
        y[i] = x[i] ^ (rand() <= (uint64_t) RAND_MAX * threshold / length);
    }
    std::cerr << "x: ";
    for (int i = 0; i < params.length; i++) {
        std::cerr << x[i];
    }
    std::cerr << std::endl;
    std::cerr << "y: ";
    for (int i = 0; i < params.length; i++) {
        std::cerr << y[i];
    }
    std::cerr << std::endl;

    int actual_distance = hamming_distance(x, y);
    std::cerr << "Actual Hamming distance: " << actual_distance << std::endl;
    // std::fill(x.begin(), x.end(), false);
    // x[0] = true;
    // // x[1] = true;
    // x[2] = true;
    // // x[4] = true;
    // std::fill(y.begin(), y.end(), false);
    // y[2] = true;
    // // y[3] = true;

    std::cerr << "attr_key_gen(sigma = A)... " << std::flush;
    double keygenA_time = measure_time_seconds([&]() {
        attr_key_gen(crs, x, pe_A, st_A);
    });
    std::cerr << keygenA_time << " seconds" << std::endl;
    std::cerr << "attr_key_gen(sigma = B)... " << std::flush;
    double keygenB_time = measure_time_seconds([&]() {
        attr_key_gen(crs, y, pe_B, st_B);
    });
    std::cerr << keygenB_time << " seconds" << std::endl;

    std::cout << "Public encoding size (Alice): " << get_public_encoding_size(crs, pe_A) << std::endl;
    std::cout << "Public encoding size (Bob): " << get_public_encoding_size(crs, pe_B) << std::endl;

    // FMPZ z_A, z_B;
    // std::cerr << "attr_key_der(sigma = A)" << std::endl;
    // attr_key_der_debug(crs, party_id::ALICE, st_A, pe_B, z_A);
    // std::cerr << "attr_key_der(sigma = B)" << std::endl;
    // attr_key_der_debug(crs, party_id::BOB, st_B, pe_A, z_B);
    // z_A.sub(z_A, z_B);
    // std::cerr << "z_A - z_B = " << z_A << std::endl;
    // return 0;

    std::cerr << "attr_key_der(sigma = A, threshold = " << threshold << " )... " << std::flush;
    std::vector<uint8_t[KEY_LENGTH]> shared_keys_A, shared_keys_B;
    double keyderA_time = measure_time_seconds([&]() {
        shared_keys_A = attr_key_der(crs, party_id::ALICE, st_A, pe_B, threshold);
    });
    std::cerr << keyderA_time << " seconds" << std::endl;
    std::cerr << "attr_key_der(sigma = B, threshold = " << threshold << " )... " << std::flush;
    double keyderB_time = measure_time_seconds([&]() {
        shared_keys_B = attr_key_der(crs, party_id::BOB, st_B, pe_A, threshold);
    });
    std::cerr << keyderB_time << " seconds" << std::endl;

    // for (int i = 0; i < threshold; i++) {
    //     std::cerr << "shared_keys_A[" << i << "] = ";
    //     for (int j = 0; j < KEY_LENGTH; j++) {
    //         std::cerr << std::hex << (int)shared_keys_A[i][j];
    //     }
    //     std::cerr << std::dec << std::endl;
    //     std::cerr << std::endl;
    // }

    // for (int i = 0; i < threshold; i++) {
    //     std::cerr << "shared_keys_B[" << i << "] = ";
    //     for (int j = 0; j < KEY_LENGTH; j++) {
    //         std::cerr << std::hex << (int)shared_keys_B[i][j];
    //     }
    //     std::cerr << std::dec << std::endl;
    //     std::cerr << std::endl;
    // }

    // Find indices where the keys match
    std::vector<int> matching_indices;
    for (int i = 0; i < threshold; i++) {
        if (memcmp(shared_keys_A[i], shared_keys_B[i], KEY_LENGTH) == 0) {
            matching_indices.push_back(i);
        }
    }

    if (matching_indices.empty()) {
        std::cerr << "No matching keys found." << std::endl;
        assert(actual_distance >= threshold);
    } else {
        std::cerr << "Matching keys found at indices: ";
        for (int index : matching_indices) {
            std::cerr << index << " ";
        }
        std::cerr << std::endl;
        assert(matching_indices.size() == 1);
        assert(matching_indices[0] == actual_distance);
    }

    std::cerr << std::endl;

    // for (int threshold = 0; threshold < crs.params.length; threshold++) {
    //     std::cerr << "# threshold = " << threshold << std::endl;

    //     std::cerr << "attr_key_der(sigma = A)" << std::endl;
    //     std::vector<uint8_t[KEY_LENGTH]> shared_keys_A = attr_key_der(crs, party_id::ALICE, st_A, pe_B, threshold);
    //     std::cerr << "attr_key_der(sigma = B)" << std::endl;
    //     std::vector<uint8_t[KEY_LENGTH]> shared_keys_B = attr_key_der(crs, party_id::BOB, st_B, pe_A, threshold);
    
    //     for (int i = 0; i < threshold; i++) {
    //         std::cerr << "shared_keys_A[" << i << "] = ";
    //         for (int j = 0; j < KEY_LENGTH; j++) {
    //             std::cerr << std::hex << (int)shared_keys_A[i][j];
    //         }
    //         std::cerr << std::endl;
    //     }
    
    //     for (int i = 0; i < threshold; i++) {
    //         std::cerr << "shared_keys_B[" << i << "] = ";
    //         for (int j = 0; j < KEY_LENGTH; j++) {
    //             std::cerr << std::hex << (int)shared_keys_B[i][j];
    //         }
    //         std::cerr << std::endl;
    //     }
    
    //     // Find indices where the keys match
    //     std::vector<int> matching_indices;
    //     for (int i = 0; i < threshold; i++) {
    //         if (memcmp(shared_keys_A[i], shared_keys_B[i], KEY_LENGTH) == 0) {
    //             matching_indices.push_back(i);
    //         }
    //     }
    
    //     if (matching_indices.empty()) {
    //         std::cerr << "No matching keys found." << std::endl;
    //     } else {
    //         std::cerr << "Matching keys found at indices: ";
    //         for (int index : matching_indices) {
    //             std::cerr << index << " ";
    //         }
    //         std::cerr << std::endl;
    //     }

    //     std::cerr << std::endl;
    // }

    output_benchmark_result(crs, keygenA_time, keygenB_time, keyderA_time, keyderB_time);

    return 0;
}
