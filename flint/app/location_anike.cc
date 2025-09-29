#include "location_anike.hh"
#include <iomanip>
#include <nlohmann/json.hpp>

#ifdef BASELINE_MKHSS
#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss
#else
#include "mkhss.hh"
#endif

using namespace anike::geolocation;

// p = 1846676733488631055471255842241675918524634969943310544305024546510335111754942467855175334391304886738089012161892806684010679369213995855689802335628273984374651322457100666803755783382146601342956837891522725255623204154321951439681667011710360653069057516378690718722316644878184359035522184553662331270365968446822008203609426500960578829207809935403657564047100423746644192595412017524710813488670533965603293464156021031202360140589859051220120924452326967
// q = 1725876987217986047624500396971859681137265214517892022764955133522145650537095430314045088828089319223467973204866689257287084472431639057265828454653668217896505768149095595788843000748135974986150770180016555436749641844269146519701313369089352948184435258565377439750312707561252924698656404987527139949369182619563623089969086836486319123231640554043654596504290380661637096995077501815992908334347140447815056438846057257331379512924564768728797003234542679
mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");

struct location_anike_profile {
    double setup_time = 0.0; // Time taken for setup
    double keygenA_time, keygenB_time, keyderA_time, keyderB_time;
    keygen_profile keygenA;
    keygen_profile keygenB;
    keyder_profile keyderA;
    keyder_profile keyderB;

    inline nlohmann::json to_json() const {
        nlohmann::json j;
        j["setup_time"] = setup_time;
        j["keygenA_time"] = keygenA_time;
        j["keygenB_time"] = keygenB_time;
        j["keyderA_time"] = keyderA_time;
        j["keyderB_time"] = keyderB_time;
        j["keygenA"] = keygenA.to_json();
        j["keygenB"] = keygenB.to_json();
        j["keyderA"] = keyderA.to_json();
        j["keyderB"] = keyderB.to_json();
        return j;
    }
};

void output_benchmark_result(const crs &crs, const location_anike_profile &profs) {
    // Output in machine-readable format
    mkhss::print_parameters(crs.mkhss_crs.params, std::cout);
    print_parameters(std::cout, crs.params);
    // #ifdef BASELINE_MKHSS
    // std::cout << "mkhss_scheme: baseline" << std::endl;
    // #else
    // std::cout << "mkhss_scheme: optimized" << std::endl;
    // #endif

    // std::cout << "anike_app: geolocation" << std::endl;
    // std::cout << "L: " << crs.params.l << std::endl;
    // std::cout << "D: " << crs.params.D << std::endl;
    // std::cout << std::endl;
    nlohmann::json j;
    j["mkhss_scheme"] = mkhss::IS_BASELINE ? "baseline" : "optimized";
    j["anike_app"] = "geolocation";
    j["L"] = crs.params.l;
    j["D"] = crs.params.D;
    j["profile"] = profs.to_json();
    std::cout << j.dump(4) << std::endl;
}

int main(int argc, char **argv) {
    srand(time(nullptr)); // Seed the random number generator
    crs crs;
    parameters params;
    public_encoding pe_A, pe_B;
    private_state st_A, st_B;
    assert(argc == 3 || argc == 4);
    int coord_length = atoi(argv[1]);
    int dimensions = atoi(argv[2]);
    bool match = true;
    if (argc == 4) {
        assert(std::string(argv[3]) == "match" || std::string(argv[3]) == "nomatch");
        match = std::string(argv[3]) == "match";
    }
    create_parameters(params, 128, 128, 3072, coord_length, dimensions);
    std::cerr << "setup..." << std::flush;
    double setup_time = measure_time_seconds([&]() {
        setup(crs, params, test_n.get_mpz_t());
    });
    std::cerr << " done" << std::endl;
    mkhss::print_parameters(crs.mkhss_crs.params, std::cerr);
    print_parameters(std::cerr, params);

    // coord_t u_A, v_A, u_B, v_B, d;
    // u_A = 100; // Example coordinates for Alice
    // v_A = 120; // Example coordinates for Alice
    // u_B = 110; // Example coordinates for Bob
    // v_B = 110; // Example coordinates for Bob
    // d = 11; // Example distance

    std::vector<coord_t> x_A, x_B, d;
    for (int i = 0; i < dimensions; i++) {
        FMPZ random, randdist, actual_dist;
        FMPZ max, maxstrict, tmp, tmp2;
        fmpz_one_2exp(max.get_fmpz(), crs.params.l); // Maximum value for coordinates
        // fmpz_one_2exp(maxstrict.get_fmpz(), crs.params.l); // Maximum value for coordinates
        max.sub(max, 1); // Ensure max is less than 2^l
        // maxstrict.sub(maxstrict, 1); // Ensure maxstrict is less than 2^(l-1)
        // std::cerr << max;
        do {
            gen_random_fmpz(random, crs.params.l - 1);
            // std::cerr << "Random coordinate: " << random.get_ui() << std::endl;
            assert(random.cmp(max) <= 0 && "Random coordinate exceeds maximum value");
            assert(random.cmp(0) >= 0 && "Random coordinate is non-positive");
            gen_random_fmpz(randdist, crs.params.l - 1); // Generate a random distance
            // std::cerr << "Random distance: " << randdist.get_ui() << std::endl;
            tmp.add(random, randdist);
            tmp2.sub(random, randdist);
        } while (tmp2.cmp(0) < 0 || tmp.cmp(max) > 0 || randdist.cmp(0) == 0);
        
        gen_random_under_fmpz(actual_dist, randdist); // Ensure d is less than x_A
        if (rand() % 2) {
            fmpz_neg(actual_dist.get_fmpz(), actual_dist.get_fmpz()); // Randomly negate the distance to allow for both positive and negative distances
        }
        x_B.push_back(random.get_ui());
        FMPZ x_A_coord;
        x_A_coord.add(random, actual_dist); // Set x_B to be close to x_A
        // if (x_B_coord.cmp(max) >= 0) {
        //     // If x_B exceeds the maximum value, adjust it to stay within bounds
        //     x_B_coord.set(max);
        // }
        // if (x_B_coord.cmp(0) <= 0) {
        //     x_B_coord.set(0);
        // }
        x_A.push_back(x_A_coord.get_ui()); // Set x_B to be close to x_A
        d.push_back(randdist.get_ui());
    }
    if (!match) {
        int rand_dim = rand() % dimensions; // Randomly choose a dimension to ensure x_B is not close to x_A
        // Ensure x_B is not close to x_A
        x_A[rand_dim] = x_B[rand_dim] + d[rand_dim];
    }

    // Print coordinates for debugging
    std::cerr << "Coordinates for Alice: ";
    for (const auto &coord : x_A) {
        std::cerr << coord << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Coordinates for Bob: ";
    for (const auto &coord : x_B) {
        std::cerr << coord << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Distance: ";
    for (const auto &dist : d) {
        std::cerr << dist << " ";
    }
    std::cerr << std::endl;

    // bool is_coord_close = std::max(std::abs((intmax_t) u_A - u_B), std::abs((intmax_t) v_A - v_B)) < (intmax_t) d;
    // std::cerr << "Coordinates are " << (is_coord_close ? "close" : "not close") << std::endl;

    std::cerr << "attr_key_gen(sigma = A)... " << std::flush;
    keygen_profile keygenprof_A, keygenprof_B;
    double keygenA_time = measure_time_seconds([&]() {
        // attr_key_gen(crs, party_id::ALICE, u_A, v_A, std::nullopt, pe_A, st_A); // d only needs to be specified by one party, here it will be Bob
        attr_key_gen(crs, party_id::ALICE, x_A, {}, pe_A, st_A, keygenprof_A); // d only needs to be specified by one party, here it will be Bob
    });
    std::cerr << keygenA_time << " seconds" << std::endl;
    std::cerr << "attr_key_gen(sigma = B)... " << std::flush;
    double keygenB_time = measure_time_seconds([&]() {
        // attr_key_gen(crs, party_id::BOB, u_B, v_B, d, pe_B, st_B); // d only needs to be specified by one party, here it will be Bob
        attr_key_gen(crs, party_id::BOB, x_B, d, pe_B, st_B, keygenprof_B); // d only needs to be specified by one party, here it will be Bob
    });
    std::cerr << keygenB_time << " seconds" << std::endl;

    // std::cerr << "Public encoding size (Alice): " << get_public_encoding_size(crs, pe_A) << std::endl;
    // std::cerr << "Public encoding size (Bob): " << get_public_encoding_size(crs, pe_B) << std::endl;

    // FMPZ z_A, z_B;
    // std::cerr << "attr_key_der(sigma = A)... " << std::flush;
    // double keyderA_time = measure_time_seconds([&]() {
    //     attr_key_der_debug(crs, party_id::ALICE, st_A, pe_B, z_A);
    // });
    // std::cerr << keyderA_time << " seconds" << std::endl;
    // std::cerr << "attr_key_der(sigma = B)... " << std::flush;
    // double keyderB_time = measure_time_seconds([&]() {
    //     attr_key_der_debug(crs, party_id::BOB, st_B, pe_A, z_B);
    // });
    // std::cerr << keyderB_time << " seconds" << std::endl;
    // z_A.sub(z_A, z_B);
    // std::cerr << "z_A - z_B = " << z_A << std::endl;
    // return 0;

    std::cerr << "attr_key_der(sigma = A)... " << std::flush;
    uint8_t k_A[KEY_LENGTH], k_B[KEY_LENGTH];
    keyder_profile keyderprof_A, keyderprof_B;
    double keyderA_time = measure_time_seconds([&]() {
        attr_key_der(crs, party_id::ALICE, st_A, pe_B, k_A, keyderprof_A);
    });
    std::cerr << keyderA_time << " seconds" << std::endl;
    std::cerr << "attr_key_der(sigma = B)... " << std::flush;
    double keyderB_time = measure_time_seconds([&]() {
        attr_key_der(crs, party_id::BOB, st_B, pe_A, k_B, keyderprof_B);
    });
    std::cerr << keyderB_time << " seconds" << std::endl;

    std::cerr << std::endl;

    std::cerr << "Key A: ";
    for (size_t i = 0; i < KEY_LENGTH; i++) {
        std::cerr << std::hex << std::setfill('0') << std::setw(2) << (int) k_A[i];
    }
    std::cerr << std::dec << std::endl;
    std::cerr << "Key B: ";
    for (size_t i = 0; i < KEY_LENGTH; i++) {
        std::cerr << std::hex << std::setfill('0') << std::setw(2) << (int) k_B[i];
    }
    std::cerr << std::dec << std::endl;

    bool are_keys_equal = memcmp(k_A, k_B, KEY_LENGTH) == 0;
    std::cerr << "Keys are " << (are_keys_equal ? "equal" : "not equal") << std::endl;
    assert(are_keys_equal == match && "Keys should match if and only if coordinates are close");

    output_benchmark_result(crs, {setup_time, keygenA_time, keygenB_time, keyderA_time, keyderB_time, keygenprof_A, keygenprof_B, keyderprof_A, keyderprof_B});

    return 0;
}
