#include "fuzzy_pake.hh"
#include <iomanip>

#ifdef BASELINE_MKHSS
#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss
#else
#include "mkhss.hh"
#endif

using namespace anike::fuzzy_pake;
using namespace std::string_literals;

// p = 1846676733488631055471255842241675918524634969943310544305024546510335111754942467855175334391304886738089012161892806684010679369213995855689802335628273984374651322457100666803755783382146601342956837891522725255623204154321951439681667011710360653069057516378690718722316644878184359035522184553662331270365968446822008203609426500960578829207809935403657564047100423746644192595412017524710813488670533965603293464156021031202360140589859051220120924452326967
// q = 1725876987217986047624500396971859681137265214517892022764955133522145650537095430314045088828089319223467973204866689257287084472431639057265828454653668217896505768149095595788843000748135974986150770180016555436749641844269146519701313369089352948184435258565377439750312707561252924698656404987527139949369182619563623089969086836486319123231640554043654596504290380661637096995077501815992908334347140447815056438846057257331379512924564768728797003234542679
mpz_class test_n("0x8c70e73ca01499a327a6733c4ebc9a6b23c1724eb2a2ba6023764ba3f7c09206f4936209cc83925c27dc15e1ee542050fe5f75d388c2e938ba76b7b56ac3131d09538dbd157e267b888d17d95053bff688744f76589d5e9cc8e48d7fb82079b2d6127e12ff00bbb89a7df722f950fb54c1f5856bfd5ac348af1847bdaeebd1398a2a34c2de9ff94a7101098714eae611eebed32f37dc1c33de4a0975b16c293973f37caa3368eeae049d0c70a15226198f07aecf051dd6e7ea303379726cf560cb444d0889714928501b83accbef555425218f28120c0758e97a857d18559556339161e185b343680052a3e7ef01e90b0c1f0339f06f01003958f3c8d5e3af0de618a6620bed48a489006a2042f206926aa12bd4b611fcf9bd6cf46c7cb9ad481caa7ac1bf6ff578967fd16bc4a0e47cef10fd3d4dc278922b6f4284959e3149f0325f06a8e531563326df1dcc34941775e75e28468e1c98fd0c12c2bb4dfd4aee3f55ec404aef48e8f5de5273645f26cd5d12726c08ef3f84845b38e54730b1");

struct fuzzy_pake_profile {
    double keygenA_time = 0.0, keygenB_time = 0.0, keyderA_time = 0.0, keyderB_time = 0.0;
    bool alice_only = false;
    keygen_profile keygenA;
    keygen_profile keygenB;
    keyder_profile keyderA;
    keyder_profile keyderB;
    double setup_time = 0.0; // Time taken for setup

    inline nlohmann::json to_json() const {
        nlohmann::json j;
        j["alice_only"] = alice_only;
        j["setup_time"] = setup_time;

        j["keygenA_time"] = keygenA_time;
        j["keyderA_time"] = keyderA_time;
        j["keygenA"] = keygenA.to_json();
        j["keyderA"] = keyderA.to_json();

        j["keygenB_time"] = keygenB_time;
        j["keygenB"] = keygenB.to_json();
        j["keyderB_time"] = keyderB_time;
        if (!alice_only) {
            j["keyderB"] = keyderB.to_json();
        }
        return j;
    }
};

void output_benchmark_result(const crs &crs, const fuzzy_pake_profile &profile) {
    // Output in machine-readable format
    mkhss::print_parameters(crs.mkhss_crs.params, std::cout);
    print_parameters(std::cout, crs.params);
    // #ifdef BASELINE_MKHSS
    // std::cout << "mkhss_scheme: baseline" << std::endl;
    // #else
    // std::cout << "mkhss_scheme: optimized" << std::endl;
    // #endif

    // std::cout << "anike_app: fuzzy_pake" << std::endl;
    // std::cout << "keygenA_time: " << keygenA_time << std::endl;
    // std::cout << "keygenB_time: " << keygenB_time << std::endl;
    // std::cout << "keyderA_time: " << keyderA_time << std::endl;
    // std::cout << "keyderB_time: " << keyderB_time << std::endl;
    // std::cout << "num_rms: " << num_rms << std::endl;
    nlohmann::json j;
    j["mkhss_scheme"] = mkhss::IS_BASELINE ? "baseline" : "optimized";
    j["anike_app"] = "fuzzy_pake";
    j["profile"] = profile.to_json();
    j["L"] = crs.params.L;
    j["W"] = crs.params.W;
    j["b"] = crs.params.b;
    j["T"] = crs.params.T;
    j["Q"] = crs.params.Q;
    std::cout << j.dump(4) << std::endl;
}

void encode_passphrase(const std::vector<std::string> &words, std::vector<bool> &encoded, int W, int b) {
    // Encode a passphrase into a binary vector
    // Each word is represented by W characters, each character is b bits
    assert(words.size() * W * b == encoded.size());
    for (size_t i = 0; i < words.size(); i++) {
        for (size_t j = 0; j < W; j++) {
            uint8_t c = (j < words[i].size()) ? words[i][j] : '\0'; // Fill with NUL if the word is shorter than W
            assert(c >> b == 0); // Ensure the character fits in b bits
            for (int k = 0; k < b; k++) {
                encoded[i * W * b + j * b + k] = (c >> (b - 1 - k)) & 1;
            }
        }
    }
}

void encode_passphrase_charset(const std::vector<std::string> &words, std::vector<bool> &encoded, const std::string &charset, int W, int b) {
    // Encode a passphrase into a binary vector using a custom charset
    // Each word is represented by W characters, each character is b bits
    assert(words.size() * W * b == encoded.size());
    for (size_t i = 0; i < words.size(); i++) {
        for (size_t j = 0; j < W; j++) {
            uint8_t c = (j < words[i].size()) ? words[i][j] : '\0'; // Fill with NUL if the word is shorter than W
            size_t index = charset.find(c);
            assert(index != std::string::npos && index < (1 << b)); // Ensure the character is in the charset and fits in b bits
            for (int k = 0; k < b; k++) {
                encoded[i * W * b + j * b + k] = (index >> (b - 1 - k)) & 1;
            }
        }
    }
}

int main(int argc, char **argv) {
    srand(time(nullptr));
    // Test the ANIKE parameters

    // Parameters suitable for eff_large_wordlist.txt
    // https://www.eff.org/deeplinks/2016/07/new-wordlists-random-passphrases
    // int length = 8; // L, number of words
    // int W = 9; // Max word length in number of characters
    // int b = 5; // b, length of each character in bits (words in eff_large_wordlist.txt contains 28 unique characters, so 5 bits is enough)
    // int threshold = 2; // T, the number of mismatching words allowed
    // int threshold_char = 2; // Q, the number of mismatching characters allowed in each word
    // // A fuzzy eight-word passphrase chosen with words from passphrase/final_filtered_eff_large_wordlist.txt
    // // provides ~51.128 bits of security against online guessing attacks
    // // roughly equivalent to four-word passphrase with exact PAKE (51.699)
    // std::string charset("\x00-abcdefghijklmnopqrstuvwxyz"s); // 28 unique characters, including NUL

    // std::cout << charset.size() << " unique characters in the charset" << std::endl;

    assert(argc == 6 || argc == 7 && "Usage: fuzzy_pake <L> <W> <b> <T> <Q> --alice-only");
    int length = atoi(argv[1]); // L, number of words
    int W = atoi(argv[2]); // Max word length in number of characters
    int b = atoi(argv[3]); // b, length of each character in bits
    int threshold = atoi(argv[4]); // T, the number of mismatching words allowed
    int threshold_char = atoi(argv[5]); // Q, the number of mismatching characters allowed in each word
    bool alice_only = (argc == 7 && std::string(argv[6]) == "--alice-only");

    fuzzy_pake_profile profile;

    crs crs;
    parameters params;
    public_encoding pe_A, pe_B;
    private_state st_A, st_B;
    create_parameters(params, 128, 128, 3072, threshold_char, threshold, W, b, length);
    std::cerr << "setup..." << std::flush;

    profile.setup_time = measure_time_seconds([&]() {
        setup(crs, params, test_n.get_mpz_t());
    });
    std::cerr << " done" << std::endl;
    mkhss::print_parameters(crs.mkhss_crs.params, std::cerr);
    print_parameters(std::cerr, params);

    std::vector<bool> x(params.total_length), y(params.total_length);
    
    for (int i = 0; i < params.total_length; i++) {
        x[i] = rand() % 2; // Random binary string
        y[i] = rand() % 2; // Random binary string
    }
    
    // encode_passphrase_charset({"cofounder", "glacial", "uninstall", "lapel", "kinship", "reggae", "dinginess", "reformat"}, x, charset, W, b);
    // encode_passphrase_charset({"cofoundre", "glacier", "reinstall", "babel", "kinsship", "beggar", "forgot", "reformed"}, y, charset, W, b);

    std::cerr << "x:";
    for (int i = 0; i < params.total_length; i++) {
        if (i % params.b == 0) {
            std::cerr << " ";
        }
        if (i % (params.b * params.W) == 0) {
            std::cerr << "| ";
        }
        std::cerr << x[i];
    }
    std::cerr << std::endl;
    std::cerr << "y:";
    for (int i = 0; i < params.total_length; i++) {
        if (i % params.b == 0) {
            std::cerr << " ";
        }
        if (i % (params.b * params.W) == 0) {
            std::cerr << "| ";
        }
        std::cerr << y[i];
    }
    std::cerr << std::endl;

    // int actual_distance = hamming_distance(x, y, length);
    // std::cerr << "Actual Hamming distance: " << actual_distance << std::endl;

    std::cerr << "attr_key_gen(sigma = A)... " << std::flush;
    double keygenA_time = measure_time_seconds([&]() {
        attr_key_gen(crs, x, pe_A, st_A, profile.keygenA);
    });
    std::cerr << keygenA_time << " seconds" << std::endl;
    std::cerr << "attr_key_gen(sigma = B)... " << std::flush;
    double keygenB_time = measure_time_seconds([&]() {
        attr_key_gen(crs, y, pe_B, st_B, profile.keygenB);
    });
    std::cerr << keygenB_time << " seconds" << std::endl;

    // std::cout << "Public encoding size (Alice): " << get_public_encoding_size(crs, pe_A) << std::endl;
    // std::cout << "Public encoding size (Bob): " << get_public_encoding_size(crs, pe_B) << std::endl;

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

    keyder_profile &keyder_prof_A = profile.keyderA, &keyder_prof_B = profile.keyderB;
    std::cerr << "attr_key_der(sigma = A)... " << std::flush;
    uint8_t k_A[KEY_LENGTH], k_B[KEY_LENGTH];
    uint64_t num_rms;
    double keyderA_time = measure_time_seconds([&]() {
        attr_key_der(crs, party_id::ALICE, st_A, pe_B, k_A, keyder_prof_A);
    });
    std::cerr << keyderA_time << " seconds, " << keyder_prof_A.num_rms_multiplications << " RMS multiplications" << std::endl;
    std::cerr << "attr_key_der(sigma = B)... " << std::flush;
    double keyderB_time = 0.0 / 0.0; // Initialize to NaN
    if (alice_only) {
        std::cerr << "Skipping key derivation for Bob (alice-only mode)" << std::endl;
    } else {
        keyderB_time = measure_time_seconds([&]() {
            attr_key_der(crs, party_id::BOB, st_B, pe_A, k_B, keyder_prof_B);
        });
        std::cerr << keyderB_time << " seconds, " << keyder_prof_A.num_rms_multiplications << " RMS multiplications" << std::endl;
    }

    std::cerr << std::endl;

    std::cerr << "Key A: ";
    for (size_t i = 0; i < KEY_LENGTH; i++) {
        std::cerr << std::hex << std::setfill('0') << std::setw(2) << (int) k_A[i];
    }
    std::cerr << std::dec << std::endl;
    if (!alice_only) {
        std::cerr << "Key B: ";
        for (size_t i = 0; i < KEY_LENGTH; i++) {
            std::cerr << std::hex << std::setfill('0') << std::setw(2) << (int) k_B[i];
        }
        std::cerr << std::dec << std::endl;

        bool are_keys_equal = memcmp(k_A, k_B, KEY_LENGTH) == 0;
        std::cerr << "Keys are " << (are_keys_equal ? "equal" : "not equal") << std::endl;
        // TODO: check against expected value
    }

    profile.keygenA_time = keygenA_time;
    profile.keygenB_time = keygenB_time;
    profile.keyderA_time = keyderA_time;
    profile.keyderB_time = keyderB_time;
    profile.alice_only = alice_only;
    output_benchmark_result(crs, profile);

    return 0;
}
