#include <vector>
#include <iostream>
#include <cstdint>

// Precompute a small table for wheel factorization w.r.t. bases 2, 3, 5
std::vector<uint64_t> wheel_factors;

constexpr uint64_t product_of_bases(const uint64_t* bases, size_t size) {
    uint64_t product = 1;
    for (size_t i = 0; i < size; ++i) {
        product *= bases[i];
    }
    return product;
}

// constexpr uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23}; // Takes ~1s
constexpr uint64_t bases[] = {2, 3}; // Bases for wheel factorization
constexpr uint64_t product = product_of_bases(bases, sizeof(bases) / sizeof(bases[0])); // Product of bases

void precompute_wheel_factors() {
    for (uint64_t i = 1; i < product; ++i) {
        bool coprime = true;
        for (uint64_t base : bases) {
            if (i % base == 0 || i % base == ((base - 1) / 2)) {
                coprime = false;
                break;
            }
        }
        if (coprime) {
            wheel_factors.push_back(i);
        }
    }
}

// inline void escape(void *p) {
//     asm volatile("" : : "g"(p) : "memory");
// }

int main() {
    precompute_wheel_factors();
    // std::cout << "Size of wheel_factors: " << wheel_factors.size() << std::endl;
    // escape(&wheel_factors);
    // Print the precomputed wheel factors as a constexpr array
    std::cout << "constexpr uint64_t wheel[] = {";
    for (size_t i = 0; i < wheel_factors.size(); ++i) {
        std::cout << wheel_factors[i];
        if (i < wheel_factors.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "};" << std::endl;
    std::cout << "constexpr size_t wheel_size = " << wheel_factors.size() << ";" << std::endl;
    std::cout << "constexpr uint64_t wheel_product = " << product << ";" << std::endl;
    std::cout << std::endl;
    return 0;
}
