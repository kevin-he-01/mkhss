#pragma once
#include <iostream>
#include <fstream>
#include <stdint.h>
typedef unsigned __int128 uint128_t;
#include "vector/avx512ifma.hpp"

class PrecomputeIO {
    private:
        std::ifstream file;

    public:
        int read_counter;

        PrecomputeIO(const char* filename) : read_counter(0) {
            file = std::ifstream(filename);
            // Check for errors
            if (!file.is_open()) {
                std::cerr << "Error: Could not open precompute file: " << filename << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        template <typename t>
        void scalar(t &d) {
            if (file.eof()) {
                std::cerr << "Error: EOF reached while reading precompute file." << std::endl;
                exit(EXIT_FAILURE);
            }
            file >> d;
            read_counter++;
        }

        template <typename t, int limbs>
        void vector(t (&d)[limbs]) {
            for (int i = 0; i < limbs; i++) {
                if (file.eof()) {
                    std::cerr << "Error: EOF reached while reading precompute file." << std::endl;
                    exit(EXIT_FAILURE);
                }
                file >> d[i]; 
                read_counter++;
            }
        }

        template <int limbs>
        void read(AVXVector<limbs> &dst) {
            uint64_t v[nonzero(dst.LIMBS_PER_VEC * dst.VEC_LIMBS)];
            for (int j = 0; j < limbs; j++) {
                if (file.eof()) {
                    std::cerr << "Error: EOF reached while reading precompute file." << std::endl;
                    exit(EXIT_FAILURE);
                }
                file >> v[j];
                read_counter++;
            }
            dst.load((uint64_t*)&v);
        }
};
