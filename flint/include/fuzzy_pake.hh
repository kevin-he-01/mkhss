#ifndef FUZZY_PAKE_HH
#define FUZZY_PAKE_HH

// Attribute-based non-interactive key exchange (ANIKE) protocol for Fuzzy PAKE
#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>
#include <optional>
#include <functional>

#include "util.hh"
#include "random_oracle.hh"
#include "anike.hh"
#include "nike.hh"

#ifdef BASELINE_MKHSS
#include "baseline_mkhss.hh"
#define mkhss baseline_mkhss
#else
#include "mkhss.hh"
#endif

namespace anike::fuzzy_pake {
    using namespace anike;

    typedef mkhss::party_id party_id;

    struct parameters {
        int stat_sec_param; // Tau. Statistical security parameter (for correctness)
        int sec_param; // Lambda. Computational security parameter (against offline attacks)
        int n_width; // When lambda >= 128, need n_width >= 3072
        int Q, T, W, b, L;
        int total_length;
    };

    inline void create_parameters(parameters &params,
                                   int stat_sec_param,
                                   int sec_param,
                                   int min_n_width,
                                   int Q,
                                   int T,
                                   int W,
                                   int b,
                                   int L
                                ) {
        assert(stat_sec_param > 0);
        assert(sec_param > 0);
        assert(min_n_width > 0);
        params.stat_sec_param = stat_sec_param;
        params.sec_param = sec_param;
        params.n_width = min_n_width;
        params.Q = Q;
        params.T = T;
        params.W = W;
        params.b = b;
        params.L = L;
        params.total_length = L * W * b; // Total length of the binary string
    }

    inline void print_parameters(std::ostream &out, const parameters &params) {
        out << "# ===== ANIKE Parameters =====" << std::endl;
        out << "# stat_sec_param: " << params.stat_sec_param << std::endl;
        out << "# sec_param: " << params.sec_param << std::endl;
        out << "# n_width: " << params.n_width << std::endl;
        out << "# Q (character typo threshold): " << params.Q << std::endl;
        out << "# T (word mismatch threshold): " << params.T << std::endl;
        out << "# W (# characters in a word): " << params.W << std::endl;
        out << "# b (character length in bits): " << params.b << std::endl;
        out << "# L (Number of words): " << params.L << std::endl;
        out << "# total_length (Total length of the binary string): " << params.total_length << std::endl;
        out << "# ============================" << std::endl << std::endl;
    }

    struct crs {
        parameters params;
        mkhss::crs mkhss_crs; // MKHSS CRS
        nike::crs nike_crs; // NIKE CRS
    };

    inline void setup(crs &crs, const parameters &params, mpz_srcptr n = nullptr) {
        mkhss::parameters mkhss_params;
        crs.params = params;
        mkhss::create_parameters(mkhss_params,
                                   true, // short_exponent_assumption
                                   params.stat_sec_param,
                                   params.sec_param,
                                   params.n_width,
                                   1); // Max integer size : 1
        mkhss::setup(crs.mkhss_crs, mkhss_params, n);
        nike::setup(crs.nike_crs);
    }

    struct public_encoding { // pe_sigma
        mkhss::public_key mkhss_pk; // MKHSS public key
        std::vector<mkhss::public_share> x_i; // Shares of the bits of the binary string

        nike::public_key nike_pk; // NIKE public key
    };

    struct private_state { // st_sigma
        mkhss::private_key mkhss_sk; // MKHSS private key
        std::vector<mkhss::private_share> x_i; // Shares of the bits of the binary string

        nike::private_key nike_sk; // NIKE private key
    };

    // Synchronized shares for one party's strings
    struct synced_shares {
        std::vector<mkhss::input_share> x;
    };

    inline void get_pointers(const crs &crs, synced_shares &shares, std::vector<mkhss::input_share *> &pointers) {
        size_t expected_size = crs.params.total_length; // Total length of the binary string

        pointers.clear();
        pointers.reserve(expected_size);
        shares.x.resize(expected_size);
        for (auto &share : shares.x) {
            pointers.push_back(&share);
        }

        assert(pointers.size() == expected_size);
    }

    struct keygen_profile {
        double keygen_time = 0.0; // Time taken for MKHSS::KeyGen
        double precompute_time = 0.0; // Time taken for precomputation
        double share_time = 0.0; // Time taken for calling MKHSS::Share
        double nike_time = 0.0; // Time taken for NIKE key generation

        inline nlohmann::json to_json() const {
            nlohmann::json j;
            j["keygen_time"] = keygen_time;
            j["precompute_time"] = precompute_time;
            j["share_time"] = share_time;
            j["nike_time"] = nike_time;
            return j;
        }
    };

    inline void attr_key_gen(const crs &crs, const std::vector<bool> &x, public_encoding &pe, private_state &st,
                            keygen_profile &profile) {
        assert(x.size() == crs.params.total_length);
        FMPZ x_i;
        profile.keygen_time = measure_time_seconds([&]() {
            mkhss::keygen(crs.mkhss_crs, pe.mkhss_pk, st.mkhss_sk);
        });
        profile.precompute_time = measure_time_seconds([&]() {
            mkhss::precompute_share(crs.mkhss_crs, pe.mkhss_pk, crs.params.total_length); // Precompute to speed up MKHSS::Share
        });
        // Generate the shares of the binary string
        pe.x_i.clear();
        pe.x_i.reserve(crs.params.total_length);
        st.x_i.clear();
        st.x_i.reserve(crs.params.total_length);
        profile.share_time = measure_time_seconds([&]() {
            for (int i = 0; i < crs.params.total_length; i++) {
                x_i.set(x[i] ? 1 : 0);
                pe.x_i.emplace_back();
                st.x_i.emplace_back();
                mkhss::public_share &x_i_pk = pe.x_i.back();
                mkhss::private_share &x_i_sk = st.x_i.back();
                mkhss::share(crs.mkhss_crs, pe.mkhss_pk, x_i, x_i_sk, x_i_pk);
            }
        });
        // NIKE
        profile.nike_time = measure_time_seconds([&]() {
            nike::keygen(crs.nike_crs, st.nike_sk, pe.nike_pk);
        });
    }

    // length: L
    // threshold: T
    // To generate predictions of how many times each input share will be used (for precomputation), we need to do a dry run
    inline void dry_run(size_t start, size_t end, size_t length, size_t threshold,
                                const std::function<void(size_t, size_t)> &differ) {
        size_t char_length = end - start; // [start, end-1] and 0-based indexing
        assert(char_length % length == 0);
        size_t stride = char_length / length; // Number of bits per character
        for (int i = 1; i <= length; i++) {
            for (int j = 0; j <= threshold; j++) {
                if (j <= i - 1 && j >= i - ((int) length - (int) threshold)) {
                    differ(start + (i - 1) * stride, start + i * stride);
                }
            }
        }
    }

    typedef size_t node_index;

    struct node {
        // Parents are dependencies of this node in the computation graph
        std::vector<std::pair<int, node_index>> rms_mult_parents; // int is index of the input share
        std::vector<std::pair<int, node_index>> linear_parents; // int is the weight of the linear parent
        int in_degree = 0; // In-degree of the node, i.e., number of parents
        std::vector<node_index> children; // Children are nodes whose values are affected by this node
        int type = 2; // Type of the node, 2 means normal node, 0 means zero node, 1 means one node, 3 means RMS multiplication node
        mkhss::memory_share value; // Value of the node
        bool processed = false; // Whether the node has been processed in the evaluation phase
    };

    struct keyder_profile {
        double rms_init_time = 0.0; // Time taken for RMSInit
        mkhss::sync_shares_profile sync_shares_profile; // Time taken for MKHSS::SyncShares
        double rms_isub = 0.0; // Time taken for RMS subtraction of input shares
        double rms_mult = 0.0; // Time taken for RMS multiplications
        double rms_precompute = 0.0; // Time taken for RMS precomputation
        double memory_management = 0.0; // Time taken for memory management (e.g., clearing precomputed data)
        double key_derivation_time = 0.0; // Time taken for key derivation + hashing
        uint64_t num_rms_multiplications = -1; // Number of RMS multiplications performed

        inline nlohmann::json to_json() const {
            nlohmann::json j;
            j["rms_init_time"] = rms_init_time;
            j["sync_shares_profile"] = sync_shares_profile.to_json();
            j["rms_isub"] = rms_isub;
            j["rms_mult"] = rms_mult;
            j["rms_precompute"] = rms_precompute;
            j["memory_management"] = memory_management;
            j["key_derivation_time"] = key_derivation_time;

            j["num_rms_muls"] = num_rms_multiplications;
            return j;
        }
    };

    class ComputationGraph {
    // private:
    public:
        std::vector<node> nodes; // Nodes in the computation graph
        std::vector<node_index> ready_stack; // Stack of nodes that are ready to be processed (i.e., in-degree is 0)
        std::vector<node_index> blocked_stack; // Stack of nodes that are blocked. To be processed by the next RMS batch multiplication on the same input share
        int last_input_index = -1; // Last input index used in mult() to save precomputation

    public:
        ComputationGraph() {
            // Initialize the computation graph with a zero node and a one node
            node &zero_node = nodes.emplace_back();
            zero_node.type = 0; // Zero node
            ready_stack.push_back(0);
            node &one_node = nodes.emplace_back();
            one_node.type = 1; // One node
            ready_stack.push_back(1);
        }

        inline node_index zero() {
            return 0; // Zero node is the first node
        }

        inline node_index one() {
            return 1; // One node is the second node
        }

        inline node_index sum(node_index parent_a, node_index parent_b) {
            node_index result_index = nodes.size(); // Index of the new node
            node &result = nodes.emplace_back();
            result.linear_parents.emplace_back(1, parent_a); // Add parent_a with weight 1
            result.linear_parents.emplace_back(1, parent_b); // Add parent_b with weight 1
            result.in_degree = 2; // In-degree is 2 since it has two parents
            nodes[parent_a].children.push_back(result_index); // Add result as a child of parent_a
            nodes[parent_b].children.push_back(result_index); // Add result as a child of parent_b
            return result_index; // Return the new node
        }

        inline node_index sub(node_index parent_a, node_index parent_b) {
            node_index result_index = nodes.size(); // Index of the new node
            node &result = nodes.emplace_back();
            result.linear_parents.emplace_back(1, parent_a); // Add parent_a with weight 1
            result.linear_parents.emplace_back(-1, parent_b); // Add parent_b with weight -1
            result.in_degree = 2; // In-degree is 2 since it has two parents
            nodes[parent_a].children.push_back(result_index); // Add result as a child of parent_a
            nodes[parent_b].children.push_back(result_index); // Add result as a child of parent_b
            return result_index; // Return the new node
        }

        inline node_index mult(int input_index, node_index parent) {
            // std::cerr << "mult(" << input_index << ", " << parent << ")" << std::endl;
            node_index result_index = nodes.size(); // Index of the new node
            node &result = nodes.emplace_back();
            result.rms_mult_parents.emplace_back(input_index, parent); // Add parent with input index
            result.in_degree = 1; // In-degree is 1 since it has one parent
            result.type = 3; // RMS multiplication node
            nodes[parent].children.push_back(result_index); // Add result as a child of parent
            return result_index; // Return the new node
        }

        bool process(const mkhss::crs &crs, mkhss::rms_context &ctx, node_index node) {
            if (nodes[node].processed) {
                return false; // If already processed, skip
            }
            bool success = false; // Success flag for processing
            switch(nodes[node].type) {
                case 0:
                    mkhss::rms_mzero(nodes[node].value); // Zero node, set value to zero
                    success = true;
                    break;
                case 1:
                    mkhss::rms_mset(nodes[node].value, ctx.get_one()); // One node, set value to one
                    success = true;
                    break;
                case 2:
                    // Normal node, compute the value based on linear parents
                    mkhss::rms_mzero(nodes[node].value); // Initialize value to zero
                    for (const auto &[weight, parent] : nodes[node].linear_parents) {
                        mkhss::memory_share tmp;
                        mkhss::rms_mcmult(nodes[parent].value, weight, tmp); // Multiply the parent's value by the weight
                        mkhss::rms_madd(nodes[node].value, tmp, nodes[node].value); // Add the result to the node's value
                    }
                    success = true;
                    break;
                case 3:
                    blocked_stack.push_back(node); // RMS multiplication node, add to blocked stack
                    success = false;
                    break;
                default:
                    assert(false && "Unknown node type"); // Unknown node type, assert
            }
            nodes[node].processed = success;
            return success;
        }

        bool done() {
            return ready_stack.empty() && blocked_stack.empty(); // Check if both stacks are empty
        }

        void evaluate(const mkhss::crs &crs, mkhss::rms_context &ctx) {
            // Process nodes in the ready stack
            while (!ready_stack.empty()) {
                node_index current = ready_stack.back();
                ready_stack.pop_back();
                bool success = process(crs, ctx, current); // Process the current node
                if (!success) {
                    continue;
                }
                // Process all children of the current node
                for (node_index child : nodes[current].children) {
                    if (nodes[child].processed) {
                        continue; // If child is already processed, skip
                    }
                    // Decrease in-degree of the child
                    assert(nodes[child].in_degree > 0); // Ensure in-degree is positive
                    nodes[child].in_degree--;
                    if (nodes[child].in_degree == 0) {
                        // If in-degree is 0, add to ready stack
                        ready_stack.push_back(child);
                    }
                }
                // Clear children of the current node to free memory
                nodes[current].children.clear();
            }
        }

        bool process_blocked(const crs &anike_crs, const mkhss::crs &crs, mkhss::rms_context &ctx, std::vector<mkhss::input_share> &x_A_minus_x_B, const std::vector<int> &usage_counts, keyder_profile& profile) {
            // Print the blocked stack for debugging
            // auto &g = *this;
            // std::cout << "Size of blocked stack: " << g.blocked_stack.size() << std::endl;
            // for (node_index node : g.blocked_stack) {
            //     assert(g.nodes[node].type == 3); // Ensure the node is an RMS multiplication node
            //     for (const auto &[input_index, parent] : g.nodes[node].rms_mult_parents) {
            //         std::cerr << input_index << g.nodes[node].processed << ", ";
            //         // std::cerr << "Blocked RMS multiplication on input share " << input_index << std::endl;
            //     }
            // }
            
            if (blocked_stack.empty()) {
                return false; // No blocked nodes to process
            }
            node_index first = blocked_stack.front(); // Get the first blocked node
            assert(nodes[first].type == 3); // Ensure the node is an RMS multiplication node
            assert(!nodes[first].rms_mult_parents.empty()); // Ensure there are input shares
            int share_index = nodes[first].rms_mult_parents.front().first; // Get the input share index
            if (share_index != last_input_index) {
                if (last_input_index != -1) {
                    assert(share_index > last_input_index); // Ensure the share index is increasing
                    // Clear precomputed data for the previous input share
                    profile.memory_management += measure_time_seconds([&]() {
                        x_A_minus_x_B[last_input_index].clear_precomputed_data();
                    });
                }
                profile.rms_precompute += measure_time_seconds([&]() {
                    // Precompute the input share for the current share index
                    x_A_minus_x_B[share_index].precompute(crs, usage_counts[share_index], bitcnt_n(anike_crs.params.total_length) + 1);
                });
                last_input_index = share_index; // Update the last input index
            }

            for (node_index node : blocked_stack) {
                assert(nodes[node].type == 3); // Ensure the node is an RMS multiplication node
                assert(!nodes[node].rms_mult_parents.empty()); // Ensure there are input shares
                assert(nodes[node].rms_mult_parents.front().first == share_index); // Ensure all nodes in the blocked stack use the same input share
                if (nodes[node].processed) {
                    continue; // If the node is already processed, skip
                }
                mkhss::memory_share &value = nodes[node].value;
                mkhss::memory_share &parent_value = nodes[nodes[node].rms_mult_parents.front().second].value;
                profile.rms_mult += measure_time_seconds([&]() {
                    // Perform the RMS multiplication
                    ctx.mult_serial(crs, x_A_minus_x_B[share_index], parent_value, value);
                });
                nodes[node].processed = true; // Mark the node as processed
                for (node_index child : nodes[node].children) {
                    if (nodes[child].processed) {
                        continue; // If child is already processed, skip
                    }
                    // Decrease in-degree of the child
                    assert(nodes[child].in_degree > 0); // Ensure in-degree is positive
                    nodes[child].in_degree--;
                    if (nodes[child].in_degree == 0) {
                        // If in-degree is 0, add to ready stack
                        ready_stack.push_back(child);
                    }
                }
            }

            blocked_stack.clear(); // Clear the blocked stack after processing
            return true;
        }
    };

    inline node_index differ_bit_g(ComputationGraph &g, int x_minus_y, node_index z_i) {
        // Base case supporting single bit comparisons
        node_index tmp = g.mult(x_minus_y, z_i); // tmp = x_minus_y * z_i
        return g.mult(x_minus_y, tmp); // z_o = x_minus_y * tmp
    }

    inline node_index exceeds_threshold_g(ComputationGraph &g, size_t start, size_t end, size_t length, size_t threshold,
                                    node_index z_i, const std::function<node_index (size_t, size_t, node_index )> &differ) {
        size_t char_length = end - start; // [start, end-1] and 0-based indexing
        assert(char_length % length == 0);
        size_t stride = char_length / length; // Number of bits per character
        std::vector<node_index > acc_i;
        acc_i.resize(threshold + 1);
        acc_i[0] = z_i; // Initialize acc_i[0] with z_i
        for (size_t j = 1; j <= threshold; j++) {
            acc_i[j] = g.zero(); // Initialize acc_i[j] with zero
        }
        node_index acc_overflow_flag = g.zero(); // Initialize overflow flag
        for (int i = 1; i <= length; i++) {
            node_index z_prev = g.zero(); // Previous result
            for (int j = 0; j <= threshold; j++) {
                node_index z;
                if (j <= i - 1 && j >= i - ((int) length - (int) threshold)) {
                    z = differ(start + (i - 1) * stride, start + i * stride, acc_i[j]);
                } else {
                    z = g.zero(); // If not in the range, set z to zero
                }
                node_index zbar = g.sub(acc_i[j], z);
                acc_i[j] = g.sum(z_prev, zbar); // acc_i[j] = z_prev + zbar
                z_prev = z;
            }
            acc_overflow_flag = g.sum(acc_overflow_flag, z_prev); // acc_overflow_flag += z_prev
        }
        return acc_overflow_flag; // Return the overflow flag
    }

    inline void differ_bit(const mkhss::crs &crs, mkhss::rms_context &ctx, const mkhss::input_share &x_minus_y, const mkhss::memory_share &z_i, mkhss::memory_share &z_o) {
        // Base case supporting single bit comparisons
        mkhss::memory_share tmp;
        ctx.mult_serial(crs, x_minus_y, z_i, tmp);
        ctx.mult_serial(crs, x_minus_y, tmp, z_o);
    }

    // length: L
    // threshold: T
    inline void exceeds_threshold(size_t start, size_t end, size_t length, size_t threshold,
                                  const mkhss::memory_share &z_i, mkhss::memory_share &z_o,
                                  const std::function<void(size_t, size_t, const mkhss::memory_share &, mkhss::memory_share &)> &differ) {
        size_t char_length = end - start; // [start, end-1] and 0-based indexing
        assert(char_length % length == 0);
        size_t stride = char_length / length; // Number of bits per character
        std::vector<mkhss::memory_share> acc_i;
        acc_i.resize(threshold + 1);
        mkhss::rms_mset(acc_i[0], z_i); // Initialize acc_i[0] with z_i
        for (size_t j = 1; j <= threshold; j++) {
            mkhss::rms_mzero(acc_i[j]);
        }
        mkhss::memory_share acc_overflow_flag;
        mkhss::rms_mzero(acc_overflow_flag);
        for (int i = 1; i <= length; i++) {
            mkhss::memory_share z_prev, z, zbar;
            mkhss::rms_mzero(z_prev);
            for (int j = 0; j <= threshold; j++) {
                if (j <= i - 1 && j >= i - ((int) length - (int) threshold)) {
                    differ(start + (i - 1) * stride, start + i * stride, acc_i[j], z);
                } else {
                    mkhss::rms_mzero(z);
                }
                mkhss::rms_msub(acc_i[j], z, zbar);
                mkhss::rms_madd(z_prev, zbar, acc_i[j]);
                mkhss::rms_mset(z_prev, z);
            }
            mkhss::rms_madd(acc_overflow_flag, z_prev, acc_overflow_flag);
        }
        mkhss::rms_mset(z_o, acc_overflow_flag);
    }

    inline void attr_key_der(const crs &crs, party_id sigma,
                                    const private_state &st_self, const public_encoding &pe_other,
                                    uint8_t key[KEY_LENGTH], keyder_profile& profile) {
        mkhss::rms_context ctx;
        profile.rms_init_time = measure_time_seconds([&]() {
            mkhss::rms_bootstrap(crs.mkhss_crs, sigma, st_self.mkhss_sk, pe_other.mkhss_pk, ctx);
        });
        
        synced_shares x_A, x_B;
        std::vector<mkhss::input_share *> inputs_A_pointers, inputs_B_pointers;
        get_pointers(crs, x_A, inputs_A_pointers);
        get_pointers(crs, x_B, inputs_B_pointers);
        profile.sync_shares_profile = mkhss::sync_shares(crs.mkhss_crs, sigma, ctx, st_self.mkhss_sk, pe_other.mkhss_pk,
                           st_self.x_i,
                           inputs_A_pointers,
                           pe_other.x_i,
                           inputs_B_pointers);
        
        mkhss::memory_share z;
        mkhss::rms_mset(z, ctx.get_one());

        std::vector<int> usage_counts(crs.params.total_length, 0);
        dry_run(0, crs.params.W * crs.params.b * crs.params.L, crs.params.L, crs.params.T,
                        [&](size_t start, size_t end) {
            assert(end - start == crs.params.W * crs.params.b);
            dry_run(start, end, crs.params.W, crs.params.Q,
                            [&](size_t start, size_t end) {
                assert(end - start == crs.params.b);
                dry_run(start, end, crs.params.b, 0,
                                 [&](size_t start, size_t end) {
                    assert(end - start == 1);
                    usage_counts[start] += 2; // differ_bit uses two RMS multiplications
                });
            });
        });
        
        std::vector<mkhss::input_share> x_A_minus_x_B;
        x_A_minus_x_B.reserve(crs.params.total_length);
        
        profile.rms_isub = measure_time_seconds([&]() {
            for (size_t i = 0; i < crs.params.total_length; i++) {
                auto &share = x_A_minus_x_B.emplace_back();
                mkhss::rms_isub(crs.mkhss_crs, x_A.x[i], x_B.x[i], share);
                // share.precompute(crs.mkhss_crs, usage_counts[i], bitcnt_n(crs.params.total_length) + 1); // Precompute FBE for the expected usage count
            }
        });

        // double rms_mult_total = 0.0;
        // double rms_precompute_total = 0.0;
        // double memory_magagement_total = 0.0;

        // TODO: can significantly save precomputed memory by computing a computation graph (DAG) of the branching program,
        // and then batch all RMS multiplications involving the same input share together
        // A possibly more natural alternative that avoids significant code restructuring is to use C++ 20 coroutines
        // intmax_t prev_start = -1;
        // intmax_t prev_end = -1;
        // exceeds_threshold(0, crs.params.W * crs.params.b * crs.params.L, crs.params.L, crs.params.T,
        //                 z, z, [&](size_t start, size_t end, const mkhss::memory_share &z_i, mkhss::memory_share &z_o) {
        //     assert(end - start == crs.params.W * crs.params.b);
        //     assert(((intmax_t) start) >= prev_start); // Due to iteration order inside exceeds_threshold, not true for more nested calls
        //     if (((intmax_t) start) != prev_start) {
        //         if (prev_start != -1) {
        //             // Uncompute FBE for the previous segment (to save memory)
        //             // std::cerr << "Uncomputing segment [" << prev_start << ", " << prev_end << ")" << std::endl;
        //             memory_magagement_total += measure_time_seconds([&] {
        //                 for (size_t i = prev_start; i < prev_end; i++) {
        //                     x_A_minus_x_B[i].clear_precomputed_data(); // Clear precomputed data
        //                 }
        //             });
        //         }
        //         // std::cerr << "Precomputing segment [" << start << ", " << end << ")" << std::endl;
        //         rms_precompute_total += measure_time_seconds([&] {
        //             for (size_t i = start; i < end; i++) {
        //                 x_A_minus_x_B[i].precompute(crs.mkhss_crs, usage_counts[i], bitcnt_n(crs.params.total_length) + 1); // Precompute FBE for the expected usage count
        //             }
        //         });
        //     }
        //     prev_start = start;
        //     prev_end = end;
        //     rms_mult_total += measure_time_seconds([&] {
        //         exceeds_threshold(start, end, crs.params.W, crs.params.Q,
        //                         z_i, z_o, [&](size_t start, size_t end, const mkhss::memory_share &z_i, mkhss::memory_share &z_o) {
        //             assert(end - start == crs.params.b);
        //             exceeds_threshold(start, end, crs.params.b, 0,
        //                             z_i, z_o, [&](size_t start, size_t end, const mkhss::memory_share &z_i, mkhss::memory_share &z_o) {
        //                 assert(end - start == 1);
        //                 differ_bit(crs.mkhss_crs, ctx, x_A_minus_x_B[start], z_i, z_o);
        //             });
        //         });
        //     });
        // });

        ComputationGraph g;
        // Create a computation graph for the RMS multiplications
        node_index z_o = exceeds_threshold_g(g, 0, crs.params.W * crs.params.b * crs.params.L, crs.params.L, crs.params.T,
                        g.one(), [&](size_t start, size_t end, node_index z_i) {
            assert(end - start == crs.params.W * crs.params.b);
            return exceeds_threshold_g(g, start, end, crs.params.W, crs.params.Q,
                            z_i, [&](size_t start, size_t end, node_index z_i) {
                assert(end - start == crs.params.b);
                return exceeds_threshold_g(g, start, end, crs.params.b, 0,
                                z_i, [&](size_t start, size_t end, node_index z_i) {
                    assert(end - start == 1);
                    // std::cout << "[" << start << ", " << end << ")";
                    return differ_bit_g(g, start, z_i);
                });
            });
        });
        
        while (!g.nodes[z_o].processed) {
            // Process the ready nodes in the computation graph
            g.evaluate(crs.mkhss_crs, ctx);
            // Process the blocked nodes in the computation graph
            g.process_blocked(crs, crs.mkhss_crs, ctx, x_A_minus_x_B, usage_counts, profile);
        }
        mkhss::rms_mset(z, g.nodes[z_o].value); // Set the final value of z

        // Print usage counts of all input shares
        // for (size_t i = 0; i < crs.params.total_length; i++) {
        //     assert(x_A_minus_x_B[i].usage_count == usage_counts[i]);
        //     std::cerr << x_A_minus_x_B[i].usage_count << " ";
        // }
        // std::cerr << std::endl;
        // std::cerr << "DEBUG: # RMS multiplications: " << ctx.get_next_inst_id() << std::endl;

        profile.key_derivation_time = measure_time_seconds([&]() {
            uint8_t k_nike[nike::SHARED_KEY_LEN];
            nike::keyder(crs.nike_crs, st_self.nike_sk, pe_other.nike_pk, k_nike);
            derive_key(k_nike, z, key);
        });

        profile.num_rms_multiplications = ctx.get_next_inst_id(); // Return the number of RMS multiplications performed
        // profile.rms_mult = rms_mult_total;
        // profile.rms_precompute = rms_precompute_total;
        // profile.memory_management = memory_magagement_total;
    }
}

#undef mkhss

#endif 
