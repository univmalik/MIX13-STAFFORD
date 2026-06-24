// swirl512_round2_extended.cpp
//
// Extended round-2 statistical analysis with multiple complementary tests.
// All large statistics arrays are heap-allocated to avoid stack overflow.
//
// Tests performed (both unconditioned and conditioned on round-1 cluster event):
//   T1. Exact 512-bit delta_2 collision count (via 64-bit hash + sort + dedup)
//   T2. Per-lane delta_2 collision count (per-lane value + sort + dedup)
//   T3. Per-bit bias (512 bit counters, compare to N/2 with sigma = sqrt(N/4))
//   T4. Per-byte chi-square (64 byte positions × 256 values, 255 dof)
//   T5. Hamming weight distribution
//   T6. Top-K truncated differentials per lane (top byte, top 16-bit)
//
// Build:
//   g++ -O3 -std=c++20 -march=native -pthread swirl512_round2_extended.cpp -o r2ext
//
// Run:
//   ./r2ext [N_exp]  default N_exp = 26 (2^26 = 67M samples)

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <thread>
#include <vector>

using u8  = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using State = std::array<u64, 8>;

static constexpr std::array<unsigned, 8> MDS = {1, 1, 4, 1, 8, 5, 2, 9};

static inline u8 gf_mul2(u8 x) {
    return static_cast<u8>((x << 1) ^ ((x & 0x80U) ? 0x1dU : 0x00U));
}
static inline u8 gf_mul(u8 x, unsigned n) {
    switch (n) {
        case 1: return x;
        case 2: return gf_mul2(x);
        case 4: return gf_mul2(gf_mul2(x));
        case 5: return static_cast<u8>(gf_mul(x, 4) ^ x);
        case 8: return gf_mul2(gf_mul(x, 4));
        case 9: return static_cast<u8>(gf_mul(x, 8) ^ x);
        default: return 0;
    }
}

static inline u64 mix64(u64 z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline u8 byte_be(u64 x, int i) {
    return static_cast<u8>((x >> (56 - 8 * i)) & 0xffULL);
}

static inline u64 splitmix64(u64& x) {
    u64 z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline State linear_layer(const State& s) {
    State out{};
    for (int j = 0; j < 8; ++j) {
        u8 ob[8] = {0,0,0,0,0,0,0,0};
        for (int k = 0; k < 8; ++k) {
            u8 v = byte_be(s[k], j);
            if (v == 0) continue;
            for (int i = 0; i < 8; ++i) {
                ob[i] ^= gf_mul(v, MDS[(i - k) & 7]);
            }
        }
        out[j] = (static_cast<u64>(ob[0]) << 56) |
                 (static_cast<u64>(ob[1]) << 48) |
                 (static_cast<u64>(ob[2]) << 40) |
                 (static_cast<u64>(ob[3]) << 32) |
                 (static_cast<u64>(ob[4]) << 24) |
                 (static_cast<u64>(ob[5]) << 16) |
                 (static_cast<u64>(ob[6]) <<  8) |
                  static_cast<u64>(ob[7]);
    }
    return out;
}

static inline State round_F(const State& s) {
    State m;
    for (int i = 0; i < 8; ++i) m[i] = mix64(s[i]);
    return linear_layer(m);
}

static inline u64 pack8_scaled(u8 b) {
    return (static_cast<u64>(gf_mul(b, MDS[0])) << 56) |
           (static_cast<u64>(gf_mul(b, MDS[1])) << 48) |
           (static_cast<u64>(gf_mul(b, MDS[2])) << 40) |
           (static_cast<u64>(gf_mul(b, MDS[3])) << 32) |
           (static_cast<u64>(gf_mul(b, MDS[4])) << 24) |
           (static_cast<u64>(gf_mul(b, MDS[5])) << 16) |
           (static_cast<u64>(gf_mul(b, MDS[6])) <<  8) |
            static_cast<u64>(gf_mul(b, MDS[7]));
}

static State round1_output_for_beta(u64 beta_mix) {
    State out;
    for (int j = 0; j < 8; ++j) {
        out[j] = pack8_scaled(byte_be(beta_mix, j));
    }
    return out;
}

static inline u64 hash_state(const State& s) {
    u64 h = 0xcbf29ce484222325ULL;
    for (int i = 0; i < 8; ++i) {
        h ^= s[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}

static inline int popcount_state(const State& s) {
    int hw = 0;
    for (int i = 0; i < 8; ++i) hw += __builtin_popcountll(s[i]);
    return hw;
}

static constexpr u64 ALPHA_LANE0 = 0x8000000200000008ULL;

static const std::array<u64, 8> CLUSTER_BETAS = {
    0x53112150A62242A0ULL,
    0xD3112151A62242A0ULL,
    0xB3112151662242A0ULL,
    0xCD1121519A2242A0ULL,
    0x4D1121509A2242A0ULL,
    0xD5112151AA2242A0ULL,
    0xB51121516A2242A0ULL,
    0x55112150AA2242A0ULL,
};

// Use raw arrays wrapped in unique_ptr for large per-thread statistics.
// Smaller arrays are kept inline as std::array.
struct ThreadStats {
    u64 samples = 0;
    u64 cluster_events = 0;
    std::array<u64, 8> per_rank{};

    // 512 bits → 4KB (inline OK)
    std::array<u64, 512> bit_counts{};
    std::array<u64, 512> cond_bit_counts{};

    // 64 × 256 = 128KB (inline OK)
    std::array<std::array<u64, 256>, 64> byte_counts{};
    std::array<std::array<u64, 256>, 64> cond_byte_counts{};

    // 513 entries (small)
    std::array<u64, 513> hw_counts{};
    std::array<u64, 513> cond_hw_counts{};

    // 8 × 256 = 16KB (inline OK)
    std::array<std::array<u64, 256>, 8> lane_top8_counts{};
    std::array<std::array<u64, 256>, 8> cond_lane_top8_counts{};

    // 8 × 65536 × 8 = 4MB - heap allocate
    std::unique_ptr<std::array<std::array<u64, 65536>, 8>> lane_top16_counts;
    std::unique_ptr<std::array<std::array<u64, 65536>, 8>> cond_lane_top16_counts;

    // For collision detection
    std::vector<u64> all_hashes;
    std::vector<u64> cond_hashes;
    std::array<std::vector<u64>, 8> lane_values;
    std::array<std::vector<u64>, 8> cond_lane_values;

    ThreadStats()
        : lane_top16_counts(std::make_unique<std::array<std::array<u64, 65536>, 8>>()),
          cond_lane_top16_counts(std::make_unique<std::array<std::array<u64, 65536>, 8>>()) {
        // zero-init
        lane_top16_counts->fill({});
        cond_lane_top16_counts->fill({});
    }

    void reserve(u64 per_thread) {
        all_hashes.reserve(per_thread);
        cond_hashes.reserve(per_thread / 64 + 1);
        for (int i = 0; i < 8; ++i) {
            lane_values[i].reserve(per_thread);
            cond_lane_values[i].reserve(per_thread / 64 + 1);
        }
    }
};

int main(int argc, char** argv) {
    int N_exp = 26;
    if (argc >= 2) N_exp = std::atoi(argv[1]);
    if (N_exp < 18 || N_exp > 30) {
        std::fprintf(stderr, "N_exp out of range (18..30)\n");
        return 1;
    }
    const u64 N = 1ULL << N_exp;

    unsigned num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    if (num_threads > 8) num_threads = 8;

    std::printf("SWIRL512 Round-2 Extended Statistical Analysis\n");
    std::printf("===============================================\n");
    std::printf("Input alpha: lane 0 = 0x%016llX, lanes 1..7 = 0\n",
                (unsigned long long)ALPHA_LANE0);
    std::printf("Samples: 2^%d = %llu\n", N_exp, (unsigned long long)N);
    std::printf("Threads: %u\n\n", num_threads);

    // Precompute cluster output differences
    std::array<State, 8> cluster_outputs;
    for (int k = 0; k < 8; ++k) {
        cluster_outputs[k] = round1_output_for_beta(CLUSTER_BETAS[k]);
    }

    // Per-thread stats (heap allocated)
    std::vector<std::unique_ptr<ThreadStats>> thread_stats(num_threads);
    for (unsigned t = 0; t < num_threads; ++t) {
        thread_stats[t] = std::make_unique<ThreadStats>();
    }

    u64 per_thread = N / num_threads;
    u64 remainder = N % num_threads;

    auto worker = [&](unsigned tid) {
        ThreadStats& ts = *thread_stats[tid];
        u64 my_count = per_thread + (tid == 0 ? remainder : 0);
        ts.reserve(my_count);

        u64 sm = 0xDEADBEEF12345678ULL + tid * 0x1234567890ULL;

        for (u64 n = 0; n < my_count; ++n) {
            State x, xa;
            for (int i = 0; i < 8; ++i) x[i] = splitmix64(sm);
            xa[0] = x[0] ^ ALPHA_LANE0;
            for (int i = 1; i < 8; ++i) xa[i] = x[i];

            State F1x = round_F(x);
            State F1xa = round_F(xa);
            State delta1;
            for (int i = 0; i < 8; ++i) delta1[i] = F1x[i] ^ F1xa[i];

            int cluster_rank = -1;
            for (int k = 0; k < 8; ++k) {
                if (delta1 == cluster_outputs[k]) {
                    cluster_rank = k;
                    break;
                }
            }
            bool is_cluster = (cluster_rank >= 0);

            ts.samples++;
            if (is_cluster) {
                ts.cluster_events++;
                ts.per_rank[cluster_rank]++;
            }

            State F2x = round_F(F1x);
            State F2xa = round_F(F1xa);
            State delta2;
            for (int i = 0; i < 8; ++i) delta2[i] = F2x[i] ^ F2xa[i];

            u64 h = hash_state(delta2);
            ts.all_hashes.push_back(h);

            for (int i = 0; i < 8; ++i) {
                ts.lane_values[i].push_back(delta2[i]);
                ts.lane_top8_counts[i][(delta2[i] >> 56) & 0xFF]++;
                (*ts.lane_top16_counts)[i][(delta2[i] >> 48) & 0xFFFF]++;
            }

            for (int lane = 0; lane < 8; ++lane) {
                u64 v = delta2[lane];
                for (int b = 0; b < 64; ++b) {
                    if ((v >> (63 - b)) & 1ULL) {
                        ts.bit_counts[lane * 64 + b]++;
                    }
                }
            }

            for (int lane = 0; lane < 8; ++lane) {
                for (int bidx = 0; bidx < 8; ++bidx) {
                    u8 bv = byte_be(delta2[lane], bidx);
                    ts.byte_counts[lane * 8 + bidx][bv]++;
                }
            }

            int hw = popcount_state(delta2);
            ts.hw_counts[hw]++;

            if (is_cluster) {
                ts.cond_hashes.push_back(h);
                for (int i = 0; i < 8; ++i) {
                    ts.cond_lane_values[i].push_back(delta2[i]);
                    ts.cond_lane_top8_counts[i][(delta2[i] >> 56) & 0xFF]++;
                    (*ts.cond_lane_top16_counts)[i][(delta2[i] >> 48) & 0xFFFF]++;
                }
                for (int lane = 0; lane < 8; ++lane) {
                    u64 v = delta2[lane];
                    for (int b = 0; b < 64; ++b) {
                        if ((v >> (63 - b)) & 1ULL) {
                            ts.cond_bit_counts[lane * 64 + b]++;
                        }
                    }
                }
                for (int lane = 0; lane < 8; ++lane) {
                    for (int bidx = 0; bidx < 8; ++bidx) {
                        u8 bv = byte_be(delta2[lane], bidx);
                        ts.cond_byte_counts[lane * 8 + bidx][bv]++;
                    }
                }
                ts.cond_hw_counts[hw]++;
            }
        }
    };

    std::vector<std::thread> threads;
    auto t0 = std::chrono::steady_clock::now();
    for (unsigned t = 0; t < num_threads; ++t) {
        threads.emplace_back(worker, t);
    }
    for (auto& th : threads) th.join();
    auto t1 = std::chrono::steady_clock::now();
    double gen_elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::printf("Sample generation done in %.1fs\n\n", gen_elapsed);

    // Merge stats
    u64 total_samples = 0;
    u64 total_cluster = 0;
    std::array<u64, 8> total_per_rank{};
    std::array<u64, 512> total_bit_counts{};
    std::array<u64, 512> total_cond_bit_counts{};
    std::array<std::array<u64, 256>, 64> total_byte_counts{};
    std::array<std::array<u64, 256>, 64> total_cond_byte_counts{};
    std::array<u64, 513> total_hw_counts{};
    std::array<u64, 513> total_cond_hw_counts{};
    std::array<std::array<u64, 256>, 8> total_lane_top8{};
    std::array<std::array<u64, 256>, 8> total_cond_lane_top8{};
    auto total_lane_top16 = std::make_unique<std::array<std::array<u64, 65536>, 8>>();
    auto total_cond_lane_top16 = std::make_unique<std::array<std::array<u64, 65536>, 8>>();
    total_lane_top16->fill({});
    total_cond_lane_top16->fill({});

    for (auto& ts_ptr : thread_stats) {
        ThreadStats& ts = *ts_ptr;
        total_samples += ts.samples;
        total_cluster += ts.cluster_events;
        for (int k = 0; k < 8; ++k) total_per_rank[k] += ts.per_rank[k];
        for (int b = 0; b < 512; ++b) {
            total_bit_counts[b] += ts.bit_counts[b];
            total_cond_bit_counts[b] += ts.cond_bit_counts[b];
        }
        for (int p = 0; p < 64; ++p) {
            for (int v = 0; v < 256; ++v) {
                total_byte_counts[p][v] += ts.byte_counts[p][v];
                total_cond_byte_counts[p][v] += ts.cond_byte_counts[p][v];
            }
        }
        for (int h = 0; h <= 512; ++h) {
            total_hw_counts[h] += ts.hw_counts[h];
            total_cond_hw_counts[h] += ts.cond_hw_counts[h];
        }
        for (int lane = 0; lane < 8; ++lane) {
            for (int v = 0; v < 256; ++v) {
                total_lane_top8[lane][v] += ts.lane_top8_counts[lane][v];
                total_cond_lane_top8[lane][v] += ts.cond_lane_top8_counts[lane][v];
            }
            for (int v = 0; v < 65536; ++v) {
                (*total_lane_top16)[lane][v] += (*ts.lane_top16_counts)[lane][v];
                (*total_cond_lane_top16)[lane][v] += (*ts.cond_lane_top16_counts)[lane][v];
            }
        }
    }

    // ============ T1: Full delta_2 collision count ============
    std::printf("=== T1: Exact 512-bit Δ2 collision count (via 64-bit hash) ===\n\n");

    std::vector<u64> all_hashes;
    all_hashes.reserve(total_samples);
    for (auto& ts_ptr : thread_stats) {
        ThreadStats& ts = *ts_ptr;
        for (u64 h : ts.all_hashes) all_hashes.push_back(h);
        std::vector<u64>().swap(ts.all_hashes);
    }

    auto t2 = std::chrono::steady_clock::now();
    std::sort(all_hashes.begin(), all_hashes.end());
    auto t3 = std::chrono::steady_clock::now();
    std::printf("  Sort time: %.1fs\n", std::chrono::duration<double>(t3 - t2).count());

    u64 unique_count = 1;
    u64 max_count = 1;
    u64 current_count = 1;
    u64 collision_pairs = 0;
    u64 max_collision_hash = 0;
    u64 max_collision_count = 0;
    for (size_t i = 1; i < all_hashes.size(); ++i) {
        if (all_hashes[i] == all_hashes[i-1]) {
            current_count++;
            if (current_count > max_count) {
                max_count = current_count;
                max_collision_hash = all_hashes[i];
                max_collision_count = current_count;
            }
        } else {
            if (current_count > 1) collision_pairs += current_count - 1;
            current_count = 1;
            unique_count++;
        }
    }
    if (current_count > 1) collision_pairs += current_count - 1;

    std::printf("  Unconditioned:\n");
    std::printf("    Samples:                %llu\n", (unsigned long long)total_samples);
    std::printf("    Distinct hashes:        %llu\n", (unsigned long long)unique_count);
    std::printf("    Hash collisions:        %llu\n", (unsigned long long)collision_pairs);
    std::printf("    Max count for any hash: %llu\n", (unsigned long long)max_count);
    if (max_count > 1) {
        std::printf("    Most colliding hash:    0x%016llX (count=%llu)\n",
                    (unsigned long long)max_collision_hash,
                    (unsigned long long)max_collision_count);
    }
    double expected_collisions = (double)total_samples * (total_samples - 1.0) / std::pow(2.0, 64);
    std::printf("    Expected (uniform 64-bit hash): %.4f\n", expected_collisions);

    std::vector<u64> cond_hashes;
    cond_hashes.reserve(total_cluster);
    for (auto& ts_ptr : thread_stats) {
        ThreadStats& ts = *ts_ptr;
        for (u64 h : ts.cond_hashes) cond_hashes.push_back(h);
        std::vector<u64>().swap(ts.cond_hashes);
    }
    std::sort(cond_hashes.begin(), cond_hashes.end());

    u64 cond_unique = 1;
    u64 cond_max_count = 1;
    u64 cond_current = 1;
    u64 cond_collisions = 0;
    for (size_t i = 1; i < cond_hashes.size(); ++i) {
        if (cond_hashes[i] == cond_hashes[i-1]) {
            cond_current++;
            if (cond_current > cond_max_count) cond_max_count = cond_current;
        } else {
            if (cond_current > 1) cond_collisions += cond_current - 1;
            cond_current = 1;
            cond_unique++;
        }
    }
    if (cond_current > 1) cond_collisions += cond_current - 1;

    std::printf("  Conditioned on round-1 cluster:\n");
    std::printf("    Samples:                %llu (cluster events)\n", (unsigned long long)total_cluster);
    std::printf("    Distinct hashes:        %llu\n", (unsigned long long)cond_unique);
    std::printf("    Hash collisions:        %llu\n", (unsigned long long)cond_collisions);
    std::printf("    Max count for any hash: %llu\n", (unsigned long long)cond_max_count);
    double cond_expected = (double)total_cluster * (total_cluster - 1.0) / std::pow(2.0, 64);
    std::printf("    Expected (uniform):     %.6f\n", cond_expected);
    std::printf("\n");

    // ============ T2: Per-lane collision count ============
    std::printf("=== T2: Per-lane Δ2 collision count ===\n\n");

    double exp_col_uncond = (double)total_samples * (total_samples - 1.0) / std::pow(2.0, 64);
    double exp_col_cond = (double)total_cluster * (total_cluster - 1.0) / std::pow(2.0, 64);

    for (int lane = 0; lane < 8; ++lane) {
        std::vector<u64> lane_all;
        lane_all.reserve(total_samples);
        for (auto& ts_ptr : thread_stats) {
            ThreadStats& ts = *ts_ptr;
            for (u64 v : ts.lane_values[lane]) lane_all.push_back(v);
        }
        std::sort(lane_all.begin(), lane_all.end());

        u64 lu = 1, lm = 1, lc = 1, lcol = 0;
        for (size_t i = 1; i < lane_all.size(); ++i) {
            if (lane_all[i] == lane_all[i-1]) {
                lc++;
                if (lc > lm) lm = lc;
            } else {
                if (lc > 1) lcol += lc - 1;
                lc = 1;
                lu++;
            }
        }
        if (lc > 1) lcol += lc - 1;

        std::vector<u64> lane_cond;
        lane_cond.reserve(total_cluster);
        for (auto& ts_ptr : thread_stats) {
            ThreadStats& ts = *ts_ptr;
            for (u64 v : ts.cond_lane_values[lane]) lane_cond.push_back(v);
        }
        std::sort(lane_cond.begin(), lane_cond.end());

        u64 cu = 1, cm = 1, cc = 1, ccol = 0;
        for (size_t i = 1; i < lane_cond.size(); ++i) {
            if (lane_cond[i] == lane_cond[i-1]) {
                cc++;
                if (cc > cm) cm = cc;
            } else {
                if (cc > 1) ccol += cc - 1;
                cc = 1;
                cu++;
            }
        }
        if (cc > 1) ccol += cc - 1;

        std::printf("  Lane %d: uncond collisions=%llu (expected %.4f), max=%llu;  "
                    "cond collisions=%llu (expected %.6f), max=%llu\n",
                    lane,
                    (unsigned long long)lcol, exp_col_uncond, (unsigned long long)lm,
                    (unsigned long long)ccol, exp_col_cond, (unsigned long long)cm);

        std::vector<u64>().swap(lane_all);
        std::vector<u64>().swap(lane_cond);
    }
    for (auto& ts_ptr : thread_stats) {
        ThreadStats& ts = *ts_ptr;
        for (int lane = 0; lane < 8; ++lane) {
            std::vector<u64>().swap(ts.lane_values[lane]);
            std::vector<u64>().swap(ts.cond_lane_values[lane]);
        }
    }
    std::printf("\n");

    // ============ T3: Per-bit bias ============
    std::printf("=== T3: Per-bit bias (512 bits) ===\n\n");

    double bit_mean = total_samples / 2.0;
    double bit_sigma = std::sqrt(total_samples / 4.0);
    double bit_max_z = 0;
    int bit_max_idx = 0;
    double cond_bit_max_z = 0;
    int cond_bit_max_idx = 0;
    double cond_bit_mean = total_cluster / 2.0;
    double cond_bit_sigma = total_cluster > 0 ? std::sqrt(total_cluster / 4.0) : 1.0;

    for (int b = 0; b < 512; ++b) {
        double z = std::abs((double)total_bit_counts[b] - bit_mean) / bit_sigma;
        if (z > bit_max_z) { bit_max_z = z; bit_max_idx = b; }

        if (total_cluster > 0) {
            double cz = std::abs((double)total_cond_bit_counts[b] - cond_bit_mean) / cond_bit_sigma;
            if (cz > cond_bit_max_z) { cond_bit_max_z = cz; cond_bit_max_idx = b; }
        }
    }

    std::printf("  Unconditioned (N=%llu):\n", (unsigned long long)total_samples);
    std::printf("    Expected count per bit: %.1f\n", bit_mean);
    std::printf("    Sigma:                  %.1f\n", bit_sigma);
    std::printf("    Max |z| over 512 bits:  %.4f at bit %d (lane %d, bit %d)\n",
                bit_max_z, bit_max_idx, bit_max_idx / 64, bit_max_idx % 64);
    std::printf("    4-sigma threshold:      4.0  (significant if exceeded)\n");

    std::printf("  Conditioned on cluster (N=%llu):\n", (unsigned long long)total_cluster);
    std::printf("    Expected count per bit: %.1f\n", cond_bit_mean);
    std::printf("    Sigma:                  %.1f\n", cond_bit_sigma);
    std::printf("    Max |z| over 512 bits:  %.4f at bit %d (lane %d, bit %d)\n",
                cond_bit_max_z, cond_bit_max_idx, cond_bit_max_idx / 64, cond_bit_max_idx % 64);
    std::printf("\n");

    // ============ T4: Per-byte chi-square ============
    std::printf("=== T4: Per-byte chi-square (64 positions × 256 values, 255 dof) ===\n\n");

    double crit_255_001 = 330.52;  // chi-square critical value, 255 dof, p=0.001
    double max_chi2 = 0;
    int max_chi2_pos = 0;
    double cond_max_chi2 = 0;
    int cond_max_chi2_pos = 0;

    for (int p = 0; p < 64; ++p) {
        double expected = (double)total_samples / 256.0;
        double chi2 = 0;
        for (int v = 0; v < 256; ++v) {
            double diff = (double)total_byte_counts[p][v] - expected;
            chi2 += diff * diff / expected;
        }
        if (chi2 > max_chi2) { max_chi2 = chi2; max_chi2_pos = p; }

        if (total_cluster > 0) {
            double cond_expected = (double)total_cluster / 256.0;
            double cond_chi2 = 0;
            for (int v = 0; v < 256; ++v) {
                double diff = (double)total_cond_byte_counts[p][v] - cond_expected;
                cond_chi2 += diff * diff / cond_expected;
            }
            if (cond_chi2 > cond_max_chi2) { cond_max_chi2 = cond_chi2; cond_max_chi2_pos = p; }
        }
    }

    std::printf("  Unconditioned (N=%llu):\n", (unsigned long long)total_samples);
    std::printf("    Expected per bucket:    %.1f\n", total_samples / 256.0);
    std::printf("    Max chi2 over 64 bytes: %.2f at byte position %d (lane %d, byte %d)\n",
                max_chi2, max_chi2_pos, max_chi2_pos / 8, max_chi2_pos % 8);
    std::printf("    Critical (p=0.001):     %.2f\n", crit_255_001);
    std::printf("    Verdict:                %s\n", max_chi2 > crit_255_001 ? "BIAS DETECTED" : "no significant bias");

    std::printf("  Conditioned on cluster (N=%llu):\n", (unsigned long long)total_cluster);
    std::printf("    Expected per bucket:    %.1f\n", total_cluster / 256.0);
    std::printf("    Max chi2 over 64 bytes: %.2f at byte position %d (lane %d, byte %d)\n",
                cond_max_chi2, cond_max_chi2_pos, cond_max_chi2_pos / 8, cond_max_chi2_pos % 8);
    std::printf("    Critical (p=0.001):     %.2f\n", crit_255_001);
    std::printf("    Verdict:                %s\n", cond_max_chi2 > crit_255_001 ? "BIAS DETECTED" : "no significant bias");
    std::printf("\n");

    // ============ T5: Hamming weight distribution ============
    std::printf("=== T5: Hamming weight distribution ===\n\n");

    double obs_sum = 0;
    for (int h = 0; h <= 512; ++h) obs_sum += (double)h * total_hw_counts[h];
    double obs_mean = obs_sum / total_samples;
    double obs_var = 0;
    for (int h = 0; h <= 512; ++h) {
        double d = h - obs_mean;
        obs_var += d * d * total_hw_counts[h];
    }
    obs_var /= total_samples;
    double obs_std = std::sqrt(obs_var);

    int obs_min = -1, obs_max = -1;
    for (int h = 0; h <= 512; ++h) {
        if (total_hw_counts[h] > 0) {
            if (obs_min < 0) obs_min = h;
            obs_max = h;
        }
    }

    std::printf("  Unconditioned (N=%llu):\n", (unsigned long long)total_samples);
    std::printf("    Expected mean:   256.000  (Binomial(512, 0.5))\n");
    std::printf("    Expected std:    11.314\n");
    std::printf("    Observed mean:   %.3f\n", obs_mean);
    std::printf("    Observed std:    %.3f\n", obs_std);
    std::printf("    Observed range:  [%d, %d]\n", obs_min, obs_max);

    if (total_cluster > 0) {
        double cond_obs_sum = 0;
        for (int h = 0; h <= 512; ++h) cond_obs_sum += (double)h * total_cond_hw_counts[h];
        double cond_obs_mean = cond_obs_sum / total_cluster;
        double cond_obs_var = 0;
        for (int h = 0; h <= 512; ++h) {
            double d = h - cond_obs_mean;
            cond_obs_var += d * d * total_cond_hw_counts[h];
        }
        cond_obs_var /= total_cluster;
        double cond_obs_std = std::sqrt(cond_obs_var);
        int cond_min = -1, cond_max = -1;
        for (int h = 0; h <= 512; ++h) {
            if (total_cond_hw_counts[h] > 0) {
                if (cond_min < 0) cond_min = h;
                cond_max = h;
            }
        }
        std::printf("  Conditioned on cluster (N=%llu):\n", (unsigned long long)total_cluster);
        std::printf("    Expected mean:   256.000\n");
        std::printf("    Expected std:    11.314\n");
        std::printf("    Observed mean:   %.3f\n", cond_obs_mean);
        std::printf("    Observed std:    %.3f\n", cond_obs_std);
        std::printf("    Observed range:  [%d, %d]\n", cond_min, cond_max);
    }
    std::printf("\n");

    // ============ T6: Top-K truncated differentials per lane ============
    std::printf("=== T6: Truncated differentials per lane (8-bit and 16-bit prefixes) ===\n\n");

    double ideal_top8_log2 = -8.0;
    double ideal_top16_log2 = -16.0;

    std::printf("  Unconditioned (N=%llu):\n", (unsigned long long)total_samples);
    std::printf("    Lane | top-1 8-bit (count, log2 p) | top-1 16-bit (count, log2 p)\n");
    std::printf("    ---------------------------------------------------------------\n");
    for (int lane = 0; lane < 8; ++lane) {
        u64 max8 = 0; int max8v = 0;
        for (int v = 0; v < 256; ++v) {
            if (total_lane_top8[lane][v] > max8) {
                max8 = total_lane_top8[lane][v];
                max8v = v;
            }
        }
        double p8 = (double)max8 / total_samples;
        double log2p8 = std::log2(p8);

        u64 max16 = 0; int max16v = 0;
        for (int v = 0; v < 65536; ++v) {
            if ((*total_lane_top16)[lane][v] > max16) {
                max16 = (*total_lane_top16)[lane][v];
                max16v = v;
            }
        }
        double p16 = (double)max16 / total_samples;
        double log2p16 = std::log2(p16);

        std::printf("    %4d | 0x%02X (%llu, %.4f vs ideal %.3f) | 0x%04X (%llu, %.4f vs ideal %.3f)\n",
                    lane, max8v, (unsigned long long)max8, log2p8, ideal_top8_log2,
                    max16v, (unsigned long long)max16, log2p16, ideal_top16_log2);
    }

    if (total_cluster > 0) {
        std::printf("\n  Conditioned on cluster (N=%llu):\n", (unsigned long long)total_cluster);
        std::printf("    Lane | top-1 8-bit (count, log2 p) | top-1 16-bit (count, log2 p)\n");
        std::printf("    ---------------------------------------------------------------\n");
        for (int lane = 0; lane < 8; ++lane) {
            u64 max8 = 0; int max8v = 0;
            for (int v = 0; v < 256; ++v) {
                if (total_cond_lane_top8[lane][v] > max8) {
                    max8 = total_cond_lane_top8[lane][v];
                    max8v = v;
                }
            }
            double p8 = (double)max8 / total_cluster;
            double log2p8 = std::log2(p8);

            u64 max16 = 0; int max16v = 0;
            for (int v = 0; v < 65536; ++v) {
                if ((*total_cond_lane_top16)[lane][v] > max16) {
                    max16 = (*total_cond_lane_top16)[lane][v];
                    max16v = v;
                }
            }
            double p16 = (double)max16 / total_cluster;
            double log2p16 = std::log2(p16);

            std::printf("    %4d | 0x%02X (%llu, %.4f vs ideal %.3f) | 0x%04X (%llu, %.4f vs ideal %.3f)\n",
                        lane, max8v, (unsigned long long)max8, log2p8, ideal_top8_log2,
                        max16v, (unsigned long long)max16, log2p16, ideal_top16_log2);
        }
    }
    std::printf("\n");

    // ============ Final summary ============
    std::printf("=== Final Summary ===\n\n");
    std::printf("Total samples:           %llu (2^%d)\n", (unsigned long long)total_samples, N_exp);
    std::printf("Round-1 cluster events:  %llu (expected ~%.0f, ratio %.6f vs ideal 2^-7=%.6f)\n",
                (unsigned long long)total_cluster,
                (double)total_samples / 128.0,
                (double)total_cluster / total_samples,
                1.0 / 128.0);
    std::printf("  Per rank: ");
    for (int k = 0; k < 8; ++k) std::printf("%llu ", (unsigned long long)total_per_rank[k]);
    std::printf("\n\n");

    std::printf("Test results:\n");
    std::printf("  T1 Exact Δ2 collisions (uncond):  %llu (expected %.4f from 64-bit hash)\n",
                (unsigned long long)collision_pairs, expected_collisions);
    std::printf("  T1 Exact Δ2 collisions (cond):    %llu (expected %.4f)\n",
                (unsigned long long)cond_collisions, cond_expected);
    std::printf("  T3 Max bit bias |z| (uncond):     %.4f  (4-sigma = significant)\n", bit_max_z);
    std::printf("  T3 Max bit bias |z| (cond):       %.4f\n", cond_bit_max_z);
    std::printf("  T4 Max byte chi2 (uncond):        %.2f  (crit p=0.001: %.2f)\n", max_chi2, crit_255_001);
    std::printf("  T4 Max byte chi2 (cond):          %.2f  (crit p=0.001: %.2f)\n", cond_max_chi2, crit_255_001);
    std::printf("  T5 HW mean (uncond):              %.3f  (expected 256.000)\n", obs_mean);
    std::printf("  T5 HW std (uncond):               %.3f  (expected 11.314)\n", obs_std);

    return 0;
}
