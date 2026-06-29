// swirl_round2_distinguisher_attack.cpp
//
// Standalone verifier for the reduced two-round SWIRL differential-linear
// distinguisher found from Mix13 structure.
//
// Reduced round model:
//   round_F(state):
//     1. Apply Stafford Mix13 independently to each 64-bit lane.
//     2. Apply the byte-transpose / Whirlpool-style MDS linear layer.
//
// Deliberately omitted:
//   - round constants
//   - round keys
//
// For XOR-difference experiments, fixed XOR constants cancel. Fixed XOR round
// keys also preserve the paired input difference; under uniform input sampling,
// they translate both paired states by the same value and do not change the
// measured differential distribution.
//
// Distinguisher:
//   alpha lane 0 = 0x0200000008000000, lanes 1..7 = 0
//   Delta2 = F(F(x)) xor F(F(x xor alpha))
//
//   Test whether:
//     parity((Delta2[2] & 0x72ae08c4e7e81a1a) ^
//            (Delta2[6] & 0xe55c1188cfd03434)) == 0
//
// Expected for a random 512-bit permutation:
//   Pr[test passes] = 0.5
//
// Observed in large test:
//   Pr[test passes] ~= 0.51826613 over 16,000,000 pairs.
//
// Build:
//   g++ -O3 -std=c++20 -march=native swirl_round2_distinguisher_attack.cpp -o r2dist
//
// Run:
//   ./r2dist                 # default N = 2^24 = 16,777,216 pairs
//   ./r2dist 23              # N = 2^23 pairs
//   ./r2dist 24 --round3     # also test carrying the same mask to round 3

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using u8 = std::uint8_t;
using u64 = std::uint64_t;
using State = std::array<u64, 8>;

static constexpr std::array<unsigned, 8> MDS = {1, 1, 4, 1, 8, 5, 2, 9};

static constexpr u64 ALPHA_LANE0 = 0x0200000008000000ULL;
static constexpr u64 FINAL_MASK_LANE2 = 0x72ae08c4e7e81a1aULL;
static constexpr u64 FINAL_MASK_LANE6 = 0xe55c1188cfd03434ULL;

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

static inline u64 mix13(u64 z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline u64 splitmix64_next(u64& x) {
    u64 z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline u8 byte_be(u64 x, int i) {
    return static_cast<u8>((x >> (56 - 8 * i)) & 0xffULL);
}

static inline State linear_layer(const State& s) {
    State out{};

    for (int j = 0; j < 8; ++j) {
        u8 ob[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        for (int k = 0; k < 8; ++k) {
            const u8 v = byte_be(s[k], j);
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
    State mixed{};
    for (int i = 0; i < 8; ++i) {
        mixed[i] = mix13(s[i]);
    }
    return linear_layer(mixed);
}

static inline State apply_rounds(State s, int rounds) {
    for (int r = 0; r < rounds; ++r) {
        s = round_F(s);
    }
    return s;
}

static inline State xor_state(const State& a, const State& b) {
    State d{};
    for (int i = 0; i < 8; ++i) {
        d[i] = a[i] ^ b[i];
    }
    return d;
}

static inline int parity64(u64 x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_parityll(x);
#else
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    x &= 0xf;
    return (0x6996 >> x) & 1;
#endif
}

static inline int distinguisher_bit(const State& delta) {
    const u64 selected =
        (delta[2] & FINAL_MASK_LANE2) ^
        (delta[6] & FINAL_MASK_LANE6);
    return parity64(selected);
}

static inline double log2_or_neg_inf(double x) {
    return x > 0.0 ? std::log2(x) : -999.0;
}

struct Result {
    u64 samples = 0;
    u64 pass_zero = 0;

    double probability() const {
        return static_cast<double>(pass_zero) / static_cast<double>(samples);
    }

    double epsilon() const {
        return probability() - 0.5;
    }

    double correlation() const {
        return std::fabs(2.0 * epsilon());
    }

    double z_score() const {
        const double sigma = std::sqrt(static_cast<double>(samples) / 4.0);
        return (static_cast<double>(pass_zero) - static_cast<double>(samples) / 2.0) / sigma;
    }
};

static Result measure(u64 n, int rounds) {
    Result result;
    result.samples = n;

    u64 sm = 0xDEADBEEF12345678ULL;

    for (u64 i = 0; i < n; ++i) {
        State x{};
        State xa{};

        for (int lane = 0; lane < 8; ++lane) {
            x[lane] = splitmix64_next(sm);
            xa[lane] = x[lane];
        }

        xa[0] ^= ALPHA_LANE0;

        const State y = apply_rounds(x, rounds);
        const State ya = apply_rounds(xa, rounds);
        const State delta = xor_state(y, ya);

        if (distinguisher_bit(delta) == 0) {
            ++result.pass_zero;
        }
    }

    return result;
}

static void print_result(const char* label, const Result& r) {
    const double p = r.probability();
    const double eps = r.epsilon();
    const double corr = r.correlation();
    const double pairs_corr = corr > 0.0 ? 1.0 / (corr * corr) : INFINITY;
    const double pairs_eps = eps != 0.0 ? 1.0 / (eps * eps) : INFINITY;

    std::printf("%s\n", label);
    std::printf("  samples:             %llu\n", static_cast<unsigned long long>(r.samples));
    std::printf("  mask parity == 0:    %llu\n", static_cast<unsigned long long>(r.pass_zero));
    std::printf("  probability:         %.10f\n", p);
    std::printf("  epsilon from 1/2:    %.10e\n", eps);
    std::printf("  correlation:         %.10e\n", corr);
    std::printf("  log2(correlation):   %.4f\n", log2_or_neg_inf(corr));
    std::printf("  z-score:             %.2f\n", r.z_score());
    std::printf("  approx 1/corr^2:     %.0f pairs\n", pairs_corr);
    std::printf("  approx 1/epsilon^2:  %.0f pairs\n", pairs_eps);
}

int main(int argc, char** argv) {
    int n_exp = 24;
    bool run_round3 = false;

    if (argc >= 2) {
        n_exp = std::atoi(argv[1]);
    }
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--round3") == 0) {
            run_round3 = true;
        }
    }

    if (n_exp < 10 || n_exp > 32) {
        std::fprintf(stderr, "N exponent must be in range 10..32\n");
        return 1;
    }

    const u64 n = 1ULL << n_exp;

    std::printf("Reduced SWIRL Round-2 Differential-Linear Distinguisher\n");
    std::printf("=======================================================\n\n");
    std::printf("Reduced model: Mix13 per lane, then byte/MDS linear layer.\n");
    std::printf("Round constants: omitted.\n");
    std::printf("Round keys:      omitted.\n\n");

    std::printf("Input difference alpha:\n");
    std::printf("  lane 0 = 0x%016llx\n", static_cast<unsigned long long>(ALPHA_LANE0));
    std::printf("  lanes 1..7 = 0\n\n");

    std::printf("Final-output parity mask:\n");
    std::printf("  lane 2 mask = 0x%016llx\n", static_cast<unsigned long long>(FINAL_MASK_LANE2));
    std::printf("  lane 6 mask = 0x%016llx\n\n", static_cast<unsigned long long>(FINAL_MASK_LANE6));

    std::printf("N = 2^%d = %llu input pairs\n\n", n_exp, static_cast<unsigned long long>(n));

    const Result r2 = measure(n, 2);
    print_result("Two reduced rounds:", r2);

    if (run_round3) {
        std::printf("\n");
        const Result r3 = measure(n, 3);
        print_result("Three reduced rounds, same mask carried forward:", r3);
    }

    return 0;
}

