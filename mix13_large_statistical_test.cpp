/*
 * Large Statistical Differential Test for Stafford Mix13 / SplitMix64 finalizer
 *
 * Purpose
 * -------
 * 1. Discover the most frequent exact output differences beta for each alpha
 *    using exact sorting (no lossy hash table and no dropped entries).
 * 2. Validate those fixed beta candidates over much larger campaigns.
 * 3. Report counts, probabilities, log2 probabilities, and 95% Wilson intervals.
 *
 * Default alphas
 * --------------
 *   original report: 0x0008000000200000
 *   strongest 2-bit candidate found so far: 0x0200000008000000
 *   stronger 3-bit candidate found so far: 0x8000000200000008
 *
 * Build (GCC / MinGW / Clang):
 *   g++ -O3 -march=native -std=c++17 -pthread \
 *       mix13_large_statistical_test.cpp -o mix13_large_test
 *
 * Typical runs:
 *   ./mix13_large_test
 *
 *   # Larger run: discovery=2^25, validation=2^30 per campaign,
 *   # 4 campaigns, 8 worker threads:
 *   ./mix13_large_test --discovery-exp 25 --validate-exp 30 \
 *       --campaigns 4 --threads 8 --top 16
 *
 *   # Test one custom alpha:
 *   ./mix13_large_test --alpha 8000000200000008 \
 *       --discovery-exp 25 --validate-exp 31 --campaigns 4 --threads 8
 *
 * Notes
 * -----
 * - Multiplication is unsigned uint64_t and therefore wraps modulo 2^64.
 * - Discovery stores 2^discovery_exp uint64_t values. Memory is approximately:
 *       2^(discovery_exp + 3) bytes
 *   For example, exp=24 uses 128 MiB; exp=25 uses 256 MiB.
 * - Validation does not store samples and therefore uses little memory.
 * - Candidate beta values are fixed before validation, avoiding selection bias
 *   in the large validation campaigns.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// Canonical Stafford Mix13
// -----------------------------------------------------------------------------
static inline uint64_t mix13(uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

// Independent bijective 64-bit permutation used to produce distinct test inputs
// from distinct counters. This is MurmurHash3 fmix64; all multipliers are odd
// and all xor-shifts are invertible, so it is a permutation of 64-bit words.
static inline uint64_t fmix64(uint64_t x) noexcept {
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    return x;
}

static inline uint64_t rotl64(uint64_t x, unsigned r) noexcept {
    return (x << r) | (x >> (64U - r));
}

// A second bijective input permutation, used for source-diverse validation.
// Composition: xor-shifts, rotations, and odd multiplications.
static inline uint64_t permute_b(uint64_t x) noexcept {
    x ^= x >> 29;
    x *= UINT64_C(0x9fb21c651e98df25);
    x = rotl64(x, 31);
    x ^= x >> 26;
    x *= UINT64_C(0xd6e8feb86659fd93);
    x ^= x >> 32;
    return x;
}

static inline int popcount64(uint64_t x) noexcept {
#if defined(__GNUG__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
#endif
}

static std::string hex64(uint64_t x) {
    std::ostringstream ss;
    ss << "0x" << std::hex << std::setfill('0') << std::setw(16) << x;
    return ss.str();
}

static uint64_t parse_u64_hex(std::string s) {
    if (s.rfind("0x", 0) == 0 || s.rfind("0X", 0) == 0) s.erase(0, 2);
    if (s.empty() || s.size() > 16) throw std::invalid_argument("invalid 64-bit hex value");
    size_t used = 0;
    uint64_t v = std::stoull(s, &used, 16);
    if (used != s.size()) throw std::invalid_argument("invalid hex character in alpha");
    return v;
}

struct Config {
    unsigned discovery_exp = 24;   // 2^24 = 16,777,216; 128 MiB
    unsigned validate_exp = 28;    // per campaign
    unsigned campaigns = 4;
    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    unsigned top_k = 16;
    std::string output_file = "mix13_large_test_results.txt";
    bool custom_alpha = false;
    uint64_t alpha = 0;
};

static void print_usage(const char* exe) {
    std::cout
        << "Usage: " << exe << " [options]\n\n"
        << "Options:\n"
        << "  --alpha HEX           Test only one custom 64-bit input difference\n"
        << "  --discovery-exp N     Discovery samples = 2^N (default 24)\n"
        << "  --validate-exp N      Samples per validation campaign = 2^N (default 28)\n"
        << "  --campaigns N         Number of validation campaigns (default 4)\n"
        << "  --threads N           Worker threads (default hardware concurrency)\n"
        << "  --top N               Number of candidate beta values retained (default 16)\n"
        << "  --output FILE         Save complete text output to FILE\n"
        << "  --help                Show this message\n\n"
        << "Memory for discovery is approximately 2^(discovery-exp + 3) bytes.\n";
}

static unsigned parse_unsigned(const char* s, const char* name) {
    char* end = nullptr;
    unsigned long v = std::strtoul(s, &end, 10);
    if (!s[0] || !end || *end != '\0' || v > std::numeric_limits<unsigned>::max()) {
        throw std::invalid_argument(std::string("invalid value for ") + name);
    }
    return static_cast<unsigned>(v);
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value after ") + name);
            return argv[++i];
        };

        if (arg == "--alpha") {
            cfg.alpha = parse_u64_hex(require_value("--alpha"));
            cfg.custom_alpha = true;
        } else if (arg == "--discovery-exp") {
            cfg.discovery_exp = parse_unsigned(require_value("--discovery-exp"), "--discovery-exp");
        } else if (arg == "--validate-exp") {
            cfg.validate_exp = parse_unsigned(require_value("--validate-exp"), "--validate-exp");
        } else if (arg == "--campaigns") {
            cfg.campaigns = parse_unsigned(require_value("--campaigns"), "--campaigns");
        } else if (arg == "--threads") {
            cfg.threads = parse_unsigned(require_value("--threads"), "--threads");
        } else if (arg == "--top") {
            cfg.top_k = parse_unsigned(require_value("--top"), "--top");
        } else if (arg == "--output") {
            cfg.output_file = require_value("--output");
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    if (cfg.discovery_exp < 10 || cfg.discovery_exp > 32)
        throw std::invalid_argument("--discovery-exp must be between 10 and 32");
    if (cfg.validate_exp < 10 || cfg.validate_exp > 40)
        throw std::invalid_argument("--validate-exp must be between 10 and 40");
    if (cfg.campaigns == 0 || cfg.campaigns > 64)
        throw std::invalid_argument("--campaigns must be between 1 and 64");
    if (cfg.threads == 0 || cfg.threads > 256)
        throw std::invalid_argument("--threads must be between 1 and 256");
    if (cfg.top_k == 0 || cfg.top_k > 64)
        throw std::invalid_argument("--top must be between 1 and 64");
    return cfg;
}

// Mirror output to console and file.
class Tee {
public:
    explicit Tee(const std::string& path) : file_(path, std::ios::out | std::ios::trunc) {
        if (!file_) throw std::runtime_error("cannot open output file: " + path);
    }

    template <typename T>
    Tee& operator<<(const T& value) {
        std::cout << value;
        file_ << value;
        return *this;
    }

    Tee& operator<<(std::ostream& (*manip)(std::ostream&)) {
        manip(std::cout);
        manip(file_);
        return *this;
    }

private:
    std::ofstream file_;
};

struct Candidate {
    uint64_t beta = 0;
    uint64_t discovery_count = 0;
};

// Generate exact discovery beta values into a vector, then sort and run-count.
static std::vector<Candidate> discover_candidates(
    uint64_t alpha,
    uint64_t samples,
    unsigned top_k,
    unsigned threads,
    Tee& out)
{
    out << "\n[Discovery] alpha=" << hex64(alpha)
        << "  samples=" << samples
        << "  memory=" << std::fixed << std::setprecision(1)
        << (samples * sizeof(uint64_t) / 1048576.0) << " MiB\n";

    std::vector<uint64_t> diffs;
    try {
        diffs.resize(static_cast<size_t>(samples));
    } catch (const std::bad_alloc&) {
        throw std::runtime_error("not enough memory for discovery vector; reduce --discovery-exp");
    }

    const auto start_gen = std::chrono::steady_clock::now();
    const unsigned workers = static_cast<unsigned>(std::min<uint64_t>(threads, samples));
    std::vector<std::thread> pool;
    pool.reserve(workers);

    for (unsigned t = 0; t < workers; ++t) {
        uint64_t begin = (samples * t) / workers;
        uint64_t end   = (samples * (t + 1)) / workers;
        pool.emplace_back([&, begin, end]() {
            // Discovery source is a disjoint-counter fmix64 permutation.
            constexpr uint64_t OFFSET = UINT64_C(0x13579bdf2468ace0);
            for (uint64_t i = begin; i < end; ++i) {
                const uint64_t x = fmix64(i + OFFSET);
                diffs[static_cast<size_t>(i)] = mix13(x) ^ mix13(x ^ alpha);
            }
        });
    }
    for (auto& th : pool) th.join();

    const auto end_gen = std::chrono::steady_clock::now();
    const double gen_sec = std::chrono::duration<double>(end_gen - start_gen).count();
    out << "  generated in " << std::setprecision(3) << gen_sec << " s ("
        << std::setprecision(2) << (samples / gen_sec / 1e6) << " M pairs/s)\n";

    out << "  sorting exact 64-bit differences...\n";
    const auto start_sort = std::chrono::steady_clock::now();
    std::sort(diffs.begin(), diffs.end());
    const auto end_sort = std::chrono::steady_clock::now();
    const double sort_sec = std::chrono::duration<double>(end_sort - start_sort).count();
    out << "  sorted in " << std::setprecision(3) << sort_sec << " s\n";

    std::vector<Candidate> top;
    top.reserve(top_k);

    auto consider = [&](uint64_t beta, uint64_t count) {
        Candidate c{beta, count};
        auto pos = std::lower_bound(top.begin(), top.end(), c,
            [](const Candidate& a, const Candidate& b) {
                if (a.discovery_count != b.discovery_count)
                    return a.discovery_count > b.discovery_count;
                return a.beta < b.beta;
            });
        top.insert(pos, c);
        if (top.size() > top_k) top.pop_back();
    };

    if (!diffs.empty()) {
        uint64_t current = diffs[0];
        uint64_t count = 1;
        for (size_t i = 1; i < diffs.size(); ++i) {
            if (diffs[i] == current) {
                ++count;
            } else {
                consider(current, count);
                current = diffs[i];
                count = 1;
            }
        }
        consider(current, count);
    }

    out << "  top exact beta candidates:\n";
    for (size_t i = 0; i < top.size(); ++i) {
        const double p = static_cast<double>(top[i].discovery_count) / static_cast<double>(samples);
        out << "    #" << std::setw(2) << (i + 1)
            << " beta=" << hex64(top[i].beta)
            << " hw=" << std::setw(2) << popcount64(top[i].beta)
            << " count=" << std::setw(10) << top[i].discovery_count
            << " p=" << std::scientific << std::setprecision(9) << p
            << " log2(p)=" << std::fixed << std::setprecision(5) << std::log2(p)
            << "\n";
    }
    return top;
}

// Fixed-size lookup table for up to 64 beta candidates. No allocation in the hot loop.
class CandidateTable {
public:
    explicit CandidateTable(const std::vector<Candidate>& candidates) {
        size_ = 1;
        while (size_ < candidates.size() * 4) size_ <<= 1;
        if (size_ < 16) size_ = 16;
        keys_.assign(size_, 0);
        indices_.assign(size_, -1);
        mask_ = size_ - 1;

        for (size_t i = 0; i < candidates.size(); ++i) {
            size_t slot = hash(candidates[i].beta) & mask_;
            while (indices_[slot] != -1) slot = (slot + 1) & mask_;
            keys_[slot] = candidates[i].beta;
            indices_[slot] = static_cast<int>(i);
        }
    }

    inline int find(uint64_t key) const noexcept {
        size_t slot = hash(key) & mask_;
        for (;;) {
            const int idx = indices_[slot];
            if (idx == -1) return -1;
            if (keys_[slot] == key) return idx;
            slot = (slot + 1) & mask_;
        }
    }

private:
    static inline uint64_t hash(uint64_t x) noexcept {
        x ^= x >> 30;
        x *= UINT64_C(0xbf58476d1ce4e5b9);
        x ^= x >> 27;
        x *= UINT64_C(0x94d049bb133111eb);
        x ^= x >> 31;
        return x;
    }

    size_t size_ = 0;
    size_t mask_ = 0;
    std::vector<uint64_t> keys_;
    std::vector<int> indices_;
};

struct ValidationResult {
    std::vector<uint64_t> counts;
    double seconds = 0.0;
};

enum class Source : unsigned {
    FmixDistinct = 0,
    PermuteBDistinct = 1,
    AffineDistinct = 2,
    FmixAlternate = 3
};

static const char* source_name(Source s) {
    switch (s) {
        case Source::FmixDistinct:     return "fmix64(counter+offset)";
        case Source::PermuteBDistinct: return "permute_b(counter+offset)";
        case Source::AffineDistinct:   return "odd-affine counter";
        case Source::FmixAlternate:    return "fmix64(counter^mask)";
    }
    return "unknown";
}

static inline uint64_t make_input(Source source, uint64_t index, uint64_t campaign) noexcept {
    switch (source) {
        case Source::FmixDistinct:
            return fmix64(index + UINT64_C(0x243f6a8885a308d3)
                         + campaign * UINT64_C(0x9e3779b97f4a7c15));
        case Source::PermuteBDistinct:
            return permute_b(index + UINT64_C(0x13198a2e03707344)
                            + campaign * UINT64_C(0xd1b54a32d192ed03));
        case Source::AffineDistinct:
            return (index + campaign * UINT64_C(0x100000001b3))
                 * UINT64_C(0xd1342543de82ef95)
                 + UINT64_C(0xa4093822299f31d0);
        case Source::FmixAlternate:
            return fmix64((index ^ UINT64_C(0x082efa98ec4e6c89))
                         + campaign * UINT64_C(0x94d049bb133111eb));
    }
    return index;
}

static ValidationResult validate_campaign(
    uint64_t alpha,
    const std::vector<Candidate>& candidates,
    uint64_t samples,
    unsigned threads,
    Source source,
    uint64_t campaign_index)
{
    CandidateTable table(candidates);
    const unsigned workers = static_cast<unsigned>(std::min<uint64_t>(threads, samples));
    std::vector<std::vector<uint64_t>> local(
        workers, std::vector<uint64_t>(candidates.size(), 0));
    std::vector<std::thread> pool;
    pool.reserve(workers);

    const auto start = std::chrono::steady_clock::now();
    for (unsigned t = 0; t < workers; ++t) {
        const uint64_t begin = (samples * t) / workers;
        const uint64_t end   = (samples * (t + 1)) / workers;
        pool.emplace_back([&, t, begin, end]() {
            auto& counts = local[t];
            for (uint64_t i = begin; i < end; ++i) {
                const uint64_t x = make_input(source, i, campaign_index);
                const uint64_t beta = mix13(x) ^ mix13(x ^ alpha);
                const int idx = table.find(beta);
                if (idx >= 0) ++counts[static_cast<size_t>(idx)];
            }
        });
    }
    for (auto& th : pool) th.join();
    const auto stop = std::chrono::steady_clock::now();

    ValidationResult result;
    result.counts.assign(candidates.size(), 0);
    for (const auto& v : local) {
        for (size_t i = 0; i < v.size(); ++i) result.counts[i] += v[i];
    }
    result.seconds = std::chrono::duration<double>(stop - start).count();
    return result;
}

struct Interval {
    double low;
    double high;
};

static Interval wilson95(uint64_t hits, uint64_t n) {
    if (n == 0) return {0.0, 0.0};
    constexpr double z = 1.959963984540054;
    const double nn = static_cast<double>(n);
    const double p = static_cast<double>(hits) / nn;
    const double z2 = z * z;
    const double denom = 1.0 + z2 / nn;
    const double center = (p + z2 / (2.0 * nn)) / denom;
    const double radius = z * std::sqrt((p * (1.0 - p) / nn) + z2 / (4.0 * nn * nn)) / denom;
    return {std::max(0.0, center - radius), std::min(1.0, center + radius)};
}

static void test_alpha(const Config& cfg, uint64_t alpha, const std::string& label, Tee& out) {
    out << "\n============================================================\n";
    out << "ALPHA: " << hex64(alpha) << "  label=" << label
        << "  Hamming weight=" << popcount64(alpha) << "\n";
    out << "After first xor-shift: alpha ^ (alpha >> 30) = "
        << hex64(alpha ^ (alpha >> 30))
        << "  HW=" << popcount64(alpha ^ (alpha >> 30)) << "\n";
    out << "============================================================\n";

    const uint64_t discovery_samples = UINT64_C(1) << cfg.discovery_exp;
    const uint64_t validation_samples = UINT64_C(1) << cfg.validate_exp;

    auto candidates = discover_candidates(
        alpha, discovery_samples, cfg.top_k, cfg.threads, out);

    std::vector<uint64_t> totals(candidates.size(), 0);
    uint64_t total_samples = 0;

    out << "\n[Validation] fixed candidates; no candidate reselection\n";
    for (unsigned c = 0; c < cfg.campaigns; ++c) {
        const Source source = static_cast<Source>(c % 4);
        auto r = validate_campaign(
            alpha, candidates, validation_samples, cfg.threads, source, c + 1);

        total_samples += validation_samples;
        for (size_t i = 0; i < totals.size(); ++i) totals[i] += r.counts[i];

        out << "\n  Campaign " << (c + 1)
            << "  source=" << source_name(source)
            << "  N=" << validation_samples
            << "  time=" << std::fixed << std::setprecision(3) << r.seconds << " s"
            << "  rate=" << std::setprecision(2)
            << (validation_samples / r.seconds / 1e6) << " M pairs/s\n";

        for (size_t i = 0; i < candidates.size(); ++i) {
            const uint64_t h = r.counts[i];
            const double p = static_cast<double>(h) / static_cast<double>(validation_samples);
            out << "    #" << std::setw(2) << (i + 1)
                << " beta=" << hex64(candidates[i].beta)
                << " hits=" << std::setw(10) << h
                << " p=" << std::scientific << std::setprecision(9) << p;
            if (h > 0) out << " log2(p)=" << std::fixed << std::setprecision(5) << std::log2(p);
            else out << " log2(p)=-inf";
            out << "\n";
        }
    }

    out << "\n[Combined validation] N=" << total_samples << "\n";
    for (size_t i = 0; i < candidates.size(); ++i) {
        const uint64_t h = totals[i];
        const double p = static_cast<double>(h) / static_cast<double>(total_samples);
        const Interval ci = wilson95(h, total_samples);
        out << "  #" << std::setw(2) << (i + 1)
            << " beta=" << hex64(candidates[i].beta)
            << " hw=" << std::setw(2) << popcount64(candidates[i].beta)
            << " hits=" << std::setw(12) << h
            << " p=" << std::scientific << std::setprecision(12) << p;
        if (h > 0) out << " log2(p)=" << std::fixed << std::setprecision(6) << std::log2(p);
        else out << " log2(p)=-inf";
        out << "  95% CI=[" << std::scientific << std::setprecision(12)
            << ci.low << ", " << ci.high << "]\n";
    }

    // Sort by validation total for a clear final ranking.
    std::vector<size_t> order(candidates.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (totals[a] != totals[b]) return totals[a] > totals[b];
        return candidates[a].beta < candidates[b].beta;
    });

    out << "\n[Final ranking by combined validation]\n";
    for (size_t rank = 0; rank < order.size(); ++rank) {
        const size_t i = order[rank];
        const double p = static_cast<double>(totals[i]) / static_cast<double>(total_samples);
        out << "  rank " << std::setw(2) << (rank + 1)
            << " beta=" << hex64(candidates[i].beta)
            << " hits=" << std::setw(12) << totals[i]
            << " p=" << std::scientific << std::setprecision(12) << p;
        if (totals[i] > 0)
            out << " log2(p)=" << std::fixed << std::setprecision(6) << std::log2(p);
        out << "\n";
    }
}

int main(int argc, char** argv) {
    try {
        const Config cfg = parse_args(argc, argv);
        Tee out(cfg.output_file);

        out << "Mix13 Large Statistical Differential Test\n";
        out << "=========================================\n";
        out << "discovery samples per alpha = 2^" << cfg.discovery_exp
            << " = " << (UINT64_C(1) << cfg.discovery_exp) << "\n";
        out << "validation samples/campaign = 2^" << cfg.validate_exp
            << " = " << (UINT64_C(1) << cfg.validate_exp) << "\n";
        out << "campaigns = " << cfg.campaigns << "\n";
        out << "validation total/alpha = "
            << ((UINT64_C(1) << cfg.validate_exp) * cfg.campaigns) << "\n";
        out << "threads = " << cfg.threads << "\n";
        out << "top beta candidates = " << cfg.top_k << "\n";
        out << "output file = " << cfg.output_file << "\n";
        out << "mix13(0) = " << hex64(mix13(0)) << " (must be 0)\n";
        out << "mix13(1) = " << hex64(mix13(1)) << "\n";

        const auto overall_start = std::chrono::steady_clock::now();

        if (cfg.custom_alpha) {
            test_alpha(cfg, cfg.alpha, "custom", out);
        } else {
            test_alpha(cfg, UINT64_C(0x0008000000200000), "original report alpha", out);
            test_alpha(cfg, UINT64_C(0x0200000008000000), "strong 2-bit candidate", out);
            test_alpha(cfg, UINT64_C(0x8000000200000008), "stronger 3-bit candidate", out);
        }

        const auto overall_stop = std::chrono::steady_clock::now();
        out << "\nTotal elapsed time: " << std::fixed << std::setprecision(3)
            << std::chrono::duration<double>(overall_stop - overall_start).count()
            << " seconds\n";
        out << "Results saved to: " << cfg.output_file << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }
}
