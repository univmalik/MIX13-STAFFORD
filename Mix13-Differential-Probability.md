# Confirmed High-Probability XOR Differentials in Stafford Mix13

**Date:** 22 June 2026  
**Primitive:** Stafford Mix13 / SplitMix64 finalizer  
**Word size:** 64 bits  
**Analysis type:** Exact differential-candidate discovery followed by large fixed-candidate statistical validation  
**Status:** Statistically confirmed over $2^{32}$ validation inputs per tested input difference; not an exhaustive full-domain differential distribution

---

## 1. Executive Summary

This report documents high-probability XOR differentials in the canonical 64-bit Stafford Mix13 permutation.

Three input differences were evaluated:

```text
0x0008000000200000    two active input bits
0x0200000008000000    two active input bits
0x8000000200000008    three active input bits
```

For each input difference $\alpha$, the accompanying C++ program first discovered the most frequent exact output differences $\beta$ by generating and sorting $2^{25}$ complete 64-bit differential outputs. It then fixed the 16 leading candidates and validated them in four large campaigns of $2^{30}$ inputs each.

The total validation size for each tested $\alpha$ was:

```math
4\cdot2^{30}=2^{32}=4,294,967,296
```
input pairs.

The strongest confirmed differential in the experiment was:

```text
alpha = 0x8000000200000008
beta  = 0x53112150a62242a0
```

with:

```math
4,197,603
```
hits among:

```math
2^{32}
```
tested inputs. The measured probability was:

```math
\widehat{p}
=
\frac{4,197,603}{2^{32}}
=
9.773306082934\times10^{-4}
\approx2^{-9.998866}.
```
This is extremely close to:

```math
2^{-10},
```
or approximately one occurrence per 1,023 tested inputs.

The same input difference produced eight distinct exact output differences, each with probability close to $2^{-10}$. Their combined measured probability was:

```math
\frac{33,557,348}{2^{32}}
=
0.007813178934
\approx2^{-6.999875}.
```
Thus approximately one input pair in 128 produced an output difference belonging to this particular eight-member cluster.

The strongest two-bit input difference tested was:

```text
alpha = 0x0200000008000000
beta  = 0x88408901d0811201
```

with measured probability:

```math
4.203915596008\times10^{-4}
\approx2^{-11.215979}.
```
The results demonstrate that Mix13 contains large and highly structured XOR-differential biases. They do not prove that the identified three-bit differential is the global maximum over all nonzero 64-bit input differences.

---

## 2. Reference Implementation

The statistical experiment is implemented in:

```text
mix13_large_statistical_test.cpp
```

The program contains:

1. the canonical Mix13 implementation;
2. exact differential-candidate discovery by sorting;
3. fixed-candidate validation;
4. four source-diverse input campaigns;
5. multithreaded execution;
6. combined hit counting;
7. $\log_2$-probability reporting;
8. 95% Wilson confidence intervals;
9. automatic text-result export.

The canonical Mix13 function used by the program is:

```cpp
static inline uint64_t mix13(uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
```

All arithmetic uses `uint64_t`. Consequently, multiplication wraps modulo $2^{64}$, as required.

### 2.1 Compilation

Using GCC, MinGW-w64, or Clang:

```bash
g++ -O3 -march=native -std=c++17 -pthread \
    mix13_large_statistical_test.cpp \
    -o mix13_large_test
```

On Windows PowerShell:

```powershell
g++ -O3 -march=native -std=c++17 -pthread `
    mix13_large_statistical_test.cpp `
    -o mix13_large_test.exe
```

### 2.2 Reproduction command

The reported experiment used:

```powershell
.\mix13_large_test.exe `
  --discovery-exp 25 `
  --validate-exp 30 `
  --campaigns 4 `
  --threads 8 `
  --top 16 `
  --output mix13_large_test_results_large.txt
```

The parameters mean:

```text
Discovery inputs per alpha          = 2^25
Validation inputs per campaign      = 2^30
Validation campaigns                = 4
Total validation inputs per alpha   = 2^32
Retained output-difference candidates = 16
Worker threads                      = 8
```

---

## 3. Differential Definition

For a function:

```math
F:\{0,1\}^{64}\rightarrow\{0,1\}^{64},
```
a fixed input XOR difference $\alpha$, and a fixed output XOR difference $\beta$, define:

```math
D_\alpha F(x)
=
F(x)\oplus F(x\oplus\alpha).
```
The probability of the differential pair $(\alpha,\beta)$ is:

```math
\Pr_x[D_\alpha F(x)=\beta].
```
It is important to distinguish between:

```math
\Pr_x[D_\alpha F(x)=\beta]
```
for a particular beta and:

```math
\max_{\beta\ne0}
\Pr_x[D_\alpha F(x)=\beta],
```
which is the largest differential probability associated with a fixed alpha.

A valid differential pair is not necessarily the most probable pair for that input difference.

---

## 4. Ideal-Permutation Reference Scale

Mix13 is a permutation because:

1. each map $x\mapsto x\oplus(x\gg s)$ is invertible;
2. each multiplication constant is odd and therefore invertible modulo $2^{64}$;
3. a composition of invertible maps is invertible.

For a random 64-bit permutation and fixed nonzero $\alpha$ and $\beta$, the average fixed-pair probability is approximately:

```math
\frac{1}{2^{64}-1}
\approx2^{-64}.
```
This ideal reference is not a claim that every differential-table entry of a random permutation is identical. It is the natural fixed-pair comparison scale.

Against this reference:

- a probability near $2^{-17.84}$ is already enormously elevated;
- a probability near $2^{-11.22}$ is stronger still;
- a probability near $2^{-10}$ is an extreme component-level differential bias.

---

## 5. Structural Origin of the Input Differences

The first Mix13 operation is:

```math
L_{30}(x)=x\oplus(x\gg30).
```
Because this operation is linear with respect to XOR:

```math
D_\alpha L_{30}(x)
=
\alpha\oplus(\alpha\gg30).
```
This allows external input differences to be selected so that several active bits cancel before the first multiplication.

### 5.1 Original two-bit input difference

```text
alpha = 0x0008000000200000
```

Its active bits are 21 and 51. Since the positions are separated by 30:

```math
\alpha\oplus(\alpha\gg30)=2^{51}.
```
Numerically:

```text
alpha                 = 0x0008000000200000
alpha >> 30           = 0x0000000000200000
result                 = 0x0008000000000000
```

The two external active bits therefore reduce to one internal active bit at position 51.

### 5.2 Stronger two-bit input difference

```text
alpha = 0x0200000008000000
```

Its active bits are 27 and 57. The first XOR-shift gives:

```math
\alpha\oplus(\alpha\gg30)=2^{57}.
```
Numerically:

```text
alpha                 = 0x0200000008000000
alpha >> 30           = 0x0000000008000000
result                 = 0x0200000000000000
```

This leaves a one-bit difference at position 57.

Because the surviving difference is higher in the word than bit 51, fewer lower product positions remain available for carry propagation. The first multiplication therefore produces a more concentrated differential distribution.

### 5.3 Strongest presently confirmed three-bit input difference

```text
alpha = 0x8000000200000008
```

Its active bits are 3, 33, and 63. The first XOR-shift gives:

```text
alpha                 = 0x8000000200000008
alpha >> 30           = 0x0000000200000008
result                 = 0x8000000000000000
```

Therefore:

```math
\alpha\oplus(\alpha\gg30)=2^{63}.
```
The external three-bit input difference collapses exactly to a single most-significant-bit difference.

This alpha can also be derived using the inverse of the XOR-shift. For:

```math
L_{30}(x)=x\oplus(x\gg30),
```
the inverse difference transform is:

```math
L_{30}^{-1}(y)
=
y\oplus(y\gg30)\oplus(y\gg60).
```
Selecting:

```math
y=2^{63}
```
gives:

```math
L_{30}^{-1}(2^{63})
=
2^{63}\oplus2^{33}\oplus2^3,
```
which is exactly:

```text
0x8000000200000008
```

---

## 6. Exact Preservation of the MSB Difference

The three-bit alpha is particularly strong because its post-shift difference is $2^{63}$.

Let the first multiplication constant be the odd integer:

```math
C_1=\texttt{0xbf58476d1ce4e5b9}.
```
For any 64-bit value $u$, toggling bit 63 is equivalent modulo $2^{64}$ to adding $2^{63}$:

```math
u\oplus2^{63}
\equiv u+2^{63}\pmod{2^{64}}.
```
Multiplication gives:

```math
C_1(u+2^{63})
\equiv C_1u+C_1 2^{63}\pmod{2^{64}}.
```
Since $C_1$ is odd:

```math
C_1 2^{63}
\equiv2^{63}\pmod{2^{64}}.
```
Therefore, the two multiplication outputs differ exactly in the most significant bit:

```math
(C_1u)
\oplus
(C_1(u\oplus2^{63}))
=
2^{63}.
```
This exact property is special to the most-significant-bit difference. Multiplication by an odd integer is not generally linear with respect to XOR differences; for other differences, carry behavior depends on the input value.

After the following XOR-shift by 27, the difference becomes:

```math
2^{63}\oplus2^{36}.
```
Numerically:

```text
after first xor-shift = 0x8000000000000000
after first multiply  = 0x8000000000000000
after xor-shift by 27 = 0x8000001000000000
```

The second multiplication is where the main probabilistic carry-dependent branching begins.

---

## 7. Experimental Methodology

### 7.1 Phase 1: Exact candidate discovery

For each alpha, the program generated:

```math
2^{25}=33,554,432
```
values of:

```math
\beta_x
=
\operatorname{Mix13}(x)
\oplus
\operatorname{Mix13}(x\oplus\alpha).
```
Every complete 64-bit beta was stored in memory.

The vector was then sorted. Equal adjacent values were run-counted to obtain exact frequencies within the discovery dataset.

This method avoids the weaknesses of a small or lossy hash table:

- no beta is silently discarded;
- no probe limit is used;
- no bucket collision alters the count;
- the complete 64-bit difference is compared;
- the most frequent discovery candidates are identified exactly for that dataset.

The 16 most frequent beta values were retained.

Discovery memory usage was approximately:

```math
2^{25}\cdot8
=
2^{28}\text{ bytes}
=
256\text{ MiB}.
```
### 7.2 Phase 2: Fixed-candidate validation

After discovery, the candidate beta values were fixed.

The validation phase did not search for new candidates and did not change beta values based on validation counts. This separation reduces candidate-selection bias.

For every validation input, the program computed:

```math
\operatorname{Mix13}(x)
\oplus
\operatorname{Mix13}(x\oplus\alpha)
```
and incremented a counter only when the result matched one of the 16 preselected beta values.

### 7.3 Validation campaigns

Four campaigns were used:

```text
Campaign 1: fmix64(counter + offset)
Campaign 2: independent bijective permutation of counter + offset
Campaign 3: odd-affine permutation of counter
Campaign 4: fmix64(counter XOR mask)
```

Each campaign processed:

```math
2^{30}=1,073,741,824
```
inputs.

Within each campaign, the input-generation map is bijective over 64-bit words, so distinct counters produce distinct inputs.

The four campaigns use different deterministic input mappings. Their purpose is to test whether the measured rate remains stable under different sample enumerations.

Global disjointness between all four differently generated input sets is not formally asserted. The campaigns should therefore be viewed as source-diverse large deterministic samples rather than as a proof based on $2^{32}$ globally unique inputs.

### 7.4 Reported statistics

For each candidate beta, the program reports:

- hit count;
- empirical probability;
- base-2 logarithm of the probability;
- combined count across campaigns;
- 95% Wilson interval.

The interval is descriptive of the observed rate under a binomial-sampling interpretation. Because the sequences are deterministic, it should not be interpreted as a substitute for a full-domain exact count.

---

## 8. Results

### 8.1 Summary of strongest validated pairs

| Input difference $\alpha$ | Input weight | Strongest validated $\beta$ | Hits in $2^{32}$ | Probability | $\log_2(p)$ |
|---|---:|---|---:|---:|---:|
| `0x0008000000200000` | 2 | `0x40a2420804448411` | 18,260 | $4.251487553120\times10^{-6}$ | $-17.843601$ |
| `0x0200000008000000` | 2 | `0x88408901d0811201` | 1,805,568 | $4.203915596008\times10^{-4}$ | $-11.215979$ |
| `0x8000000200000008` | 3 | `0x53112150a62242a0` | 4,197,603 | $9.773306082934\times10^{-4}$ | $-9.998866$ |

### 8.2 Original-alpha result

For:

```text
alpha = 0x0008000000200000
```

the strongest beta among the 16 discovered and validated candidates was:

```text
beta = 0x40a2420804448411
```

The combined result was:

```text
hits       = 18,260
samples    = 4,294,967,296
p          = 4.251487553120e-06
log2(p)    = -17.843601
```

The 95% Wilson interval was:

```math
4.190268280500\times10^{-6}
\le p\le
4.313601227463\times10^{-6}.
```
The four campaign counts were:

```text
4,573
4,628
4,484
4,575
```

The agreement across campaigns supports the reproducibility of the result.

### 8.3 Strong two-bit result

For:

```text
alpha = 0x0200000008000000
```

the leading beta was:

```text
beta = 0x88408901d0811201
```

The combined result was:

```text
hits       = 1,805,568
samples    = 4,294,967,296
p          = 4.203915596008e-04
log2(p)    = -11.215979
```

The 95% Wilson interval was:

```math
4.197789450425\times10^{-4}
\le p\le
4.210050678165\times10^{-4}.
```
The four campaign counts were:

```text
451,411
451,366
451,129
451,662
```

A second beta was almost equally probable:

```text
beta       = 0x8840890050811202
hits       = 1,804,008
p          = 4.200283437967e-04
log2(p)    = -11.217226
```

The next fourteen validated candidates formed a cluster close to:

```math
2^{-12.214}.
```
This repeated probability structure indicates that the distribution is organized into carry-dependent families rather than containing one isolated anomalous beta.

### 8.4 Strong three-bit result

For:

```text
alpha = 0x8000000200000008
```

the highest combined count belonged to:

```text
beta = 0x53112150a62242a0
```

The combined result was:

```text
hits       = 4,197,603
samples    = 4,294,967,296
p          = 9.773306082934e-04
log2(p)    = -9.998866
```

The 95% Wilson interval was:

```math
9.763965603782\times10^{-4}
\le p\le
9.782655488696\times10^{-4}.
```
The four campaign counts were:

```text
1,050,185
1,049,539
1,049,092
1,048,787
```

The result is extremely close to $2^{-10}$.

---

## 9. Eight-Member Differential Cluster

The leading eight exact output differences for:

```text
alpha = 0x8000000200000008
```

were:

| Rank | Beta | Hits | Probability | $\log_2(p)$ |
|---:|---|---:|---:|---:|
| 1 | `0x53112150a62242a0` | 4,197,603 | $9.773306082934\times10^{-4}$ | $-9.998866$ |
| 2 | `0xd3112151a62242a0` | 4,196,674 | $9.771143086255\times10^{-4}$ | $-9.999185$ |
| 3 | `0xb3112151662242a0` | 4,196,448 | $9.770616889000\times10^{-4}$ | $-9.999263$ |
| 4 | `0xcd1121519a2242a0` | 4,195,099 | $9.767476003617\times10^{-4}$ | $-9.999727$ |
| 5 | `0x4d1121509a2242a0` | 4,194,033 | $9.764994028956\times10^{-4}$ | $-10.000093$ |
| 6 | `0xd5112151aa2242a0` | 4,193,164 | $9.762970730662\times10^{-4}$ | $-10.000392$ |
| 7 | `0xb51121516a2242a0` | 4,192,244 | $9.760828688741\times10^{-4}$ | $-10.000709$ |
| 8 | `0x55112150aa2242a0` | 4,192,083 | $9.760453831404\times10^{-4}$ | $-10.000764$ |

The combined hit count was:

```math
33,557,348.
```
The combined probability was:

```math
\frac{33,557,348}{2^{32}}
=
0.007813178934
\approx2^{-6.999875}.
```
This is effectively:

```math
2^{-7}.
```
Consequently, approximately one input pair in 128 produced one of these eight exact 64-bit output differences.

A second cluster of candidates appeared near $2^{-11}$, showing another approximate factor-of-two probability tier.

---

## 10. Relative Strength

Using the strongest validated beta for each tested alpha:

### 10.1 Three-bit versus strong two-bit alpha

```math
\frac{9.773306082934\times10^{-4}}
     {4.203915596008\times10^{-4}}
\approx2.325.
```
The leading three-bit differential is approximately 2.3 times more probable than the leading two-bit differential.

### 10.2 Strong two-bit versus original alpha

```math
\frac{4.203915596008\times10^{-4}}
     {4.251487553120\times10^{-6}}
\approx98.88.
```
The stronger two-bit differential is approximately 99 times more probable than the strongest validated candidate for the original alpha.

### 10.3 Three-bit versus original alpha

```math
\frac{9.773306082934\times10^{-4}}
     {4.251487553120\times10^{-6}}
\approx229.87.
```
The leading three-bit differential is approximately 230 times more probable than the strongest validated candidate for the original alpha.

---

## 11. Why Avalanche Testing Does Not Detect This

A standard avalanche test typically measures whether changing one or more input bits causes approximately half of the output bits to change.

This is not the same as measuring exact differential probabilities.

Two output differences may both have Hamming weight close to 32 while representing completely different 64-bit values. Avalanche testing merges all such patterns into summary statistics.

Differential analysis instead asks:

```math
\text{Does the same exact 64-bit output difference recur unusually often?}
```
Mix13 can produce broad bit diffusion and still exhibit highly concentrated exact XOR-difference distributions.

The present results therefore do not contradict good avalanche measurements. They show that avalanche behavior is insufficient as a test of differential uniformity.

---

## 12. Interpretation of the Finding

The appropriate description is:

> Confirmed high-probability XOR differentials and differential clustering in Stafford Mix13.

The experiment establishes that:

1. exact output differences recur at probabilities vastly above the ideal fixed-pair scale;
2. the bias is reproducible across several large input enumerations;
3. the bias follows structural cancellation through the first XOR-shift;
4. the strongest tested alpha forces an exact MSB difference through the first multiplication;
5. multiple preferred beta values form regular probability clusters;
6. Mix13 should not be modeled as an ideal cryptographic 64-bit permutation under XOR differentials.

Because Mix13 is bijective, a nonzero alpha cannot produce beta zero. Therefore, these findings are not ordinary collisions inside Mix13.

The result is a differential distinguishability weakness in the permutation.

---

## 13. Scope and Limitations

### 13.1 No global maximum proof

The strongest presently confirmed pair is:

```text
0x8000000200000008 -> 0x53112150a62242a0
```

near $2^{-10}$.

However, the program did not search all:

```math
2^{64}-1
```
nonzero input differences. It therefore does not prove that $2^{-10}$ is the global maximum differential probability of Mix13.

The defensible wording is:

> strongest differential presently confirmed among the tested candidates.

### 13.2 Candidate limitation

Only the 16 most frequent beta values from the $2^{25}$-sample discovery phase were validated.

A beta omitted from the discovery top 16 was not measured during large validation. For the highly concentrated two-bit and three-bit distributions, the leading candidates occurred thousands or tens of thousands of times in discovery, making their identification stable.

For weaker input differences, a larger `--top` value or larger discovery set may be appropriate.

### 13.3 Statistical rather than exact full-domain counts

The experiment measures large deterministic samples. It does not enumerate all $2^{64}$ possible x values.

The displayed probabilities are therefore empirical estimates, not exact complete differential-distribution-table entries.

### 13.4 Campaign overlap

Each individual campaign contains distinct inputs because its input generator is bijective on distinct counters.

The four campaign sets are not proven mutually disjoint. This does not explain the observed high probabilities, but it means that the combined count should be described as $2^{32}$ tested inputs rather than $2^{32}$ proven globally unique inputs.

### 13.5 Confidence-interval interpretation

The Wilson intervals quantify uncertainty under a binomial-sampling model. The test sequences are deterministic and permutation-generated.

The intervals are useful descriptions of rate stability but are not a formal proof about the complete 64-bit domain.

---

## 14. Recommended Additional Mix13 Tests

The following work would strengthen the analysis further.

### 14.1 Repeat with larger discovery

For example:

```powershell
.\mix13_large_test.exe `
  --alpha 8000000200000008 `
  --discovery-exp 27 `
  --validate-exp 30 `
  --campaigns 4 `
  --threads 8 `
  --top 64
```

This uses approximately 1 GiB for discovery and retains more candidate beta values.

### 14.2 Validate selected fixed betas directly

A specialized fixed-pair counter can test only:

```text
alpha = 0x8000000200000008
beta  = 0x53112150a62242a0
```

over larger input ranges without storing discovery values.

### 14.3 Search all low-weight alphas

A systematic search can include:

- all 64 one-bit alphas;
- all 2,016 two-bit alphas;
- all 41,664 three-bit alphas;
- structured four-bit and higher-weight families;
- all inverse images of selected low-weight internal differences under the first XOR-shift.

Candidate discovery must use a non-lossy counting method or a reliable heavy-hitter algorithm followed by independent validation.

### 14.4 Study intermediate differential distributions

For each strong alpha, record exact differences:

1. after `x ^= x >> 30`;
2. after the first multiplication;
3. after `x ^= x >> 27`;
4. after the second multiplication;
5. after `x ^= x >> 31`.

This would identify the exact carry branches responsible for the $2^{-10}$, $2^{-11}$, and $2^{-12}$ probability tiers.

### 14.5 Derive exact probabilities analytically

The near-power-of-two values strongly suggest a finite set of carry conditions.

An analytical model could classify the relevant carry patterns and possibly derive exact differential probabilities without full-domain enumeration.

---

## 15. Reproducibility Checklist

A reproduction should verify all of the following:

- [ ] the Mix13 constants exactly match the canonical values;
- [ ] all arithmetic uses unsigned 64-bit wraparound;
- [ ] the first XOR difference is computed as `mix13(x) ^ mix13(x ^ alpha)`;
- [ ] discovery counts complete 64-bit beta values;
- [ ] discovery and validation use separate sequences;
- [ ] beta candidates are fixed before validation;
- [ ] validation does not reselect the winner;
- [ ] per-campaign counts are retained;
- [ ] total sample counts use 64-bit integers;
- [ ] probabilities are calculated as `hits / samples`;
- [ ] the result file records command-line parameters;
- [ ] compiler optimization does not alter unsigned wraparound semantics.

---

## 16. Final Assessment

The large statistical experiment confirms that Stafford Mix13 has substantial exact XOR-differential biases.

The strongest presently confirmed pair is:

```text
Input difference:
0x8000000200000008

Output difference:
0x53112150a62242a0
```

with measured probability:

```math
9.773306082934\times10^{-4}
\approx2^{-9.998866}.
```
Seven additional exact output differences for the same input difference also occur individually near $2^{-10}$. Together, the leading eight-member cluster occurs with probability approximately:

```math
2^{-7}.
```
The strongest tested two-bit input difference is:

```text
0x0200000008000000
```

with a leading exact differential near:

```math
2^{-11.216}.
```
These findings are far above the approximately $2^{-64}$ fixed-pair reference scale of an ideal 64-bit permutation. They establish a clear differential distinguishability weakness in Mix13.

The analysis is statistically compelling but should not be overstated as a proof of the global maximum differential probability. A complete global result would require either an exact analytical derivation or a substantially broader exhaustive strategy.

---

## Appendix A — Files

### A.1 Statistical test source

```text
mix13_large_statistical_test.cpp
```

Purpose:

- discover exact beta candidates;
- validate fixed candidates over large campaigns;
- calculate probabilities and confidence intervals;
- save complete textual output.

### A.2 Result log

```text
mix13_large_test_results_large.txt
```

The documented run used:

```text
discovery samples per alpha = 2^25
validation samples/campaign = 2^30
campaigns                    = 4
validation total/alpha       = 2^32
threads                      = 8
top beta candidates          = 16
```

---

## Appendix B — Primary Confirmed Values

```text
Original alpha:
alpha = 0x0008000000200000
beta  = 0x40a2420804448411
p     = 4.251487553120e-06
log2  = -17.843601

Strong two-bit alpha:
alpha = 0x0200000008000000
beta  = 0x88408901d0811201
p     = 4.203915596008e-04
log2  = -11.215979

Strong three-bit alpha:
alpha = 0x8000000200000008
beta  = 0x53112150a62242a0
p     = 9.773306082934e-04
log2  = -9.998866
```
