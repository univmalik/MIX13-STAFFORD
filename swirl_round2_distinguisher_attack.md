# Reduced SWIRL Round-2 Differential-Linear Distinguisher

**Date:** 29 June 2026  
**Target:** Reduced two-round SWIRL512 data path  
**Core component:** Stafford Mix13 / SplitMix64 finalizer  
**Attack type:** Differential-linear distinguisher  
**Code:** `swirl_round2_distinguisher_attack.cpp`

---

## 1. Summary

This document describes a reduced-round distinguisher for two rounds of the
SWIRL data path. It does not find a collision, preimage, or full-round break.
It finds a statistical fingerprint: one carefully selected parity of the
two-round output difference is biased.

For a random 512-bit permutation, the tested parity should be zero with
probability:

```text
0.5
```

For the reduced two-round SWIRL path, the measured probability is:

```text
0.51826613
```

That means the selected parity bit is zero about:

```text
51.8266% of the time
```

instead of the ideal:

```text
50% of the time
```

The bias is small, but it is statistically very strong when repeated over many
chosen input pairs.

---

## 2. Reduced Model

The tested reduced SWIRL round is:

```text
round_F(state):
  for each 64-bit lane:
    lane = Mix13(lane)

  state = byte/MDS linear diffusion layer
```

The code deliberately omits:

```text
round constants
round keys
```

For XOR-difference experiments, fixed XOR round constants cancel. Fixed XOR
round keys also preserve the paired input difference. Under uniform input
sampling, they translate both paired states by the same value and do not change
the measured differential distribution.

So this is a data-path distinguisher for the reduced round function, not a claim
about key recovery.

---

## 3. Mix13

SWIRL uses the Stafford Mix13 permutation as its 64-bit lane mixer:

```cpp
z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9;
z = (z ^ (z >> 27)) * 0x94D049BB133111EB;
z = z ^ (z >> 31);
```

Mix13 has strong local XOR-differential structure. The strongest standalone
effects include:

```text
deterministic linear invariants for selected shift-30 alphas
high-probability exact differential clusters
conditional carry-free events
```

SWIRL's diffusion layer largely destroys exact Mix13 differential trails after
one full round. However, a weaker two-round projected signal remains.

---

## 4. Input Difference

The distinguisher uses the following input difference in lane 0:

```text
alpha = 0x0200000008000000
```

All other lanes are inactive:

```text
alpha_state:
  lane 0 = 0x0200000008000000
  lane 1 = 0
  lane 2 = 0
  lane 3 = 0
  lane 4 = 0
  lane 5 = 0
  lane 6 = 0
  lane 7 = 0
```

This is a shift-30 Mix13-aligned two-bit difference:

```text
bits 27 and 57 are active
```

The first Mix13 xor-shift can collapse this kind of difference into a more
structured internal difference. SWIRL's linear layer then spreads it, so the
final attack is not an exact differential. It is a differential-linear parity
bias.

---

## 5. Output Difference

For each chosen input pair:

```text
x
x xor alpha_state
```

compute two reduced SWIRL rounds:

```text
y  = F(F(x))
y' = F(F(x xor alpha_state))
```

Then compute:

```text
Delta2 = y xor y'
```

`Delta2` is a 512-bit output difference split into eight 64-bit lanes:

```text
Delta2_lane0
Delta2_lane1
Delta2_lane2
Delta2_lane3
Delta2_lane4
Delta2_lane5
Delta2_lane6
Delta2_lane7
```

---

## 6. The Parity Test

The distinguisher checks only selected bits in output lanes 2 and 6.

Masks:

```text
M2 = 0x72ae08c4e7e81a1a
M6 = 0xe55c1188cfd03434
```

Test:

```text
parity(
  (Delta2_lane2 & M2) ^
  (Delta2_lane6 & M6)
) == 0
```

In C++:

```cpp
const uint64_t selected =
    (delta[2] & 0x72ae08c4e7e81a1aULL) ^
    (delta[6] & 0xe55c1188cfd03434ULL);

bool pass = parity64(selected) == 0;
```

For a random permutation, this parity should be a fair bit:

```text
0 half the time
1 half the time
```

For reduced two-round SWIRL, it leans toward zero.

---

## 7. Measured Result

Large verification run:

```text
samples = 16,000,000 chosen input pairs
passes  = 8,292,258
```

Probability:

```text
p = 8,292,258 / 16,000,000
  = 0.51826613
```

Bias from ideal:

```text
epsilon = p - 0.5
        = 0.01826613
```

Correlation:

```text
corr = |2p - 1|
     = 0.03653225
     ~= 2^-4.775
```

Statistical strength:

```text
z-score ~= 146
```

This is far too large to be sampling noise.

---

## 8. Plain-English Meaning

The full 512-bit output difference does not repeat.

The attack does not say:

```text
Delta2 is always equal to some fixed value.
```

Instead, it says:

```text
If we select a particular set of bits from Delta2 lanes 2 and 6,
then XOR those selected bits down to one bit,
that final bit is slightly more likely to be 0 than 1.
```

For a random function:

```text
0 about 50% of the time
1 about 50% of the time
```

For reduced two-round SWIRL:

```text
0 about 51.82% of the time
1 about 48.18% of the time
```

This is a statistical fingerprint.

---

## 9. Complexity Estimate

Let:

```text
epsilon = 0.01826613
corr    = 0.03653225
```

A simple pass-rate distinguisher needs roughly:

```text
1 / epsilon^2 ~= 2,997 pairs
```

For correlation-style estimates:

```text
1 / corr^2 ~= 749 pairs
```

Practical rough guide:

```text
2^10 pairs: weak visible signal
2^12 pairs: clear signal
2^15 pairs: strong signal
```

The exact confidence depends on the statistical decision rule, but the measured
bias is large enough to be detected with far fewer than the 16 million pairs used
for verification.

---

## 10. Why This Survives Two Rounds

Mix13 has local structure for shift-aligned differences. In isolation, some of
that structure is deterministic. SWIRL's diffusion layer prevents the exact
one-lane differential from continuing cleanly.

However, after one round, the dispersed difference still has enough structure
that the second Mix13 layer produces a biased parity relation.

The original internal relation was found before the final second-round linear
layer:

```text
third lane of the second Mix13 layer:
  bit13 xor bit44 is biased toward 0
```

Because the final MDS layer is linear over GF(2), that internal relation can be
converted into an equivalent final-output mask. That conversion produced:

```text
lane 2 mask = 0x72ae08c4e7e81a1a
lane 6 mask = 0xe55c1188cfd03434
```

So the final-output test is just the externally visible version of an internal
second-Mix13 parity bias.

---

## 11. What Happens at Three Rounds

The same final-output mask was carried forward and tested on three reduced
rounds.

Result:

```text
samples = 2,000,000 pairs
p       = 0.50038350
z       ~= 1.08
```

This is consistent with random noise.

Additional simple scans on the three-round output also found no meaningful
signal:

```text
maximum single-bit bias:                    |z| ~= 3.29
maximum bit[k] == bit[k+31] lane relation:  |z| ~= 3.85
Hamming weight:                             consistent with random
```

A targeted pre-final-linear scan of the third Mix13 layer also found no useful
continuation:

```text
best pre-linear relation:
  p    = 0.50111450
  corr ~= 2^-8.81
  z    ~= 3.15
```

After accounting for many tested masks, that is not significant.

So the current evidence is:

```text
1 round:
  exact Mix13 differential structure is visible

2 rounds:
  exact trail is dispersed, but a parity bias remains

3 rounds:
  tested parity bias disappears
```

---

## 12. Build and Run

Build with a C++20 compiler:

```bash
g++ -O3 -std=c++20 -march=native swirl_round2_distinguisher_attack.cpp -o r2dist
```

Run default verification:

```bash
./r2dist
```

The default sample size is:

```text
2^24 = 16,777,216 pairs
```

Run a smaller quick test:

```bash
./r2dist 20
```

Run with a three-round comparison:

```bash
./r2dist 22 --round3
```

The first argument is the exponent:

```text
./r2dist 24  means N = 2^24 pairs
```

---

## 13. Expected Output Shape

For two rounds, output should look like:

```text
Two reduced rounds:
  samples:             16777216
  mask parity == 0:    around 8695xxx
  probability:         around 0.518
  epsilon from 1/2:    around 1.8e-2
  correlation:         around 3.6e-2
  log2(correlation):   around -4.8
  z-score:             large
```

For three rounds with the same mask:

```text
Three reduced rounds:
  probability:         close to 0.5
  z-score:             noise-level
```

---

## 14. Security Interpretation

This is an impressive reduced-round distinguisher, but it is not a full attack
on the complete hash function.

It demonstrates:

```text
The Mix13 line leaves a measurable statistical footprint after two reduced
SWIRL rounds.
```

It does not demonstrate:

```text
a collision
a preimage attack
a second-preimage attack
a full-round distinguisher
a key-recovery attack
```

The important design lesson is that Mix13 has local differential-linear
structure, but SWIRL's diffusion appears to disperse it quickly. The tested
two-round signal does not visibly survive into the third reduced round.

---

## 15. Bottom Line

The reduced two-round attack is:

```text
Input difference:
  lane0 = 0x0200000008000000

Output test:
  parity((Delta2_lane2 & 0x72ae08c4e7e81a1a) ^
         (Delta2_lane6 & 0xe55c1188cfd03434)) == 0

Ideal probability:
  0.5

Measured reduced-SWIRL probability:
  0.51826613

Correlation:
  about 2^-4.78
```

This is best described as a two-round differential-linear distinguisher caused
by Mix13's local structure and partially surviving SWIRL's aggressive diffusion.

