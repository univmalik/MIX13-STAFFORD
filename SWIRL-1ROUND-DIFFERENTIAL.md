# One-Round SWIRL512 XOR Differential Candidates

Date: 23 June 2026

This report gives high-probability XOR differential candidates for one SWIRL512 data-path round only:

```text
mix64_each_lane -> byte_transpose -> Whirlpool linear diffusion
```

The model intentionally excludes:

- data-path round constants, including `S[0] ^= WMIX_RC[r]`;
- XOR RoundKey;
- key schedule;
- compression feed-forward.

The reduced one-round map is:

```text
R(S) = LinearDiffusion(ByteTranspose(mix64(S[0]), ..., mix64(S[7])))
```

The source `mix64` is the Stafford Mix13 / SplitMix64 finalizer:

```text
z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
z = (z ^ (z >> 27)) * 0x94d049bb133111eb
z =  z ^ (z >> 31)
```

All words below are 64-bit big-endian state lanes in the same lane order as `swirl512_keyed_fast.cpp`.

## Lifting Rule

The SWIRL512 data round first applies `mix64` independently to each of the eight lanes. The byte transpose and Whirlpool diffusion are linear over XOR.

Therefore, if a source-lane Mix13 differential is:

```text
alpha -> beta
```

and the 512-bit input difference places `alpha` only in source lane `i`, then the one-round 512-bit output difference is:

```text
out[j] = RAW[i][byte_j(beta)]  for j = 0..7
```

where `RAW[i][b]` is the Whirlpool linear-diffusion table generated from:

```text
MDS = {1, 1, 4, 1, 8, 5, 2, 9}
RAW[i][b] = pack8(gf_mul(b, MDS[(row - i) & 7]) for row = 0..7)
```

The 512-bit differential probability is the same as the underlying 64-bit Mix13 differential probability. The remaining SWIRL512 steps in this reduced round are deterministic linear transforms, so they relabel `beta` but do not change the event probability.

A fixed data-path round constant before `mix64` also does not change XOR differential probability: `x` is simply sampled as `x ^ constant`. It is omitted here exactly as requested.

## Primary Candidate

Strongest confirmed Mix13 seed from the supplied statistical report:

```text
alpha = 0x8000000200000008
beta  = 0x53112150a62242a0
p     = 9.773306082934e-4 ~= 2^-9.998866
```

Lifted to source lane 0:

```text
input difference:
8000000200000008 0000000000000000 0000000000000000 0000000000000000
0000000000000000 0000000000000000 0000000000000000 0000000000000000

output difference:
53535153a202a6f1 1111441188552299 2121842115a54234 50505d50ba0da0ea
a6a6a2a6590451ff 222288220daa442f 424215422a578468 a0a0baa0691a5dc9
```

This gives a one-round SWIRL512 exact 512-bit differential near probability `2^-10`, or roughly one hit per 1024 random paired states.

## Lane-0 Candidate Table

All eight cluster candidates below use the same lane-0 input difference:

```text
8000000200000008 0000000000000000 0000000000000000 0000000000000000
0000000000000000 0000000000000000 0000000000000000 0000000000000000
```

| Rank | Mix13 beta | Reported log2(p) | Lifted one-round SWIRL512 output difference |
|---:|---|---:|---|
| 1 | `53112150a62242a0` | -9.998866 | `53535153a202a6f1 1111441188552299 2121842115a54234 50505d50ba0da0ea a6a6a2a6590451ff 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 2 | `d3112151a62242a0` | -9.999185 | `d3d36bd3d6b8bb05 1111441188552299 2121842115a54234 51515951b208a2e3 a6a6a2a6590451ff 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 3 | `b3112151662242a0` | -9.999263 | `b3b3f6b3f1457b42 1111441188552299 2121842115a54234 51515951b208a2e3 6666856617e3cc71 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 4 | `cd1121519a2242a0` | -9.999727 | `cdcd13cd26de87eb 1111441188552299 2121842115a54234 51515951b208a2e3 9a9a529aa4c8293e 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 5 | `4d1121509a2242a0` | -10.000093 | `4d4d294d52649a1f 1111441188552299 2121842115a54234 50505d50ba0da0ea 9a9a529aa4c8293e 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 6 | `d5112151aa2242a0` | -10.000392 | `d5d573d5e6a6b733 1111441188552299 2121842115a54234 51515951b208a2e3 aaaa92aa39384993 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 7 | `b51121516a2242a0` | -10.000709 | `b5b5eeb5c15b7774 1111441188552299 2121842115a54234 51515951b208a2e3 6a6ab56a77dfd41d 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |
| 8 | `55112150aa2242a0` | -10.000764 | `55554955921caac7 1111441188552299 2121842115a54234 50505d50ba0da0ea aaaa92aa39384993 222288220daa442f 424215422a578468 a0a0baa0691a5dc9` |

The eight-member cluster has combined reported Mix13 probability about `2^-7`. Lifted to this reduced SWIRL512 round, this means about one pair in 128 should land in one of these eight exact 512-bit output differences.

## Strong Two-Bit Candidate

Strongest tested two-bit Mix13 seed from the supplied report:

```text
alpha = 0x0200000008000000
beta  = 0x88408901d0811201
p     = 4.203915596008e-4 ~= 2^-11.215979
```

Lifted to source lane 0:

```text
input difference:
0200000008000000 0000000000000000 0000000000000000 0000000000000000
0000000000000000 0000000000000000 0000000000000000 0000000000000000

output difference:
88881a8834920dbc 40401d403a5d807a 89891e893c970fb5 0101040108050209
d0d067d0ceb7bd1e 81813e817cbf1ffd 12124812905a2482 0101040108050209
```

## Other Source Lanes

The same Mix13 seed can be placed in any source lane `i = 0..7`. The input difference moves to that lane, and the output words are recomputed with `RAW[i]`.

For example, placing the strongest three-bit seed in source lane 3 gives:

```text
input difference:
0000000000000000 0000000000000000 0000000000000000 8000000200000008
0000000000000000 0000000000000000 0000000000000000 0000000000000000

output difference:
02a6f153535153a2 5522991111441188 a542342121842115 0da0ea50505d50ba
0451ffa6a6a2a659 aa442f222288220d 578468424215422a 1a5dc9a0a0baa069
```

Its probability remains the same as the source Mix13 differential, about `2^-10`.

## Statistical Sanity Check

A reduced-round validator was run over random 512-bit states with:

```text
samples = 2^22 = 4,194,304 per checked candidate
```

The validator computed:

```text
one_round(S) XOR one_round(S XOR input_difference)
```

and compared against the lifted 512-bit output difference.

| Candidate | Hits in 2^22 | Measured p | Expected source log2(p) |
|---|---:|---:|---:|
| strongest three-bit, lane 0 | 4195 | 0.00100016593933 | -9.998866 |
| strongest three-bit, lane 3 | 4045 | 0.000964403152466 | -9.998866 |
| strongest two-bit, lane 0 | 1817 | 0.000433206558228 | -11.215979 |

These checks are consistent with the supplied Mix13 probabilities at this sample size. They also confirm that the nonzero source lane is handled correctly by the byte transpose and Whirlpool RAW/MDS lifting.

## Caveats

These are high-probability one-round data-path candidates, not a proof of the global maximum differential probability for SWIRL512. The candidates are inherited from confirmed Mix13 statistical differentials and then lifted exactly through SWIRL512's deterministic linear layer.

They do not include the full compression function, key schedule behavior, feed-forward, or multi-round trail probability.
