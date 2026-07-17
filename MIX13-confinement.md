## Theorem: 28-dimensional MIX13 differential confinement

Let

[
F(x)=X_{31}!\left(C_2,X_{27}!\left(C_1X_{30}(x)\right)\right),
]

where all multiplications are modulo (2^{64}),

[
X_s(x)=x\oplus(x\gg s),
]

and

[
C_1=\texttt{0xBF58476D1CE4E5B9},
\qquad
C_2=\texttt{0x94D049BB133111EB}.
]

Choose

[
\alpha=\texttt{0x8000000200000008}
=2^{63}\oplus2^{33}\oplus2^3.
]

For every (x), define the MIX13 output difference

[
\beta(x)=F(x)\oplus F(x\oplus\alpha).
]

Then there is a unique (T(x)\in\mathbb F_2^{28}) such that

[
\boxed{\beta(x)=(T(x)\ll36)\oplus(T(x)\ll5)}.
]

Moreover, the span of all attainable (\beta(x)) has dimension exactly 28.

---

## 1. Difference through the first xor-shift

Because (X_{30}) is linear over XOR,

[
X_{30}(x)\oplus X_{30}(x\oplus\alpha)
=====================================

X_{30}(\alpha).
]

Now:

[
\alpha\gg30
===========

2^{33}\oplus2^3.
]

Therefore:

[
X_{30}(\alpha)
==============

# \alpha\oplus(\alpha\gg30)

2^{63}.
]

Writing

[
u=X_{30}(x),
\qquad
u'=X_{30}(x\oplus\alpha),
]

we obtain:

[
u'=u\oplus2^{63}.
]

Toggling the most-significant bit is equivalent to adding (2^{63}) modulo (2^{64}), so:

[
u'\equiv u+2^{63}\pmod{2^{64}}.
]

---

## 2. Difference through the first multiplication

Let:

[
v=C_1u,
\qquad
v'=C_1u'
\pmod{2^{64}}.
]

Since (C_1) is odd, write (C_1=2k+1). Then:

[
C_1 2^{63}
==========

# (2k+1)2^{63}

k2^{64}+2^{63}
\equiv2^{63}\pmod{2^{64}}.
]

Consequently:

[
v'
\equiv v+2^{63}
\equiv v\oplus2^{63}
\pmod{2^{64}}.
]

Thus the XOR difference remains exactly:

[
v\oplus v'=2^{63}.
]

---

## 3. Difference through the second xor-shift

Let:

[
w=X_{27}(v),
\qquad
w'=X_{27}(v').
]

Because (X_{27}) is XOR-linear,

[
w\oplus w'
==========

# X_{27}(v\oplus v')

X_{27}(2^{63}).
]

Therefore:

[
w\oplus w'
==========

2^{63}\oplus2^{36}.
]

This relation holds for every (x).

---

## 4. The second multiplication forces 36 zero low bits

Let:

[
a=C_2w,
\qquad
a'=C_2w'
\pmod{2^{64}}.
]

Since:

[
w'=w\oplus2^{63}\oplus2^{36},
]

toggling bit 63 always adds (2^{63}), while toggling bit 36 either adds or subtracts (2^{36}), depending on bit 36 of (w). Hence:

[
w'
\equiv
w+2^{63}+\varepsilon2^{36}
\pmod{2^{64}},
]

where:

[
\varepsilon=
\begin{cases}
+1,&w_{36}=0,\
-1,&w_{36}=1.
\end{cases}
]

Multiplication gives:

[
a'-a
\equiv
C_2\left(2^{63}+\varepsilon2^{36}\right)
\pmod{2^{64}}.
]

Because (C_2) is odd:

[
C_2 2^{63}\equiv2^{63}\pmod{2^{64}}.
]

Also,

[
C_2 2^{36}
\equiv
(C_2\bmod2^{28})2^{36}
\pmod{2^{64}}.
]

Thus:

[
a'-a
\equiv
2^{63}
+
\varepsilon(C_2\bmod2^{28})2^{36}
\pmod{2^{64}}.
]

Every term on the right is divisible by (2^{36}). Therefore:

[
a'\equiv a\pmod{2^{36}}.
]

So the lowest 36 bits of (a) and (a') are identical. Their XOR difference

[
\delta=a\oplus a'
]

must consequently satisfy:

[
\delta_i=0,\qquad 0\le i<36.
]

Hence there is a unique 28-bit word (T) such that:

[
\boxed{\delta=T\ll36}.
]

---

## 5. The final xor-shift duplicates (T)

The final MIX13 operation is (X_{31}). Therefore:

[
\begin{aligned}
\beta
&=X_{31}(a)\oplus X_{31}(a')\
&=X_{31}(a\oplus a')\
&=X_{31}(\delta)\
&=\delta\oplus(\delta\gg31).
\end{aligned}
]

Substituting (\delta=T\ll36):

[
\delta\gg31
===========

# (T\ll36)\gg31

T\ll5.
]

Therefore:

[
\boxed{
\beta=(T\ll36)\oplus(T\ll5)
}.
]

The two terms occupy disjoint bit ranges:

```text
T << 36: bits 36 ... 63
T <<  5: bits  5 ... 32
```

so no cancellation occurs between them.

---

# The 36 MIX13 constraints

Number bits from 0 at the least-significant bit.

From

[
\beta=(T\ll36)\oplus(T\ll5),
]

we obtain three groups of constraints.

### Eight fixed-zero bits

There are no (T)-bits in positions (0,\ldots,4):

[
\boxed{\beta_0=\beta_1=\beta_2=\beta_3=\beta_4=0.}
]

There are also no (T)-bits in positions (33,34,35):

[
\boxed{\beta_{33}=\beta_{34}=\beta_{35}=0.}
]

These give eight constraints.

### Twenty-eight equality constraints

For every (0\le i<28),

[
\beta_{5+i}=T_i
]

and:

[
\beta_{36+i}=T_i.
]

Therefore:

[
\boxed{
\beta_{5+i}\oplus\beta_{36+i}=0,
\qquad 0\le i<28.
}
]

These give 28 additional constraints.

The total is:

[
8+28=36.
]

They are linearly independent because:

* The eight zero constraints concern distinct bits.
* Each equality constraint contains its own unique pair of bits.
* The zero-bit positions do not overlap the equality positions.

Consequently, their common kernel has dimension:

[
64-36=28.
]

Define:

[
V=
\left{
(T\ll36)\oplus(T\ll5):
T\in\mathbb F_2^{28}
\right}.
]

The map:

[
A:\mathbb F_2^{28}\rightarrow\mathbb F_2^{64},
\qquad
A(T)=(T\ll36)\oplus(T\ll5)
]

is linear and injective, because (T) can be recovered from bits (36,\ldots,63). Hence:

[
\dim V=28.
]

So every MIX13 difference belongs to the explicitly defined 28-dimensional subspace (V).

---

# Why the dimension is exactly 28

The preceding argument proves containment in a 28-dimensional subspace. We must additionally prove that the attainable differences do not span some smaller subspace.

Consider the branch (w_{36}=0). In that branch:

[
a'\equiv a+D_+\pmod{2^{64}},
]

where:

[
D_+=d,2^{36},
]

and:

[
d=2^{27}+(C_2\bmod2^{28})
=\texttt{0x0B3111EB}.
]

In particular, (d) is odd.

Writing:

[
u=a\gg36\in\mathbb Z/2^{28}\mathbb Z,
]

the 28-bit value (T) is:

[
T=u\oplus(u+d\bmod2^{28}).
]

### Every (u) is attainable in this branch

Write:

[
w=r+2^{36}b,
]

where (0\le r<2^{36}), (b\in\mathbb Z/2^{28}\mathbb Z), and (b) is even because (w_{36}=0).

Let:

[
c=C_2\bmod2^{28}.
]

Then (c) is odd and therefore invertible modulo (2^{28}). The upper 28 bits of (C_2w) satisfy:

[
u\equiv q_r+cb\pmod{2^{28}},
]

where:

[
q_r=
\left\lfloor\frac{C_2r}{2^{36}}\right\rfloor
\bmod2^{28}.
]

For (r=0):

[
q_0=0.
]

Thus every even target (u) is obtained using:

[
b=c^{-1}u,
]

which is even.

For (r=1):

[
q_1=(C_2\gg36)\bmod2^{28}
=\texttt{0x094D049B},
]

which is odd. Thus every odd target (u) is obtained using:

[
b=c^{-1}(u-q_1),
]

which is again even.

Therefore every (u\in\mathbb Z/2^{28}\mathbb Z) occurs in the (w_{36}=0) branch.

### Lemma

For odd (d), the set

[
S_d=
\left{
u\oplus(u+d):
u\in\mathbb Z/2^{28}\mathbb Z
\right}
]

spans (\mathbb F_2^{28}).

#### Proof

Suppose a mask (m\in\mathbb F_2^{28}) is orthogonal to every element of (S_d). Then:

[
\langle m,u\oplus(u+d)\rangle=0
]

for every (u). Equivalently:

[
\langle m,u\rangle
==================

\langle m,u+d\rangle.
]

Define the linear Boolean function:

[
\ell_m(u)=\langle m,u\rangle.
]

The preceding equation says:

[
\ell_m(u)=\ell_m(u+d)
]

for every (u).

Because (d) is odd,

[
\gcd(d,2^{28})=1.
]

Therefore repeated addition by (d) visits every element of
(\mathbb Z/2^{28}\mathbb Z). Hence (\ell_m) must be constant over the entire space.

But:

[
\ell_m(0)=0,
]

so (\ell_m(u)=0) for every (u). This is possible only when:

[
m=0.
]

Thus (S_d^\perp={0}), and consequently:

[
\operatorname{span}(S_d)=\mathbb F_2^{28}.
]

This proves the lemma. ∎

Since the attainable (T)-values span (\mathbb F_2^{28}), and (A) is injective:

[
\boxed{
\dim\operatorname{span}{\beta(x):x\in\mathbb F_2^{64}}=28.
}
]

Therefore the 36 constraints listed above are complete: there are no additional universal homogeneous linear constraints on this MIX13 differential.

---

# Transport through one SWIRL round

Let the 512-bit state space be:

[
\mathcal S=(\mathbb F_2^{64})^8.
]

Let (e_{r,j}) denote the coordinate vector selecting bit (j) of lane (r), where:

[
0\le r<8,\qquad 0\le j<64.
]

Only lane 0 has input difference (\alpha). Immediately after the eight parallel MIX13 calls, the state difference is:

[
z=(\beta,0,0,0,0,0,0,0).
]

Define the embedding:

[
E:\mathbb F_2^{64}\rightarrow\mathcal S,
\qquad
E(\beta)=(\beta,0,\ldots,0).
]

The pre-linear difference space is:

[
U=E(V).
]

Since (E) is injective:

[
\dim U=\dim V=28.
]

---

## Invertibility of the SWIRL linear layer

Let (L) be the composition of:

1. The byte transpose.
2. The Whirlpool MDS transformation.

The byte transpose is a coordinate permutation and therefore invertible.

For each byte slice, the MDS matrix is:

[
A_{r,c}=m_{(r-c)\bmod8},
\qquad
m=(1,1,4,1,8,5,2,9).
]

Over (\operatorname{GF}(2^8)), using the Whirlpool polynomial
(x^8+x^4+x^3+x^2+1), its determinant is:

[
\det A=\texttt{0x1C}\ne0.
]

Thus the MDS map is invertible. Consequently:

[
L\in\operatorname{GL}(512,2).
]

The output-difference space is therefore:

[
W=L(U),
]

and:

[
\boxed{\dim W=28.}
]

Hence:

[
\boxed{\dim W^\perp=512-28=484.}
]

These are the 484 independent probability-1 output parity relations.

---

# Complete basis of the 484 pre-linear constraints

The orthogonal space (U^\perp) has the following basis.

## 1. Eight fixed-zero MIX13 masks

For:

[
j\in Z={0,1,2,3,4,33,34,35},
]

define:

[
q_j^{(0)}=e_{0,j}.
]

These enforce:

[
\langle q_j^{(0)},z\rangle=0.
]

Count: (8).

## 2. Twenty-eight repeated-bit masks

For (0\le i<28), define:

[
q_i^{(=)}
=========

e_{0,5+i}\oplus e_{0,36+i}.
]

These enforce:

[
z_{0,5+i}\oplus z_{0,36+i}=0.
]

Count: (28).

## 3. Four hundred forty-eight inactive-lane masks

For:

[
1\le r<8,\qquad 0\le j<64,
]

define:

[
q_{r,j}^{(\mathrm{inactive})}=e_{r,j}.
]

These enforce that lanes (1,\ldots,7) are zero before the linear layer.

Count:

[
7\cdot64=448.
]

The total number is:

[
8+28+448=484.
]

These masks are independent and form a complete basis of (U^\perp).

---

# Transporting all 484 masks through (L^{-T})

Let the observable one-round output difference be:

[
y=Lz.
]

For any pre-linear constraint mask (q\in U^\perp), define:

[
\boxed{m=L^{-T}q}.
]

Then:

[
\begin{aligned}
\langle m,y\rangle
&=\langle L^{-T}q,Lz\rangle\
&=\langle q,z\rangle\
&=0.
\end{aligned}
]

Therefore the complete externally observable mask basis is:

### Eight transported zero-bit masks

[
\boxed{
m_j^{(0)}=L^{-T}e_{0,j},
\qquad
j\in{0,1,2,3,4,33,34,35}.
}
]

### Twenty-eight transported equality masks

[
\boxed{
m_i^{(=)}
=========

L^{-T}
\left(e_{0,5+i}\oplus e_{0,36+i}\right),
\qquad 0\le i<28.
}
]

### Four hundred forty-eight transported inactive-lane masks

[
\boxed{
m_{r,j}^{(\mathrm{inactive})}
=============================

L^{-T}e_{r,j},
\qquad
1\le r<8,\quad0\le j<64.
}
]

For every chosen pair and every one of these masks:

[
\boxed{\langle m,y\rangle=0}
]

with probability exactly one.

Because (L^{-T}) is invertible, it preserves linear independence. Thus these 484 transported masks are independent.

They are also complete: if an output mask (m) vanishes on every attainable output difference, then:

[
L^Tm\in U^\perp.
]

Since the (q)-masks above form a basis of (U^\perp), (L^Tm) is a linear combination of them, and therefore (m) is the corresponding linear combination of their (L^{-T}) images.

---

## Final theorem

For the chosen lane difference

[
\alpha=\texttt{0x8000000200000008},
]

the MIX13 output difference satisfies:

[
\boxed{\beta=(T\ll36)\oplus(T\ll5)}
]

for a 28-bit (T). Its attainable values span an exactly 28-dimensional subspace of (\mathbb F_2^{64}), characterized by exactly 36 independent equations:

[
\beta_{0\ldots4}=0,
\qquad
\beta_{33\ldots35}=0,
\qquad
\beta_{5+i}=\beta_{36+i}.
]

When lifted into lane 0 of SWIRL and transported through the invertible byte-transpose/MDS layer (L), the one-round output differences span an exactly 28-dimensional subspace of (\mathbb F_2^{512}). Its orthogonal complement has dimension:

[
\boxed{512-28=484},
]

with the complete constraint basis obtained by applying:

[
\boxed{L^{-T}}
]

to the 36 MIX13 constraints and the 448 inactive-lane constraints.
