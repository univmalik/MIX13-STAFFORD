# Proof of the 28-Dimensional MIX13 Differential Theorem

## Abstract

This document proves that, for the chosen input XOR difference

$$
\alpha=\mathtt{0x8000000200000008},
$$

every output difference of Stafford's MIX13 permutation is of the form

$$
\beta=(T\ll36)\oplus(T\ll5)
$$

for a 28-bit word $T$. The attainable differences span an exactly 28-dimensional subspace of $\mathbb F_2^{64}$. This gives exactly 36 independent homogeneous linear constraints at the MIX13 output. When the difference is embedded in one lane of SWIRL and transported through its invertible byte-transpose/MDS layer $L$, the resulting one-round differences span an exactly 28-dimensional subspace of $\mathbb F_2^{512}$. Its orthogonal complement therefore has dimension $512-28=484$. A complete basis of these 484 externally observable constraints is obtained by applying $L^{-T}$ to the 36 MIX13 constraints and the 448 inactive-lane constraints.

---

## 1. Definitions and notation

All additions and multiplications of 64-bit words are performed in

$$
\mathbb Z/2^{64}\mathbb Z.
$$

XOR differences and linear relations are considered over $\mathbb F_2$.

For $s\geq 0$, define the XOR-shift map

$$
X_s(x)=x\oplus(x\gg s).
$$

The MIX13 permutation is

$$
F(x)=X_{31}\!\left(C_2X_{27}\!\left(C_1X_{30}(x)\right)\right),
$$

where

$$
C_1=\mathtt{0xBF58476D1CE4E5B9}
$$

and

$$
C_2=\mathtt{0x94D049BB133111EB}.
$$

The chosen input difference is

$$
\alpha=\mathtt{0x8000000200000008}
      =2^{63}\oplus2^{33}\oplus2^3.
$$

For every $x$, define

$$
\beta(x)=F(x)\oplus F(x\oplus\alpha).
$$

Bits are numbered from $0$ at the least-significant bit to $63$ at the most-significant bit.

---

## 2. Main theorem

### Theorem 1

For every $x\in\mathbb F_2^{64}$, there exists a unique $T(x)\in\mathbb F_2^{28}$ such that

$$
\boxed{\beta(x)=(T(x)\ll36)\oplus(T(x)\ll5).}
$$

Furthermore,

$$
\dim\operatorname{span}\{\beta(x):x\in\mathbb F_2^{64}\}=28.
$$

Consequently, the MIX13 output differences satisfy exactly 36 independent homogeneous linear constraints.

---

## 3. Propagation through the first XOR-shift

The map $X_{30}$ is linear over XOR. Therefore

$$
X_{30}(x)\oplus X_{30}(x\oplus\alpha)=X_{30}(\alpha).
$$

We have

$$
\alpha\gg30=2^{33}\oplus2^3.
$$

Hence

$$
X_{30}(\alpha)
=\alpha\oplus(\alpha\gg30)
=2^{63}.
$$

Let

$$
u=X_{30}(x)
$$

and

$$
u'=X_{30}(x\oplus\alpha).
$$

Then

$$
u'=u\oplus2^{63}.
$$

Toggling the most-significant bit is equivalent to adding $2^{63}$ modulo $2^{64}$, so

$$
u'\equiv u+2^{63}\pmod{2^{64}}.
$$

---

## 4. Propagation through the first multiplication

Let

$$
v=C_1u
$$

and

$$
v'=C_1u'
$$

modulo $2^{64}$. Since $C_1$ is odd, write

$$
C_1=2k+1.
$$

Then

$$
C_1 2^{63}
=(2k+1)2^{63}
=k2^{64}+2^{63}
\equiv2^{63}\pmod{2^{64}}.
$$

Consequently,

$$
v'\equiv v+2^{63}\equiv v\oplus2^{63}\pmod{2^{64}}.
$$

Thus the XOR difference after the first multiplication remains exactly

$$
v\oplus v'=2^{63}.
$$

---

## 5. Propagation through the second XOR-shift

Let

$$
w=X_{27}(v)
$$

and

$$
w'=X_{27}(v').
$$

Since $X_{27}$ is XOR-linear,

$$
w\oplus w'
=X_{27}(v\oplus v')
=X_{27}(2^{63}).
$$

Therefore

$$
\boxed{w\oplus w'=2^{63}\oplus2^{36}.}
$$

This equality holds for every $x$.

---

## 6. The second multiplication forces 36 zero low bits

Let

$$
a=C_2w
$$

and

$$
a'=C_2w'
$$

modulo $2^{64}$. Since

$$
w'=w\oplus2^{63}\oplus2^{36},
$$

toggling bit $63$ always adds $2^{63}$, while toggling bit $36$ either adds or subtracts $2^{36}$ according to bit $36$ of $w$. Thus

$$
w'\equiv w+2^{63}+\varepsilon2^{36}\pmod{2^{64}},
$$

where

$$
\varepsilon=
\begin{cases}
+1,&w_{36}=0,\\
-1,&w_{36}=1.
\end{cases}
$$

It follows that

$$
a'-a
\equiv C_2\left(2^{63}+\varepsilon2^{36}\right)
\pmod{2^{64}}.
$$

Because $C_2$ is odd,

$$
C_2 2^{63}\equiv2^{63}\pmod{2^{64}}.
$$

Also,

$$
C_2 2^{36}
\equiv(C_2\bmod2^{28})2^{36}
\pmod{2^{64}}.
$$

Therefore

$$
a'-a
\equiv
2^{63}+\varepsilon(C_2\bmod2^{28})2^{36}
\pmod{2^{64}}.
$$

Every term on the right-hand side is divisible by $2^{36}$. Hence

$$
a'\equiv a\pmod{2^{36}}.
$$

The lowest 36 bits of $a$ and $a'$ are consequently identical. Their XOR difference

$$
\delta=a\oplus a'
$$

must satisfy

$$
\delta_i=0
\qquad\text{for }0\leq i<36.
$$

There is therefore a unique 28-bit word $T$ such that

$$
\boxed{\delta=T\ll36.}
$$

---

## 7. The final XOR-shift duplicates the 28-bit word

The final MIX13 operation is $X_{31}$. Therefore

$$
\begin{aligned}
\beta
&=X_{31}(a)\oplus X_{31}(a')\\
&=X_{31}(a\oplus a')\\
&=X_{31}(\delta)\\
&=\delta\oplus(\delta\gg31).
\end{aligned}
$$

Substituting $\delta=T\ll36$ gives

$$
\delta\gg31=(T\ll36)\gg31=T\ll5.
$$

Hence

$$
\boxed{\beta=(T\ll36)\oplus(T\ll5).}
$$

The two terms occupy disjoint bit ranges:

$$
T\ll36:\quad 36,37,\ldots,63,
$$

and

$$
T\ll5:\quad 5,6,\ldots,32.
$$

Thus there is no cancellation between the two copies of $T$.

---

## 8. The 36 MIX13 constraints

The representation

$$
\beta=(T\ll36)\oplus(T\ll5)
$$

immediately gives 36 homogeneous linear constraints.

### 8.1 Eight fixed-zero constraints

There are no active $T$-bits in positions $0$ through $4$, so

$$
\boxed{\beta_0=\beta_1=\beta_2=\beta_3=\beta_4=0.}
$$

There are also no active $T$-bits in positions $33$ through $35$, so

$$
\boxed{\beta_{33}=\beta_{34}=\beta_{35}=0.}
$$

### 8.2 Twenty-eight repeated-bit constraints

For every $0\leq i<28$,

$$
\beta_{5+i}=T_i
$$

and

$$
\beta_{36+i}=T_i.
$$

Therefore

$$
\boxed{\beta_{5+i}\oplus\beta_{36+i}=0,
\qquad0\leq i<28.}
$$

There are $8+28=36$ constraints in total. They are independent because the eight fixed-zero constraints use distinct coordinates, while each repeated-bit constraint contains its own distinct pair of coordinates and none of those coordinates overlaps a fixed-zero coordinate.

Define

$$
V=
\left\{
(T\ll36)\oplus(T\ll5):T\in\mathbb F_2^{28}
\right\}.
$$

The map

$$
A:\mathbb F_2^{28}\longrightarrow\mathbb F_2^{64}
$$

defined by

$$
A(T)=(T\ll36)\oplus(T\ll5)
$$

is linear and injective, because $T$ can be recovered from bits $36$ through $63$ of $A(T)$. Hence

$$
\boxed{\dim V=28.}
$$

---

## 9. Proof that the dimension is exactly 28

The preceding sections prove that every output difference lies in a 28-dimensional space. We now prove that the attainable differences span that entire space rather than a smaller subspace.

Consider the branch $w_{36}=0$. In this branch,

$$
a'\equiv a+D_+\pmod{2^{64}},
$$

where

$$
D_+=d2^{36}
$$

and

$$
d=2^{27}+(C_2\bmod2^{28})
 =\mathtt{0x0B3111EB}.
$$

In particular, $d$ is odd. Writing

$$
u=a\gg36\in\mathbb Z/2^{28}\mathbb Z,
$$

the corresponding 28-bit difference is

$$
T=u\oplus(u+d\bmod2^{28}).
$$

### Lemma 1: Every $u$ is attainable in the branch $w_{36}=0$

Write

$$
w=r+2^{36}b,
$$

where $0\leq r<2^{36}$ and $b\in\mathbb Z/2^{28}\mathbb Z$. The condition $w_{36}=0$ means that $b$ is even.

Let

$$
c=C_2\bmod2^{28}.
$$

Since $c$ is odd, it is invertible modulo $2^{28}$. The upper 28 bits of $C_2w$ satisfy

$$
u\equiv q_r+cb\pmod{2^{28}},
$$

where

$$
q_r=
\left\lfloor\frac{C_2r}{2^{36}}\right\rfloor
\bmod2^{28}.
$$

For $r=0$,

$$
q_0=0.
$$

Thus every even target $u$ is obtained by choosing

$$
b=c^{-1}u,
$$

which is even.

For $r=1$,

$$
q_1=(C_2\gg36)\bmod2^{28}
=\mathtt{0x094D049B},
$$

which is odd. Every odd target $u$ is therefore obtained by choosing

$$
b=c^{-1}(u-q_1),
$$

which is also even. Hence every $u\in\mathbb Z/2^{28}\mathbb Z$ is attained in this branch. $\square$

### Lemma 2: The attainable $T$-values span $\mathbb F_2^{28}$

For odd $d$, define

$$
S_d=
\left\{
u\oplus(u+d):u\in\mathbb Z/2^{28}\mathbb Z
\right\}.
$$

Suppose that a mask $m\in\mathbb F_2^{28}$ is orthogonal to every element of $S_d$. Then

$$
\langle m,u\oplus(u+d)\rangle=0
$$

for every $u$, and therefore

$$
\langle m,u\rangle=\langle m,u+d\rangle.
$$

Define

$$
\ell_m(u)=\langle m,u\rangle.
$$

The preceding equality states that

$$
\ell_m(u)=\ell_m(u+d)
$$

for every $u$. Since $d$ is odd,

$$
\gcd(d,2^{28})=1.
$$

Repeated addition by $d$ therefore visits every element of $\mathbb Z/2^{28}\mathbb Z$. The function $\ell_m$ must consequently be constant on the entire space. But

$$
\ell_m(0)=0,
$$

so

$$
\ell_m(u)=0
$$

for every $u$. A linear functional that vanishes on every binary vector must have

$$
m=0.
$$

Thus

$$
S_d^\perp=\{0\},
$$

which implies

$$
\boxed{\operatorname{span}(S_d)=\mathbb F_2^{28}.}
$$

This proves Lemma 2. $\square$

Because the attainable $T$-values span $\mathbb F_2^{28}$ and $A$ is injective,

$$
\boxed{
\dim\operatorname{span}\{\beta(x):x\in\mathbb F_2^{64}\}=28.
}
$$

The 36 constraints in Section 8 are therefore complete: there are no additional universal homogeneous linear constraints on this MIX13 differential.

---

## 10. Transport through one SWIRL round

Let the 512-bit state space be

$$
\mathcal S=(\mathbb F_2^{64})^8.
$$

Let $e_{r,j}$ denote the coordinate vector selecting bit $j$ of lane $r$, where

$$
0\leq r<8
$$

and

$$
0\leq j<64.
$$

Only lane $0$ has the input difference $\alpha$. Immediately after the eight parallel MIX13 operations, the state difference is

$$
z=(\beta,0,0,0,0,0,0,0).
$$

Define the embedding

$$
E:\mathbb F_2^{64}\longrightarrow\mathcal S
$$

by

$$
E(\beta)=(\beta,0,\ldots,0).
$$

The pre-linear difference space is

$$
U=E(V).
$$

Since $E$ is injective,

$$
\boxed{\dim U=28.}
$$

---

## 11. Invertibility of the SWIRL linear layer

Let $L$ be the composition of the byte transpose and the Whirlpool MDS transformation.

The byte transpose is a coordinate permutation and is therefore invertible.

For each byte slice, the MDS matrix is defined by

$$
M_{r,c}=m_{(r-c)\bmod8},
$$

where

$$
m=(1,1,4,1,8,5,2,9).
$$

Over $\operatorname{GF}(2^8)$ with reduction polynomial

$$
x^8+x^4+x^3+x^2+1,
$$

its determinant is

$$
\det M=\mathtt{0x1C}\neq0.
$$

Therefore the MDS transformation is invertible. Hence

$$
\boxed{L\in\operatorname{GL}(512,2).}
$$

The one-round output-difference space is

$$
W=L(U).
$$

Since $L$ is invertible,

$$
\boxed{\dim W=\dim U=28.}
$$

It follows that

$$
\boxed{\dim W^\perp=512-28=484.}
$$

---

## 12. Complete basis of the 484 pre-linear constraints

We use the standard binary inner product

$$
\langle p,z\rangle
=\bigoplus_{k=0}^{511}p_kz_k.
$$

The following masks form a basis of $U^\perp$.

### 12.1 Eight fixed-zero MIX13 masks

Let

$$
Z=\{0,1,2,3,4,33,34,35\}.
$$

For every $j\in Z$, define

$$
q_j^{(0)}=e_{0,j}.
$$

These masks enforce

$$
\langle q_j^{(0)},z\rangle=0.
$$

Count:

$$
8.
$$

### 12.2 Twenty-eight repeated-bit masks

For every $0\leq i<28$, define

$$
q_i^{(=)}=e_{0,5+i}\oplus e_{0,36+i}.
$$

These masks enforce

$$
z_{0,5+i}\oplus z_{0,36+i}=0.
$$

Count:

$$
28.
$$

### 12.3 Four hundred forty-eight inactive-lane masks

For

$$
1\leq r<8
$$

and

$$
0\leq j<64,
$$

define

$$
q_{r,j}^{(\mathrm{inactive})}=e_{r,j}.
$$

These enforce that lanes $1$ through $7$ are zero before the linear layer.

Count:

$$
7\cdot64=448.
$$

The total number of masks is

$$
8+28+448=484.
$$

They are independent and form a complete basis of $U^\perp$.

---

## 13. Transporting all 484 constraints through $L^{-T}$

Let the observable one-round output difference be

$$
y=Lz.
$$

For every pre-linear constraint mask $q\in U^\perp$, define

$$
\boxed{p=L^{-T}q.}
$$

Then

$$
\begin{aligned}
\langle p,y\rangle
&=\langle L^{-T}q,Lz\rangle\\
&=\langle q,z\rangle\\
&=0.
\end{aligned}
$$

Thus every pre-linear constraint becomes an externally observable output parity constraint.

### 13.1 Eight transported zero-bit masks

For every

$$
j\in\{0,1,2,3,4,33,34,35\},
$$

define

$$
\boxed{p_j^{(0)}=L^{-T}e_{0,j}.}
$$

### 13.2 Twenty-eight transported equality masks

For every $0\leq i<28$, define

$$
\boxed{
p_i^{(=)}
=L^{-T}\left(e_{0,5+i}\oplus e_{0,36+i}\right).
}
$$

### 13.3 Four hundred forty-eight transported inactive-lane masks

For

$$
1\leq r<8
$$

and

$$
0\leq j<64,
$$

define

$$
\boxed{p_{r,j}^{(\mathrm{inactive})}=L^{-T}e_{r,j}.}
$$

For every chosen pair and every mask in these three families,

$$
\boxed{\langle p,y\rangle=0}
$$

with probability exactly one.

Since $L^{-T}$ is invertible, it preserves linear independence. Therefore the 484 transported masks are independent.

They are also complete. If an output mask $p$ vanishes on every attainable output difference, then

$$
L^Tp\in U^\perp.
$$

Since the 484 masks in Section 12 form a basis of $U^\perp$, the vector $L^Tp$ is a linear combination of them. Consequently, $p$ is the corresponding linear combination of their $L^{-T}$ images.

Hence

$$
\boxed{W^\perp=L^{-T}(U^\perp).}
$$

---

## 14. Equivalent joint membership test

Instead of evaluating 484 transported parity masks separately, an implementation may compute

$$
z=L^{-1}y
$$

and test the original pre-linear conditions:

$$
z_r=0
\qquad\text{for }1\leq r<8,
$$

$$
z_{0,j}=0
\qquad\text{for }j\in\{0,1,2,3,4,33,34,35\},
$$

and

$$
z_{0,5+i}=z_{0,36+i}
\qquad\text{for }0\leq i<28.
$$

This joint membership test is mathematically equivalent to testing all 484 output relations.

---

## 15. Final conclusion

For

$$
\alpha=\mathtt{0x8000000200000008},
$$

the MIX13 output difference always satisfies

$$
\boxed{\beta=(T\ll36)\oplus(T\ll5)}
$$

for a 28-bit word $T$. The attainable differences span an exactly 28-dimensional subspace of $\mathbb F_2^{64}$. This subspace is characterized by exactly 36 independent equations:

$$
\beta_j=0
\qquad
j\in\{0,1,2,3,4,33,34,35\},
$$

and

$$
\beta_{5+i}=\beta_{36+i}
\qquad
0\leq i<28.
$$

When lifted into lane $0$ of SWIRL and transported through the invertible byte-transpose/MDS layer $L$, the one-round output differences span an exactly 28-dimensional subspace of $\mathbb F_2^{512}$. Its orthogonal complement has dimension

$$
\boxed{512-28=484.}
$$

A complete basis of the 484 externally observable probability-one constraints is obtained by applying

$$
\boxed{L^{-T}}
$$

to the 36 MIX13 constraint masks and the 448 inactive-lane masks.
