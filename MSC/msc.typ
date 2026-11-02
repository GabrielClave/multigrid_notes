// finite dimension translation of "THE METHOD OF SUBSPACE CORRECTIONS
// IN HILBERT SPACE JINCHAO XU AND LUDMIL ZIKATANOV"

#import "@preview/cetz:0.5.1"
// #import "figures.typ": *

#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  numbering: "1",
)

#set text(
  font: "xCharter",
  size: 11pt,
  lang: "en",
)

#set par(justify: true, leading: 0.65em)

// Custom Math Operators
#let range = math.op("range")
#let nul = math.op("null")
#let im = math.op("Im")
#let ker = math.op("Ker")
#let Ri = $overline(R)^(-1)$
#let Rb = $overline(R)$
#let Tb = $overline(T)$
#let Tbi = $overline(T_i)$
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let Ti = $T_i$
#let lk = $lambda_k$
#let lmin = $lambda_min$
#let tb(x) = $tilde(bold(#x))$
#let iff = $arrow.l.r.double.long$

// Document Title & Metadata
#align(center)[
  #text(size: 18pt, weight: "bold")[Subspace Correction Methods] \
  #v(8pt)
  #datetime.today().display("[month repr:long] [day], [year]")
]

#v(1.5em)

// Table of Contents
#outline(indent: auto)

= Subspace Decomposition

Let $V = bb(R)^n$. We consider the linear system:
$ A u = f $
where $A in bb(R)^(n times n)$ is symmetric positive definite

The space $V$ is decomposed into $J$ subspaces $V_1, dots, V_J subset bb(R)^n$, where each $V_i$ has dimension $n_i <= n$. We have:
$ V = sum_(i = 1)^J V_i $
which means each $v in V$ can be decomposed into a (non-unique) $v = sum_i v_i quad v_i in V_i$

== Prolongation Matrices

$P_i in bb(R)^(n times n_i), quad P_i : bb(R)^(n_i) arrow bb(R)^n$ is a full column-rank matrix whose columns form a basis of $V_i$, 
representing the prolongation  operator from $V_i = bb(R)^(n_i)$ into $V = bb(R)^n$. \
The corresponding restriction operator is $P_i^T in bb(R)^(n_i times n)$.

The global space decomposition condition (A0) takes the block-matrix form:
$ op("range")(P) = V, quad "where" P = mat(P_1 |, P_2 |, dots,| P_J) in bb(R)^(n times N) "with" N = sum_(i=1)^J n_i $

== Subspace Solvers

Exact Subspace Matrix: The restriction of $A$ to $V_i$ is $A_(i, "exact") = P_i^T A P_i in bb(R)^(n_i times n_i)$, which is SPD.

Approximate Local Solver: $A_i in bb(R)^(n_i times n_i)$ is an SPD matrix approximating $A_(i, "exact")$.

= Successive Subspace Corrections

The Method of Successive Subspace Corrections corrects residuals sequentially across subspaces $1$ to $J$.

== Algorithm:

- Compute local residual: $r_(i-1)^("local") = P_i^T A e_(i-1) in bb(R)^(n_i)$.
- Solve local subproblem: $x_i = A_i^(-1) r_(i-1)^("local") in bb(R)^(n_i)$.
- Update global iterate:
$ 
  e_i &= e_(i-1) - P_i x_i \
  &= (I - P_i A_i^(-1) P_i^T A ) e_(i-1)
$

We introduce the subspace iteration matrix $T_i$:
$ T_i = P_i A_i^(-1) P_i^T A in bb(R)^(n times n) $

The local update satisfies:
$ e_i = (I - T_i) e_(i-1) $

and so the overall error propagation matrix $E in bb(R)^(n times n)$ is:
$ E = (I - T_J)(I - T_(J-1)) dots (I - T_1) $

When $A_i = A_(i, "exact")$, $T_i$ is the $A$-orthogonal projection matrix $Pi_i$:
$ Pi_i = P_i (P_i^T A P_i)^(-1) P_i^T A $


= Parallel Subspace Corrections

The method of parallel subspace corrections solves local subproblems independently and sum their corrections.

== Algorithm

$
  e^(k+1) &= (I - sum_(i=1)^J P_i A_i^(-1) P_i^T A) e^k \
  &= (I - sum_(i=1)^J T_i ) e^k \
  &= ( I - T ) e^k
$

We define $B,T in bb(R)^(n times n)$:
$ B = sum_(i=1)^J P_i A_i^(-1) P_i^T \
T = B A $

$B$ is more or less an additive Schwartz preconditioner. \
We see that when all $A_i$ are SPD, $B$ is also SPD

It can be useful to see $B$ as
$ B = P tilde(A)^(-1) P^T $

where $tilde(A) = "blockdiag"(A_1, A_2, dots, A_J) in bb(R)^(N times N)$ is SPD.

$ T = B A = sum_(i=1)^J P_i A_i^(-1) P_i^T A = sum_(i=1)^J T_i $

== Invertibility characterization

Assume $A$ and all local solvers $A_i$ are SPD. The following conditions are equivalent:
1. The space decomposition condition holds: $op("range")(P) = bb(R)^n$.
2. $B$ is SPD on $bb(R)^n$ and $T = B A$ is non-singular.

Proof

$ forall v in V$, we have:

$
  x^T B x &= sum_(i=1)^J x^T P_i A_i^(-1) P_i^T x \
  &= sum_(i=1)^J (P_i^T x)^T A_i^(-1) P_i^T x
$

because $A_i^(-1) "is SPD," x^T B x >= 0$ and if $range(P) = V$:
$ x^T B x = 0 & arrow.double.l.r forall i quad P_i^T x = 0 \
  & arrow.double.l.r P^T x = 0 \
  & arrow.double.l.r x in ker(P^T) = range(P)^perp = {0} \
  & arrow.double.l.r x = 0 \

$

== Minimal energy decomposition

If condition (A0) holds ($op("range")(P) = V$) and each $A_i$ is SPD, then for any vector $v in V$:
$ v^T B^(-1) v = inf_(sum_(i=1)^J P_i x_i = v) sum_(i=1)^J x_i^T A_i x_i $

// interpretation: there is decomposition of $v$ of "minimal local energy $(v_i,v_i)_A_i^(-1)$" and its norm is $v^T B^(-1) v$

Proof:

The minimum is achieved for $x^* in bb(R)^N$:
$ x^* = tilde(A)^(-1) P^T B^(-1) v $

$ P x^* = P (tilde(A)^(-1) P^T B^(-1) v) = (P tilde(A)^(-1) P^T) B^(-1) v = B B^(-1) v = v $

$ (x^*)^T tilde(A) x^* &= (tilde(A)^(-1) P^T B^(-1) v)^T tilde(A) (tilde(A)^(-1) P^T B^(-1) v) \
&= v^T B^(-1) (P tilde(A)^(-1) P^T) B^(-1) v \
&= v^T B^(-1) v $

For any $x in bb(R)^N$ satisfying $P x = v$. 

We decompose $x = x^* + w$.
$ P x = v & arrow.double P x^* + P w = v \
& => P w = 0 \
& => w in op("ker")(P) $

Expanding $x^T tilde(A) x$:
$ x^T tilde(A) x &= (x^* + w)^T tilde(A) (x^* + w) \
&= (x^*)^T tilde(A) x^* + 2 w^T tilde(A) x^* + w^T tilde(A) w $

We have $ w^T tilde(A) x^* = 0$ :
$ w^T tilde(A) x^* = w^T tilde(A) (tilde(A)^(-1) P^T B^(-1) v) = w^T P^T B^(-1) v = (P w)^T B^(-1) v = 0 $

And because $tilde(A)$ is SPD, $w^T tilde(A) w >= 0$ for all $w$:

$ x^T tilde(A) x >= v^T B^(-1) v $

== Spectral Bound

The minimum eigenvalue of $T = B A$ satisfies:
$ 1 / (lmin (B A)) = sup_(v != 0) inf_(sum_(i=1)^J P_i x_i = v) frac(sum_(i=1)^J x_i^T A_i x_i, v^T A v) $

= Convergence

We want to prove that
$
||E||_A < 1
$

We assume
$
  ||I - T_i||_A < 1 \
  iff forall v in V, (T_i v, T_i v)_A < 2(T_i v, v)_A
$
Meaning there is $omega in [0,2],$
$
  forall v in V, (T_i v, T_i v)_A <= omega(T_i v, v)_A
$

== Error norm reduction

At step $i$, the error norm reduction is given by:
$
  ||(I - Ti)v||_i^2 &= ((I - T_i) v, (I - T_i) v)_A \
  &= ((I - T_i)^*(I - T_i) v, v)_A
$

where $M^*$ is the adjoint of $M$ with respect to $A$ \
We have:
$
  (I - T_i)^*(I - T_i) &= I - (Ti + Ti^* - Ti^*Ti) \
  & = I - Tbi
$
$
  Tbi = Ti + Ti^* - Ti^*Ti
$

$
  T_i^* &= A^(-1) Ti^T A \
  &= A^(-1) (P_i A_i^(-1) P_i^T A)^T A \
  &= P_i (A_i^(-1))^T P_i^T A
$

which we can think of the subspace correction applying the transpose of the local solver

and so:
$
    ||(I - Ti)v||_A^2 &= ((I - Tbi) v, v)_A \
  &= ||v||_A^2 - (Tbi v,v)_A
$

$(Tbi v,v)_A$ quantifies the norm reduction of the subspace correction

== range and kernel identities

$
  Ti = P_i A_i^(-1) P_i^T A \
  Ti^* = P_i A_i^(-T) P_i^T A \
  Tbi = Ti + Ti^* - Ti^*Ti \
$

$A_i$ is non-singular, $P_i$ is surjective and $range(P_i) = V_i$

We have
$
  range(Ti) = range(Ti^*) = V_i \
  range(Tbi) subset V_i
$

$
  x in ker(Ti) &iff P_i A_i^(-1) P_i^T A x = 0 \
  & iff x in ker(P_i^T A) = V_i^(perp_A)
$
$
  ker(Ti) = ker(Ti^*) = V_i^(perp_A) \
  ker(Tbi) subset V_i^(perp_A)
$

and from
$
  & ||(I - Ti)v||_A^2 = ||v||_A^2 - (Tbi v,v)_A iff (Tbi v,v)_A = ||v||_A^2 - ||(I - Ti)v||_i^2
$

using $||Ti v||^2_A <= omega(T_i v, v)_A$

$
  (Tbi v,v)_A &= 2(T_i v, v)_A - ||Ti v||^2_A \
  (Tbi v,v)_A &>= (2 - omega)/omega ||Ti v||^2_A  
$
which implies $ker(Ti) = V_i^(perp_A) subset ker(Tbi)$

As a result:
$
  ker(Tbi) = V_i^(perp_A) \
  range(Tbi) = V_i
$

== Bounding the operator norm

We quantify successive error norm reduction using $Tbi$ \
At step 1:
$
  ||(I - T_1)v||_A^2 = ||v||_A^2 - (Tb_1 v,v)_A  
$

At step 2:
$
  ||(I - T_2)(I - T_1)v||_A^2 = ||v||_A^2 - (Tb_1 v,v)_A - (Tb_2 (I - T_1)v, (I - T_1)v)_A
$

we introduce:
$
  E_i = (I - T_i)E_(i-1) \
  E_0 = I \
  E_J = E
$

$
  ||E v||_A^2 = ||v||_A^2 - sum_(i = 1)^(J) (Tb_i E_(i-1) v, E_(i-1)v)_A \
  // ||v||_A^2 - ||E v||_A^2 = sum_(i = 1)^(J) (Tb_i E_(i-1) v, E_(i-1)v)_A
$

= Notes

== Lemma 2.4

In the paper we have $ (T^(-1) v, v) = inf_(sum v_i = v) sum_(i=1)^J (T_i^(-1) v_i, v_i) $

In matrix algebra, the local iteration matrix $T_i = P_i A_i^(-1) P_i^T A$ 
is not symmetric under the standard euclidean inner product ($u^T v$). \
it is self-adjoint under the $A$-energy inner product $( u, v )_A = v^T A u$

In $bb(R)^n$, $T_i = P_i A_i^(-1) P_i^T A in bb(R)^(n times n)$ has rank $n_i < n$, making $T_i$ singular as an $n times n$ matrix.

Restricted to the subspace $V_i = op("range")(P_i)$,$quad T_i|_V_i$ is a bijection. 

for any $v_i = P_i x_i in V_i$ ($x_i in bb(R)^(n_i)$), the equation $T_i w = v_i$ for $w = P_i y_i in V_i$ yields:

$ P_i A_i^(-1) P_i^T A (P_i y_i) = P_i x_i ==> y_i & = (P_i^T A P_i)^(-1) A_i x_i \
&= A_(i, "exact")^(-1) A_i x_i $

$ (T_i^(-1) v_i, v_i)_A = v_i^T A w = (P_i x_i)^T A [P_i (P_i^T A P_i)^(-1) A_i x_i] = x_i^T A_i x_i $