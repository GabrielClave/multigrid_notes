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
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let lk = $lambda_k$
#let lmin = $lambda_min$
#let tb(x) = $tilde(bold(#x))$

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

The space $V = bb(R)^n$ is decomposed into $J$ subspaces $V_1, dots, V_J subset bb(R)^n$, where each $V_i$ has dimension $n_i <= n$. We have:
$ V = sum_(i = 1)^J V_i $
which means each $v in V$ can be decomposed into a (non-unique) $v = sum_i v_i quad v_i in V_i$

== Prolongation Matrices

$P_i in bb(R)^(n times n_i), quad P_i : bb(R)^(n_i) arrow bb(R)^n$ is a full column-rank matrix whose columns form a basis of $V_i$, representing the prolongation  operator from $bb(R)^(n_i)$ into $bb(R)^n$. \
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

We define $B in bb(R)^(n times n)$:
$ B = sum_(i=1)^J P_i A_i^(-1) P_i^T $

which is more or less an additive Schwartz preconditioner. \
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

= Notes

== Lemma 2.4

In the paper we have $ (T^(-1) v, v) = inf_(sum v_i = v) sum_(i=1)^J (T_i^(-1) v_i, v_i) $

In matrix algebra, the local iteration matrix $T_i = P_i A_i^(-1) P_i^T A$ 
is not symmetric under the standard euclidean inner product ($u^T v$). \
it is self-adjoint under the $A$-energy inner product $( u, v )_A = v^T A u$

In $bb(R)^n$, $T_i = P_i A_i^(-1) P_i^T A in bb(R)^(n times n)$ has rank $n_i < n$, making $T_i$ singular as an $n times n$ matrix.

Restricted to the subspace $V_i = op("range")(P_i)$,$quad T_i|_V_i$ is a bijection. 

for any $v_i = P_i x_i in V_i$ ($x_i in bb(R)^(n_i)$), the equation $T_i w = v_i$ for $w = P_i y_i in V_i$ yields:

$ P_i A_i^(-1) P_i^T A (P_i y_i) = P_i x_i ==> y_i = (P_i^T A P_i)^(-1) A_i x_i $

$ (T_i^(-1) v_i, v_i)_A = v_i^T A w = (P_i x_i)^T A [P_i (P_i^T A P_i)^(-1) A_i x_i] = x_i^T A_i x_i $