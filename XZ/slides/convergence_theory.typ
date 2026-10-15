// Import the Touying package and the Metropolis theme
#import "@preview/touying:0.8.0": *
#import themes.metropolis: *
#import "@preview/cetz:0.5.1"
#import "convergence_theory_figures.typ": *
#import "../figures.typ": *

// #set text(font: "Fira Sans")

// Initialize the theme
#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Multigrid Convergence],
    author: [Gabriel Clave],
    date: datetime.today().display(),
  )
)

// #show heading.where(level: 1): set text(size: 20pt)

#let range = math.op("range")
#let nul = math.op("null")
#let im = math.op("Im")
#let ker = math.op("ker")
#let Ri = $overline(R)^(-1)$
#let Rb = $overline(R)$
#let Tb = $overline(T)$
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let lk = $lambda_k$
#let lmin = $lambda_min$

// figures

// Generate the title slide
#title-slide()

= Multigrid recap

== Core idea (picture)

#align(center)[
  #only("1")[#fig-amg(step: 1)]
  #only("2")[#fig-amg(step: 2)]
  #only("3")[#fig-amg(step: 3)]
]

== Core idea (text)

A multigrid solver relies on three core ingredients:

- A *Smoother $R$* eliminates error in a subspace $H$.
- A *Coarse Grid* & Interpolation ($P$)
- A *Coarse Solver* (can recursively be another multigrid cycle).

*Fundamental Principle:*\
The central objective when designing a multigrid method is to construct $P$ and $S$ such that $range(P)$ and $H$ are *complementary*:
$ V approx range(P) xor_A H $

Ensuring that error components missed by the smoother are effectively captured and eliminated by the coarse grid correction.

== Some notations (smoothers)
The first component of a multigrid solver is a smoother.\
It is a linear operator $R : V arrow.r.long V$ that can be used to solve the system through the  iteration:
$
  v^((k+1)) = v^((k)) + R(f - A v^((k)))
$
with associated error propagation operator
$
  S = I - R A
$

the error iterate verifying:
$
  e^((k+1)) = S e^(k)
$
// introduce symmetrized ?
== Coarse grid correction

The effect of the coarse-grid correction on the error is:

$
  e^((k+1)) &= (I - P(P^T A P)^(-1) P^T A) e^(k) \
  &= (I - Pi_c) e^(k)
$

$Pi_c : V arrow.r.long range(P)$ is the $A$-orthogonal projection onto $range(P)$

The coarse-grid correction removes the component of the error along range(P)

== Two-Grids Operator

A two-grid operator is an operator with two steps:
- Apply the subspace correction:
  $ e arrow.l (I - Pi_c) e $
- Apply the smoother:
  $ e arrow.l (I - R A)e $

The associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

== Subspace decomposition

We decompose the space into high and low frequency errors:

#align(center)[
  #fig-subspace-decomposition-SPD
]


// = Classic Convergence Theory

// == Oscillatory error

// An error is oscillatory when has a "high energy":\
// it has directions associated with large eigenvalues of A.
// #v(1cm)

// $ ||e||_A^2 = ( A e, e ) $

// $ ( A e, e ) = sum_k lambda_k c_k^2 $

// Huge if components are along the eigenvectors associated with large eigenvalues.

// == Smooth error

// An error is smooth when it cannot be reduced effectively by the smoother
// #v(1cm)
// $ S e &approx e \
//   arrow.double A e = r &approx 0
// $

// (or more rigorously when $||r||_(D^(-1)) << ||e||_A$)

// == informal 

// For the multigrid method to work we need two conditions:

// - The smoother is able to reduce high frequency errors

// - $range(P)$ is close to the space of low frequency error 

// == Smoothing property
// *The Smoothing Property*: The smoother $S$ reduces high-frequency errors.
// $ ||S e||_A^2 <= ||e||_A^2 - alpha ||A e||_(D^(-1))^2 $

// - $S$ is a contraction in the $A$-norm, and the relaxation converge.

// - The negative term quantifies the drop: 
//   - If the error is oscillatory, $||r||_(D^(-1))^2 = ||A e||_(D^(-1))^2 = sum_k lambda_k^2 c_k^2$ is huge and progress is fast.
//   - If the error is smooth, $A e approx 0$, and progress stalls.

// == The approximation property
// *The Approximation Property*: The coarse grid can approximate smooth errors.
// $ min_(u_c) ||e - P u_c||_D^2 <= beta ||e||_A^2 $

// == The approximation property
// *The Approximation Property*: The coarse grid can approximate smooth errors.
// $ min_(u_c) ||e - P u_c||_D^2 <= beta ||e||_A^2 $
// - If the error is oscillatory, $||e||_A$ is huge, and the inequality is easily satisfied, regardless of the quality of the interpolation.
// - If $e$ is a smooth error, the interpolation needs to accurately represent it: smooth errors must be near $op("range")(P)$.

// == Convergence theorem
// *Convergence theorem* \
// If the smoothing and the approximation property is true, and if $alpha < beta$,\ the iteration converges, with:
// $ ||S T G||_A <= sqrt(1 - alpha/beta) < 1 $

// #pagebreak()

// == The problems

// - developed in a geometric settings
// - In practice $alpha$ and $beta$ are out of reach
// - sufficient but not necessary conditions


= Subspace Correction Theory

// == The promise 

// - "In this paper, we try to develop a *unified framework* and theory that can be used to derive and analyze different algebraic multigrid methods in a coherent manner."

// - "Our theory applies to *most existing multigrid methods*, including the standard geometric multigrid method, the classic AMG, energy-minimization AMG, unsmoothed and smoothed aggregation AMG, and spectral AMGe."
// and they develop these method in their framework

// - they allegedly developed a unified theory ("subspace correction methods") applicable for domain decomposition and multigrid

== Informal Motivation

#align(center)[
  #fig-smoother-error
]

// #pagebreak()

// The vector $(I - Pi_c)(I - R A)e$ represents the residual error on $H$.

// The smaller $||(I - Pi_c)(I - R A)e||$ the better.

// #v(1cm)

// Informally, let's assume that $(I - R A)e_H in H$,\ meaning $(I - Pi_c)(I - R A)e_H = (I - R A)e_H$ \ 
// and let's assume that $R A$ is "symmetric".

// For a normalized error $||e_H|| = 1$, the maximum error is bounded by:

// $ sup_(e_H in H \ ||e_H|| = 1) ||(I - R A)e_H|| = lm(I - R A) = 1 - lmin(R A) $
// sup_(e_H in H \ ||e_H||_A = 1)||S e_H||

// == Two-Grids Operator

// A two-grid correction operator is an operator $B : V' -> V$ that acts on $f in V'$ with two steps:
// - Apply the subspace correction:
//   $ w = P A_c^(-1)P' f $
// - Apply the smoother:
//   $ B f = w + R(f - A w) $

// The associated error transfer operator is:
// $ E = (I - R A)(I - Pi_c) $

== Quantifying the smoother effect
#uncover("1-")[
  The effect of the smoother on any error $e in H$ is given by $S e = (I - R A)e$\

$
  ||S e||_A^2 &= ((I - R A) e, (I - R A) e)_A \
  &= ((I - R A)^*(I - R A) e, e)_A \
$

$(I - R A)^*$ the adjoint of $(I - R A)$ with respect to A: $(I - R A)^* = (I - R^T A)$

]

#uncover("2-")[
  $
  (I - R^T A)(I - R A) &= I - (R + R^T - R^T A R)A  \
    &= I - Rb A \
    "with" Rb = R& + R^T - R^T A R
  $
  acting as a symmetrized version of the smoother (equivalent to applying the smoother $R$ then its transpose $R^T$)
]
#uncover("3-")[
  $
    Rb A "is A-self adjoint:" quad (Rb A)^* = Rb A
  $
]

#pagebreak()

== Smoother convergence
$
  ||S||_A^2 &= max_(e in V) (||S e||_A^2) / (||e||_A^2) \
  #uncover("2-")[$&= max_(e in V) ((I - R A) e, (I - R A) e)_A / ((e,e)_A)$] \
  #uncover("3-")[$&= max_(e in V) (((I - Rb A)e,e)_A) / ((e,e)_A)$] \
  #uncover("4-")[$&= lm(I - Rb A), quad Rb A "is A-self adjoint"$] \
  #uncover("5-")[$&= 1 - lmin(Rb A)$]
$

#pagebreak()

#align(center)[
  #fig-X
]

$
  X = (I - Pi_c) Rb A : quad H^A -> H^A
$


== The operator X

with the same idea, for a general $v in V$, we have:

$
  ||E||_A^2 &= max_(v in V) (||(I - R A)(I - Pi_c) v||_A^2) / (||v||_A^2) \
  #uncover("2-")[$&= max_(v in V) (((I - Rb A)(I - Pi_c) v, (I - Pi_c) v)_A) / (||v||_A^2)$] \
  #uncover("3-")[$&= 1 - min_(v in V) ((Rb A (I - Pi_c) v, (I - Pi_c) v)_A) / (||(I - Pi_c) v||_A^2) quad "(not obvious)" $] \ // not obvious at all
  #uncover("4-")[$&= 1 - min_(e_H in H^A) (( Rb A e_H, e_H)_A) / (||e_H||_A^2)$]
$

#pagebreak()
$
  (Rb A e_H, e_H)_A = ((I-Pi_c) Rb A e_H, e_H)_A
$
// can probably be expressed more cleanly using the orthogonal decomposition
$Rb A e_H $, is composed of some part along $range(P)$ that will vanish when applying the scalar product with $e_H in H^A perp range(P)$: \
$ forall w in V, v in H^A, quad(w,v)_A &= ((I-Pi_c)w,v)_A + (Pi_c w,v)_A \
&= ((I-Pi_c)w,v)_A + 0 $

#pagebreak()

$
  ||E||_A^2 &= 1 - min_(e_H in H^A) (((I - Pi_c) Rb A e_H, e_H)_A) / (||e_H||_A^2) \
  &= 1 - lambda_(min)(X),
$

where $ X = (I - Pi_c) Rb A quad : quad H^A -> H^A $
// develop
$X$ represents the effect of the (symmetrized) smoother on the non-smooth errors.

#pagebreak()

#align(center)[
  #fig-X
]

$
  X = (I - Pi_c) Rb A : quad H^A -> H^A
$

#pagebreak()

// the effect of the two-grids cycle on $e in H$ is:
// $ overline(E)|_H &= (I - overline(T))( I - Pi_c ) \
//                  &= I - overline(T) quad (I - Pi_c = I "on" H) $

// The components along $H$ are given by:
// $ (I - Pi_c)overline(E)|_H &= I - (I - Pi_c)overline(T) \
//                             &= I - X $

// #pagebreak()

$ ||E||_A^2 = 1 - lambda_"min" (X) $

*Interpretation of $X$:* \
- If $X approx I$ on $H$, the smoother is effective and the error is almost entirely annihilated: \
  $lambda_"min" (X) approx 1$ and $||E||_A^2 = 1 - lambda_"min" (X) approx 0$

- If $X$ has a small eigenvalue, the associated direction in $H$ will escape both the effect of the smoother and the coarse space correction: convergence will stagnate.

== Optimal coarse space
// show => and assume it is optimal
$
  ||E||_A^2 = 1 - min_(e_H in H^A) (( Rb A e_H, e_H)_A) / (||e_H||_A^2)
$
// Orthogonal diagonalization of Rb A:
let $(mu_j, q_j)_(j=1)^n$ be the $A$-orthonormal eigenpairs of $Rb A$ ordered such that $0 < mu_1 <= dots <= mu_n$

Any high-frequency error $e_H in H^A$ can be expanded as:
$
  e_H = sum_(j=1)^n c_j q_j
$

// Rayleigh quotient expansion:
$
  (( Rb A e_H, e_H)_A) / (||e_H||_A^2) = (sum_(j=1)^n mu_j c_j^2) / (sum_(j=1)^n c_j^2)
$

#pagebreak()

$ ( Rb A e_H, e_H)_A = sum_(j=1)^n mu_j c_j^2$ is maximal when $H^A_("opt") = "span"{q_(n_c + 1), dots, q_n}$.

The minimum is attained at $e_H = q_(n_c + 1)$:
$
  min_(e_H in H^A_("opt")) (( Rb A e_H, e_H)_A) / (||e_H||_A^2) = mu_(n_c + 1) \
  quad arrow.r.double quad 
  ||E||_A^2 = 1 - mu_(n_c + 1)
$

= Applications

== The Ideal Smoother

#align(center)[
  #fig-optimal-smoother
]

#pagebreak()

An ideal smoother would make $(I - R A)$ the $A$-orthogonal projection onto $range(P)$.

If we write the eigendecomposition of $A$:
$ A = U mat(
  lambda_1(A), , , ;
  , lambda_2(A), , ;
  , , dots.down, ;
  , , , lambda_n(A)
) U^T $

Then the operators would take the form:
$ I - R A = U mat(
  1, , , , , ;
  , dots.down, , , , ;
  , , 1, , , ;
  , , , 0, , ;
  , , , , dots.down, ;
  , , , , , 0
) U^T 
quad "and" quad
R A = U mat(
  0, , , , , ;
  , dots.down, , , , ;
  , , 0, , , ;
  , , , 1, , ;
  , , , , dots.down, ;
  , , , , , 1
) U^T $

#pagebreak()

Under these conditions, the error transfer operator vanishes entirely:
$ E = (I - R A)(I - Pi_c) = 0 $

This is achieved when $R$ acts as the perfect inverse of $A$ on the high-frequency space $H$:

$ R = U mat(
  0, , , , , ;
  , dots.down, , , , ;
  , , 0, , , ;
  , , , 1/lambda_(n_(c+1)), , ;
  , , , , dots.down, ;
  , , , , , 1/lambda_n
) U^T $

and on $H$ we have indeed $lmin (R A) = 1$

== 1D Laplacian

For a 1D Laplacian operator A of size n,

we have:
$ A = mat(
  2,-1 , , ;
  -1 , 2, dots.down, ;
  ,dots.down , dots.down, -1 ;
  , , -1 , 2
)  quad
 = quad U mat(
  lambda_1(A), , , ;
  , lambda_2(A), , ;
  , , dots.down, ;
  , , , lambda_n (A)
) U^T $

// we have the eigenpairs for $ k in [| 1,n |]$:
// $
//   lk (A) = 4sin((k pi) / (2(n+1)))\
//   q_k (A) = ( sin((j k pi) / (n+1)) )_j
// $

== weighted Jacobi

$
  R = omega D^(-1) = omega /2 I 
$
for $omega = 2/3$,
$
  R = 1/3 I, quad S = I - 1/3 A \
  Rb =  1/3(2 I - 1/3 A)  
$
$A, R, S "and" Rb A$ share the same eigenvectors:
$
  lk (Rb A) = lk(Rb) lk(A) = 1/3lk(A)(2 - 1/3lk(A))
$

#pagebreak()

and we have:
$ A = U mat(
  lambda_1(A), , , ;
  , lambda_2(A), , ;
  , , dots.down, ;
  , , , lambda_n(A)
) U^T, quad
Rb A = U mat(
  lambda_1(Rb A), , , ;
  , lambda_2(Rb A), , ;
  , , dots.down, ;
  , , , lambda_n(Rb A)
) U^T $

to measure the action on $H^A$ we look at $ X = Rb A (I - Pi_c)$

$ I - Pi_c = U mat(
  0, , , , , ;
  , dots.down, , , , ;
  , , 0, , , ;
  , , , 1, , ;
  , , , , dots.down, ;
  , , , , , 1
) U^T , quad "for " range(P) = range(q_1, dots , q_(n_c) ) $

we obtain

$ X = U mat(
  0, , , , , ;
  , dots.down, , , , ;
  , , 0, , , ;
  , , , lambda_(n_c + 1)(Rb A), , ;
  , , ,  dots.down, ,;
  , , , ,lambda_n (Rb A)
) U^T quad X: V arrow V $

So 
$ X|_H = U_H mat(
   lambda_(n_c + 1)(Rb A), , ;
   dots.down , ;
   quad quad quad quad lambda_n (Rb A)
) U_H^T quad X|_H: H^A arrow H^A, quad U_H = U[ n_c + 1 : n] $

#pagebreak()

The two grid convergence is governed by
$
||E||_A^2 &= 1 - lmin(X)\ 
&= lambda_(n_c + 1)(Rb A)\
&= 1 - 1/3lambda_(n_c + 1)(A)(2 - 1/3lambda_(n_c + 1)(A))
$

= Numerical applications

== Effect of S

#align(center)[
  #fig-smoother-error
]

#pagebreak()

#align(center)[
 #image("../../pictures/plot_S_spectrum.png", width: 80%)
]

== Effect of $Rb A$
#align(center)[
 #image("../../pictures/plot_T_spectrum.png", width: 80%)
]

== Effect of X
#align(center)[
 #image("../../pictures/plot_X_opt_spectrum.png", width: 80%)
]

== two-grid convergence
#align(center)[
 #image("../../pictures/plot_convergence.png", width: 80%)
]


== Effect of X with P = linear interpolation
#align(center)[
 #image("../../pictures/plot_X_spectrum.png", width: 80%)
]

== Effect of X with P = linear interpolation
#align(center)[
 #image("../../pictures/plot_fourier_attenuation.png", width: 80%)
]

// numerical example

// Jacobi

// plot S = I - RA for Jacobi + GS + Additive Schwartz

// plot X for different P

// plot the theoretical convergence rate vs actual experiments

// supplemental: Z
// supplemental: fourier attenuation