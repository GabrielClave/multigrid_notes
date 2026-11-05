#import "@preview/touying:0.8.0": *
#import themes.metropolis: *
#import "@preview/cetz:0.5.1"
#import "concace_workshop_figures.typ": *

#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Elements of Multigrid],
    author: [Gabriel Clave],
    date: datetime.today().display(),
  )
)

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

#title-slide()

== Notations

= Smoother

The first component of a multigrid solver is a smoother.\
It is a linear operator $B : V arrow.r.long V$ that can be used to solve the system through the iteration:
$
  v^((k+1)) = v^((k)) + B(f - A v^((k)))
$
with associated error propagation operator
$
  S = I - B A
$

the error iterate verifying:
$
  e^((k+1)) = S e^(k)
$

== Smoother convergence

#align(center)[
 #image("../../pictures/jacobi_convergence.png", width: 80%)
]
// smoother alone aren't great solver and stagnates

== Smoother spectrum
#align(center)[
 #image("../../pictures/plot_S_spectrum.png", width: 80%)
]
// the cause is eigenvalues of S close to 1

// "works well" / "not great" subspace
// we can decompose the space into H = "effective smoother space" and L = "ineffective smoother space"


== title

// convergence figure, one smoother step
// show big reduction in H, small reduction in S

#align(center)[
  #only("1")[#fig-smoother(step: 1)]
  #only("2")[#fig-smoother(step: 2)]
  #only("2")[#fig-smoother(step: 3)]
]

// S is cheap to compute, but effective only on some subspace H

= Coarse grid correction

// subspace correction:
// instead of solving Au = f on V, we solve it on a smaller V_c dim(V_c) < dim(V)
// ingredients:
// transfer operator P
// coarse solver A_c^(-1) //(or local operator A_c P^TAP directly ?)
// algorithm:
// - Compute local residual: $r_(i-1)^("local") = P_i^T A e_(i-1) in bb(R)^(n_i)$.
// - Solve local subproblem: $x_i = A_i^(-1) r_(i-1)^("local") in bb(R)^(n_i)$.
// - Update global iterate:
// $ 
//   e_i &= e_(i-1) - P_i x_i \
//   &= (I - P_i A_i^(-1) P_i^T A ) e_(i-1)
// $

// We introduce the subspace iteration matrix $T_i$:
// $ T_i = P_i A_i^(-1) P_i^T A in bb(R)^(n times n) $

// The local update satisfies:
// $ e_i = (I - T_i) e_(i-1) $
// If A_c = P^TAP -> A-orthogonal projection
// 
// convergence figure, coarse grid correction step
// The coarse-grid correction removes the component of the error along range(P)

= Basic Multigrid

== Core idea (text)

*Fundamental Principle:*\
The central objective when designing a multigrid method is to construct $P$ and $S$ such that $range(P)$ and $H$ are *complementary*:
$ V approx range(P) xor_A H $

Ensuring that error components missed by the smoother are effectively captured and eliminated by the coarse grid correction.

== Core idea (picture)

// two grid iteration figure with "close to orthogonal" range(P) and H

// a word on recursive solver => multigrid cycle

= Convergence of two-grid iteration

== Two-Grids Operator

A two-grid operator is an operator with two steps:
- Apply the subspace correction:
  $ e arrow.l (I - Pi_c) e $
- Apply the smoother:
  $ e arrow.l (I - R A)e $

The associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

== Informal Motivation

#align(center)[
  #only("1")[#fig-smoother-convergence(step: 1)]
  #only("2")[#fig-smoother-convergence(step: 2)]
]

== Two-grid convergence
// "trust me bro" version

// actually we don't study X = (I - Pi_c) B A e_H but X = (I - Pi_c) Bb A e_H

// The worst possible case for convergence is the one that minimize the effect of X, ie maximize the norm of (I - X)e_H

// ||E||_A^2 = sup ||(I - X)e_H|| / ||e_H|| = lmax(I - X) = 1 - lmin(X)


// full difficult version:

// first introduce Rb somehow
// == Quantifying the smoother effect
// #uncover("1-")[
//   The effect of the smoother on any error $e in H$ is given by $S e = (I - R A)e$\

// $
//   ||S e||_A^2 &= ((I - R A) e, (I - R A) e)_A \
//   &= ((I - R A)^*(I - R A) e, e)_A \
// $

// $(I - R A)^*$ the adjoint of $(I - R A)$ with respect to $A$
// ]
// introduce Rb somehow "we can show that (I - R A)^*(I - R A) = (I - Rb A)"

// == Smoother convergence
// $
//   ||S||_A^2 &= max_(e in V) (||S e||_A^2) / (||e||_A^2) \
//   #uncover("2-")[$&= max_(e in V) ((I - R A) e, (I - R A) e)_A / ((e,e)_A)$] \
//   #uncover("3-")[$&= max_(e in V) (((I - Rb A)e,e)_A) / ((e,e)_A)$] \
//   #uncover("4-")[$&= lm(I - Rb A), quad Rb A "is A-self adjoint"$] \
//   #uncover("5-")[$&= 1 - lmin(Rb A)$] // some hypothesis needed
// $

// #pagebreak()

// #align(center)[
//   #fig-X
// ]

$
  X = (I - Pi_c) Rb A : quad H^A -> H^A
$


$ ||E||_A^2 = 1 - lambda_"min" (X) $

*Interpretation of $X$:* \
- If $X approx I$ on $H$, the smoother is effective and the error is almost entirely annihilated: \
  $lambda_"min" (X) approx 1$ and $||E||_A^2 = 1 - lambda_"min" (X) approx 0$

- If $X$ has a small eigenvalue, the associated direction in $H$ will escape both the effect of the smoother and the coarse space correction: convergence will stagnate.


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
