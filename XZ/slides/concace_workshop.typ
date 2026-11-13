#import "@preview/touying:0.8.0": *
#import themes.metropolis: *
#import "@preview/cetz:0.5.1"
#import "concace_workshop_figures.typ": *

#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Elements of Multigrid],
    subtitle: [Smoothers, Coarse Projections, and Convergence],
    author: [Gabriel Clave],
    date: datetime.today().display(),
  )
)

// #show heading.where(level: 1): set text(size: em)

#let range = math.op("range")
#let nul = math.op("null")
#let im = math.op("Im")
#let ker = math.op("ker")
#let Mb = $overline(M)$
#let Tb = $overline(T)$
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let lk = $lambda_k$
#let lmin = $lambda_min$

#title-slide()

== Notations

$A in RR^(n times n)$ is SPD and $u,v,f in V = RR^n$
$
  "Problem:"          & quad A u = f \
  "Error:"            & quad e = u - v \
  "Residual:"         & quad r = f - A v \
  "Residual Problem:" & quad A e = r
$

= The Smoother

A smoother is a linear operator $M : V arrow.r.long V$ that can be used to solve the system through the iteration:
$
  v^((k+1)) = v^((k)) + M(f - A v^((k)))
$
with associated error propagation operator
$
  S = I - M A
$

the error iterate verifying:
$
  e^((k+1)) = S e^(k) \

  => e^((k+1)) = S^k e^(0)
$

The iteration is convergent if $rho(S) < 1$

== Jacobi Smoother convergence: relative error norm per iteration

#align(center)[
 #image("../../pictures/jacobi_convergence.png", width: 80%)
]
// smoother alone aren't great solver and stagnates

== Jacobi Smoother spectrum

#align(center)[
 #image("../../pictures/plot_J_spectrum.png", width: 80%)
]
// #align(center)[
//   #fig-spectrum-split(step: 1)
// ]

// the cause is eigenvalues of S close to 1

// "works well" / "not great" subspace
// we can decompose the space into H = "effective smoother space" and L = "ineffective smoother space"


== Smoother effect on the error

// convergence figure, one smoother step
// show big reduction in H, small reduction in S

#align(center)[
  #only("1")[#fig-smoother(step: 1)]
  #only("2")[#fig-smoother(step: 2)]
  #only("2")[#fig-smoother(step: 3)]
]

// S is cheap to compute, but effective only on some subspace H

= Coarse grid correction

Instead of solving $ A u = f quad "on" V $

Solve the smaller problem

$ A_c u_c = f_c quad "on a smaller" V_c $

//  correct the iterate using a smaller subspace $V_c$ with $ "dim"(V_c) < "dim"(V) $

- Transfer operators:
$ "Prolongation operator" P : V_c arrow.r.long V  \
"Restriction operator" P^T : V arrow.r.long V_c $
- Coarse Operator:  $A_c = P^T A P$

#pagebreak()

Coarse correction algorithm:

#align(center)[
  $
    #uncover("1-")[ $text("Compute residual:")       & quad r = A e$ ] \
    #uncover("2-")[ $text("Restrict residual:")      & quad r_c = P^T r$ ] \
    #uncover("3-")[ $text("Coarse solve:")     & quad e_c = A_c^(-1) r_c = (P^T A P)^(-1) r_c$ ] \
    #uncover("4-")[ $text("Interpolate & Correct:") & quad e arrow.l e - P e_c$ ]
  $
]

#v(0.5em)
#uncover("5-")[
  Combining these steps we get:
  $
    e arrow.l (I - P (P^T A P)^(-1) P^T A) e

  $
]

#pagebreak()

Coarse correction algorithm:
$
  e  arrow.l (I - P (P^T A P)^(-1) P^T A) e \
  e arrow.l (I - Pi_c)e
$

$
  Pi_c = P (P^T A P)^(-1) P^T A
$

$Pi_c$ is the $A$-orthogonal projection on $range(P)$

$=>$ The coarse correction removes all error component along $range(P)$

// #align(center)[
//   #only("1")[#fig-coarse-correction(step: 1)]
//   #only("2")[#fig-coarse-correction(step: 2)]
// ]

= Basic multigrid

== Core idea (text)

*Fundamental Principle:*\
The central objective when designing a multigrid method is to construct $P$ and $S$ such that $range(P)$ and $H$ are *complementary*:
$ V approx range(P) xor_A H $

Ensuring that error components missed by the smoother are effectively captured and eliminated by the coarse grid correction.

== Core idea (picture)

// two grid iteration figure with "close to orthogonal" range(P) and H
#align(center)[
  #only("1")[#fig-twogrid(step: 1)]
  #only("2")[#fig-twogrid(step: 2)]
  #only("2")[#fig-twogrid(step: 3)]
  #only("2")[#fig-twogrid(step: 4)]
]

== A multigrid solver

A multigrid solver relies on three core ingredients:

- A *Smoother $M$* eliminates error in a subspace $H$.

- An *Interpolation $P$* & the associated Coarse Grid

- A *Coarse Solver* (can recursively be another multigrid cycle).

== Recursive multigrid: V-cycle

#align(center)[
  #fig-v-cycle
]
// a word on recursive solver => multigrid cycle

= Convergence of two-grids iteration

== Two-Grids Operator

A two-grid operator is an operator with two steps:
- Apply the subspace correction:
  $ e arrow.l (I - Pi_c) e $
- Apply the smoother:
  $ e arrow.l (I - M A)e $

The associated error transfer operator is:
$ E = (I - M A)(I - Pi_c) $

== Informal Motivation

// Introducing X
#align(center)[
  #only("1")[#fig-smoother-convergence(step: 1)]
  #only("2")[#fig-smoother-convergence(step: 2)]
]

== Two-grid convergence
// "trust me bro" version

// uncover
To quantify the convergence we want to compute:
$
  ||E||_A = ||(I - M A)(I - Pi_c)||_A
$

To do so we study a slightly different operator:

$
  X = (I - Pi_c) Mb A e_H : quad H^A -> H^A
$

$ "with" quad Mb = M + M^T - M^T A M $

// The worst possible case for convergence is the one that minimize the effect of $X$, \ ie maximize the norm of $(I - X)e_H$ :
// $
//   ||E||_A^2 &= sup_(||e_H||_A = 1)||(I - X)e_H||_A^2 \
//   &= lm(I - X) \
//   &= 1 - lmin(X)
// $

We can show that:
$
  ||E||_A^2 = 1 - lmin(X)
$
#pagebreak()

$
  X = (I - Pi_c) Mb A : quad H^A -> H^A
$


$ ||E||_A^2 = 1 - lambda_"min" (X) $

*Interpretation of $X$:* \
- If $X approx I$ on $H$, the smoother is well-matched to the coarse space, and the error is almost entirely eliminated: \
  $ lambda_"min" (X) approx 1 => ||E||_A^2 = 1 - lambda_"min" (X) approx 0 $

- If $X$ has a small eigenvalue, the associated direction in $H$ will escape both the effect of the smoother and the coarse space correction: convergence will stagnate.

== Designing $P$

For a given problem $(A, V)$ and a given smoother $M$, how do we construct $P$ ?

#v(1em)

#only("1")[

  *Theoretical Ideal:*

  Theory yields an optimal $P_"opt"$ maximizing $lmin(X)$, \ built from the lowest modes of a global eigenvalue problem:

  $ Mb A v = lambda v $

  #v(2em)

  $arrow.r.double$ Impractical: more expensive than solving $A u = f$ !
]

#only("2")[
  #grid(
    columns: (1fr, 1fr),
    gutter: 1.5em,
    [
      *Geometric Multigrid (GMG)*
      - $P$ is built directly from the *mesh, geometry, or physics* (e.g., linear interpolation).

      - Very cheap setup, fast execution.

      - *Problem-dependent:* requires geometric hierarchy and domain expertise.
    ],
    [
      *Algebraic Multigrid (AMG)*
      - $P$ is built purely from the matrix *$A$* and the smoother *$M$*.

      - Black-box / general-purpose (unstructured grids, complex topologies).

      - *Setup overhead:* adds operator construction costs.
    ]
  )
]

= Numerical applications

== Spectrum of classic multigrid smoother

#align(center)[
 #image("../../pictures/S_Spectrum_true_GS.png", width: 90%)
]

== Spectrum of $X = (I - Pi_c) Mb A$ for classic multigrid smoother
#align(center)[
 #image("../../pictures/plot_X_opt_spectrum_trunc.png", width: 90%)
]

// example of GMG
== Gauss-Seidel and linear interpolation Prolongation
#align(center)[
 #image("../../pictures/plot_X_spectrum_GS_gmg.png", width: 90%)
]

// example of AMG heuristic
== Gauss-Seidel and AMG Prolongation
#align(center)[
 #image("../../pictures/plot_X_spectrum_GS_sa.png", width: 90%)
]

#focus-slide[
  #text(size: 1.8em, weight: "medium")[Thank you for your attention!]

  #v(1.5em)

  #text(size: 0.7em, fill: white.darken(15%))[
    *Main Reference:* \
    Jinchao Xu and Ludmil Zikatanov, _"Algebraic multigrid methods"_, \
    Acta Numerica, vol. 26, pp. 591–721, Cambridge University Press, 2017. \
  ]
]