#import "@preview/cetz:0.5.1"

// Document setup & page styling
#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  numbering: "1",
)

// Font configuration (using XCharter / Charter equivalent)
#set text(
  font: "xCharter",
  size: 11pt,
  lang: "en",
)

// Paragraph spacing
#set par(justify: true, leading: 0.65em)

// Custom Math Operators
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
#let tb(x) = $tilde(bold(#x))$

// Document Title & Metadata
#align(center)[
  #text(size: 18pt, weight: "bold")[Abstract Multigrid notes] \
  #v(8pt)
  #text(size: 12pt)[Gabriel Clave] \
  #v(8pt)
  #datetime.today().display("[month repr:long] [day], [year]")
]

#v(1.5em)

// Table of Contents
#outline(indent: auto)

#v(2em)

= Variational Formulation

We consider the problem
$ A u = f $
with $A: V -> V'$ SSPD, $u in V$ and $f in V'$.

$u$ and $f$ do not live in the same vector spaces:

- $u$ is in $V tilde.equiv RR^n$, and represents a primal state (temperature, pressure), and the natural norm is the energy norm $\|x\|^2 = x^T A x$.
- $f$ is in $V' tilde.equiv (RR^n)'$, and represents a source (force, flux), and the natural norm is the dual norm $\|b\|^2 = b^T A^(-1) b$.

// In particular x^T b does not really make sense

It is this problem that arises naturally from the weak formulation.

We solve
$ "Find " u in V " such that " chevron.l A u, v chevron.r = chevron.l f, v chevron.r quad forall v in V $

// Noting that Au is in V', so <Au, v> is the real number (Au)(v)

== Smoother

=== general linear iteration

The first component of a multigrid solver is a smoother.\
It is a linear operator $R : V' arrow.r.long V$ that can be used to solve (equation) through the linear iteration:
$
  v^((k+1)) = v^((k)) + R(f - A v^((k)))
$

The effect of the smoother on any error $e in V$ is given by $S e = (I - R A)e$\

// when $R$ is not symmetric, $S$ is not symmetric either:
$
  ||S e||_A^2 &= ((I - R A) e, (I - R A) e)_A \
  &= ((I - R A)^*(I - R A) e, e)_A \
$

$(I - R A)^*$ the adjoint of $(I - R A)$ with respect to A:
$
 (I - R A)^* &= A^(-1)(I - R A)^T A\
 & = (I - R^T A)
$
which correspond to the error propagation operator of the smoother $R^T$.

we have:
$
 (I - R^T A)(I - R A) &= I - (R + R^T - R^T A R)A  \
  &= I - Rb A \
  "with" Rb = R& + R^T - R^T A R
$
$Rb$ is a symmetric operator, acting as a symmetrized version of the smoother. \
$Rb$ is equivalent to applying the smoother $R$ then its transpose $R^T$:

$
cases(
  u^(m - 1/2) &= u^(m - 1) + R(f - A u^(m - 1)),
  u^m &= u^(m - 1/2) + R'(f - A u^(m - 1/2))
)
$

with associated error propagation operator
$
  overline(S) = I - Rb A \
$

// The symmetrized iteration converges if $rho(overline(S)) <1$,\
// which implies $||S||_A = sqrt(rho(overline(S))) <1$ and so regular iteration is also convergent.

=== A multigrid smoother

Let $ R : V' arrow.r.long V$ be a valid smoother. \
If $e in ker(A)$, then $S e = (I - R A)e = e$: nothing happens.\
The smoother alone is unable to reduce the composant of the error that are in ker(A)

similarly, if $e_lambda in E_lambda (A)$, then $S e_lambda = (I - lambda R)e_lambda$\
for a given $R$ with bounded norm,  $S e_lambda arrow.r.long_(lambda arrow.r 0) e_lambda $

It will be difficult for the smoother to remove composant of the error in directions associated with the smallest eigenvalues of $A$.

The role of R is to act as an approximate inverse of A, but not on the entire spectrum.\
The application of the smoother will damp error composants in the directions associated with the highest eigenvalues of $A$.\
We will call this vector space $H^A$, the space of high frequency error.\
The error not in $H^A$ will need to be addressed separately.

==== smoother geometry

$Ri$ is invertible by design and symmetric: it defines a scalar product on $V$ $(dot,dot)_(Ri).$

Let $Tb = Rb A$

We have $ker(A) subset ker(Tb)$

$Tb$ acts as a bridge between the $A$ geometry and the $Ri$ geometry
$ 
  (u,v)_A = (Rb A u, v)_(Ri) = (Tb u,v)_(Ri)
$

we also have:

$
  x in ker(Tb) & arrow.r.double ||Tb x||_(Ri) = 0 \
                         & arrow.r.double ||x||_A = 0 \
                         & arrow.r.double x in ker(A)
$

which means $ker(A) = ker(Tb)$ and so $range(Tb) = W$
// W is not defined yet

== Coarse space correction

In order to deal with the error composant unaffected by the smoother, the idea is to use a coarse space $V_c$, linked with V by the injective operator $P: V_c -> V$

We will then look to approximate the solution of (equation) in the range of $P$, which means finding a $u_P in range(P)$, such that $u approx u_P = P u_c, quad u_c in V_c$.
// need smoother link: Galerkin projection allows to remove composant in range(P)

We are doing a Galerkin projection: instead of solving the variational problem on the entire space V, we restrict the search space to $range(P) in V$:

Find $u_P in range(P)$ such that
$
                        &chevron.l A u_P, v_P chevron.r = chevron.l f, v_P chevron.r quad forall v_P in range(P) \
  arrow.l.r.double.long &chevron.l A P u_c, P v_c chevron.r = chevron.l f, P v_c chevron.r quad  forall v_c in V_c \
  arrow.l.r.double.long &chevron.l A P u_c, P v_c chevron.r = chevron.l f, P v_c chevron.r quad forall v_c in V_c \
  arrow.l.r.double.long &chevron.l P' A P u_c, v_c chevron.r = chevron.l P' f, v_c chevron.r quad forall v_c in V_c
$

Noting $A_c = P' A P$ and $f_c = P'f$, we are solving the equivalent system on $V_c$
$ A_c u_c = f_c $

If we define $e_("approx") = u - u_P$, we can note that $chevron.l A e_("approx"), v_P chevron.r = 0$ for all $v_P in range(P)$, i.e., the approximation error is orthogonal to $range(P)$.

// a word on solving the coarse problem, and designing an AMG method

= Convergence Theory

Let us consider $A$ SSPD. \
// We will consider the case $V_c subset V$.\
A has a non trivial kernel $N$.\
The smoother $E = I - R A$ is blind to anything in $N$, so the coarse space correction has to deal with the error in $N$: we need
$
 N subset range(P)
$

Using W = $N^(perp_A)$ we decompose $V$ into
$
  V = N xor_A W 
$

and we define $ W_P = range(P) inter W$\
meaning we have
$
  range(P) = W_P xor_A N = P(V_c)
$

$P$ is injective from $V_c$ to $range(P)$, so $V_c$ can be decomposed into
$
  V_c = W_c xor_A N_c
$
with $W_c = P^(-1)(W_P)$ and $N_c = P^(-1)(N)$

#align(center)[
  #cetz.canvas({
    import cetz.draw: *

    let x0 = 0.0
    let xN = 4.0
    let xP = xN + 2
    let xend = 12.0

    let y0 = 0.0
    let yend = 2.0
    let margin = 0.25

    // Custom palette
    let main-purple = rgb("#7570b3")
    let light-purple = rgb("#c1bfda").transparentize(65%)
    let red-accent  = rgb("#e04040")
    let gray-bg     = rgb("#f0f0f0")

    let brace-margin = 6.5*margin

    set-style(
      content: (padding: 0.1),
      text: (size: 20pt)
    )

    // Subspace W_P = Range(P) \cap W
    rect((xN, y0), (xP, yend), fill: light-purple, stroke: 1pt)
    content(((xN + xP) / 2, (y0 + yend) / 2), $W_P$)

    // Entire space V
    rect((x0, y0), (xend, yend), stroke: 1pt, name: "V_rect")

    // Subspace N = Ker(A)
    rect((x0, y0), (xN, yend), fill: gray-bg, stroke: 1pt)
    content(((x0 + xN) / 2, (y0 + yend) / 2), $N = "Ker"(A)$)

    // Complementary space H^A
    content(((xP + xend) / 2, (y0 + yend) / 2), $H^A = "range"(P)^(perp_A)$)

    // Outer dashed box for Range(P)
    rect(
      (x0 - margin, y0 - margin),
      (xP, yend + margin),
      stroke: (dash: "dashed", paint: red-accent, thickness: 1pt),
      name: "rangeP"
    )

    // Label for Range(P) anchored north-east above the box
    content(
      (rel: (0, margin), to: "rangeP.north-east"),
      [#text(fill: red-accent, $"range"(P)$)],
      anchor: "south-east"
    )

    // Double-headed arrow / edge-line for W = N^\perp_A below the main drawing
    let yW = y0 - (2 * margin)

    line(
      (xN, yW),
      (xend, yW),
      mark: (start: "bar", end: "bar"),
      stroke: 0.8pt,
      name: "lineW"
    )

    // Label for W centered underneath the line
    content(
      (rel: (0, -margin), to: "lineW.mid"),
      [#text($W = N^(perp_A)$)],
      anchor: "north"
    )

    // low freq brace
    cetz.decorations.brace((xP - 0.25*margin, y0 - brace-margin), (x0, y0 - brace-margin), name: "b_low")
    content(
      (rel: (0, -margin), to: "b_low.center"),
      [#text("low frequency error")],
      anchor: "north"
    )

    // high freq brace
    cetz.decorations.brace((xend, y0 - brace-margin), (xP + 0.25*margin, y0 - brace-margin), name: "b_low")
    content(
      (rel: (0, -margin), to: "b_low.center"),
      [#text("High frequency error")],
      anchor: "north"
    )
  })
]
#v(1cm)

We define the following $A$-orthogonal projections:

$
  Pi_c : V &arrow.r.long range(P)\
  Pi_W : V &arrow.r.long W
$

and we will note that when $w in W_P = range(P) inter W$, we have $Pi_c (w) = Pi_W (w)$

#v(1cm)

We define the following $Ri$-orthogonal projections:

$
 Q_c: V &arrow.r.long range(P)\
 Q_W: V &arrow.r.long W
$
#v(1cm)
=== Two-Grids Operator

A two-grid correction operator is an operator $B : V' arrow.r.long V$ that acts on $f in V'$ with the two steps:
- apply the subspace correction
$
  w &= P A_c^(-1)P' f // f or r ?
$
- apply the smoother
$
  B f &= w + R(f - A w)
$

and the associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

the convergence rate of the method is given by:
$ ||E||_A^2 = 1 - 1 / (K (V_c)) $

with
$
  K (V_c) = max_(w in W) (|| (I - Q_c)v||_(Ri)^2 ) / ( ||v||_A^2) =
  max_(w in W) min_(v_P in range(P)) (|| w - v_P||_(Ri)^2 ) / ( ||v||_A^2)
$
#v(1cm)
=== Proof

For all $v in V$
$
  ||(I - R A) v||_A^2
  &= ((I - R A) v, (I - R A) v)_A \
  &= ((I - R A)^*(I - R A) v, v)_A \
  &= ((I - Rb A) v, v)_A \
$

We then have:
$
  ||E||_A^2 &= max_(w in W) (||(I - R A)(I - Pi_c) w||_A^2) / (||w||_A^2) \
  &= max_(w in W) (((I - Rb A)(I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) \
  &= 1 - min_(w in W) ((Rb A (I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) \
  &= 1 - min_(w in W) ((Rb A (I - Pi_c) w, (I - Pi_c) w)_A) / (||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2) \
  &= 1 - min_(v in H^A) (( Rb A v, v)_A) / (||v||_A^2) \
  &= 1 - min_(v in H^A) (((I - Pi_c) Rb A v, v)_A) / (||v||_A^2) \
  &= 1 - lambda_(min)(X),
$

where $ X = (I - Pi_c) Rb A quad : quad H^A arrow.r.long H^A $

#v(1cm)

additional details:

$
  ((Rb A (I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) = ((Rb A (I - Pi_c) w, (I - Pi_c) w)_A) / (||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2)
$

// $Rb A (I - Pi_c) w in range(Tb) = W$ so it is unchanged by $Q_W$ \
$||w||_A^2 = ||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2$ using the orthogonal decomposition of W

#v(1cm)

$
  (Rb A v, v)_A = ((I-Pi_c) Rb A v, v)_A
$
// can probably be expressed more cleanly using the orthogonal decomposition
$Rb A v in range(Rb A) = W$, so it is composed of some part along $W_P$ that will vanish when applying the scalar product with $v in H^A perp W_P$: \
$ forall w in W, v in H^A, quad(w,v)_A &= ((I-Pi_c)w,v)_A + (Pi_c w,v)_A \
&= ((I-Pi_c)w,v)_A + 0 $

// #align(center)[
//   #cetz.canvas({
//     import cetz.draw: *

//     let circle_radius = 0.06
//     let main-purple = rgb("#7570b3")
//     let red-accent  = rgb("#e04040")
//     let green-accent = rgb("#1b9e77")

//     // Coordinates
//     let orig = (0, 0)
//     let w-pt = (3.5, 3.2)           // Vector w in W
//     let w-P = (3.5, 0)             // Projection onto W_P (Pi_c w)
//     let w-H = (0, 3.2)             // Projection onto W_P^\perp_A ((I - Pi_c) w)
//     let v-P = (5, 0)             // Test vector v_P in W_P

//     // Axes
//     line((-0.5, 0), (5.5, 0), stroke: 1pt, name: "x-axis")
//     line((0, -0.5), (0, 4.2), stroke: 1pt, name: "y-axis")

//     content((5.6, 0), $W_P$, anchor: "west")
//     content((0, 4.3), $W_P^perp_A$, anchor: "south")

//     // --- DECOMPOSITION OF w ---
//     // Dashed projection lines for w
//     line(w-pt, w-H, stroke: (dash: "dashed", paint: gray))
//     line(w-pt, w-P, stroke: (dash: "dashed", paint: gray))

//     // Vector w
//     line(
//       orig, w-pt,
//       mark: (end: "triangle", fill: main-purple),
//       stroke: (paint: main-purple, thickness: 1.5pt)
//     )
//     circle(w-pt, radius: circle_radius, fill: black, stroke: none)
//     content(w-pt, padding: 0.15, $w$, anchor: "south-west")

//     // Component w_H = (I - Pi_c)w on y-axis
//     content(w-H, padding: 0.15, $w_H$, anchor: "east")

//     // Component w_P = Pi_c w on x-axis
//     content(w-P, padding: 0.15, $w_P$, anchor: "north")

//     // --- TEST VECTOR v_P AND PRODUCT EQUIVALENCE ---
//     // Vector v_P along W_P
//     line(
//       orig, v-P,
//       mark: (end: "triangle", fill: green-accent),
//       stroke: (paint: green-accent, thickness: 1.5pt),
//       name: "v_vector"
//     )
//     content(v-P, padding: 0.15, $v_P$, anchor: "north")

//     // Annotation for inner product equality
//     content(
//       (rel: (-2, 0.25), to: "v_vector.mid") ,
//       text(size: 8.5pt)[
//         $(w, v_P)_A = (w_P, v_P)_A$
//       ],
//       anchor: "west"
//     )
//   })
// ]

#v(1cm)

==== The operator X

$
  X = (I - Pi_c) Tb quad : quad H^A arrow.r.long H^A
$

X is self-adjoint with respect to $(dot, dot)_A$, as $I-Pi_c$ is an $A$-orthogonal projection and $Rb$ is SPD.\
// The $Q_W$ projection is redundant here as in our setting $im(Tb) = W$

X represents the effect of the (symmetrized) smoother on the non-smooth errors.

// this one is weird
Let $H^A = range(P)^(perp_A)$ be the space of these high frequency errors.\
We have $H^A = W_P^(perp_A) inter W$, and the effect of the two-grids cycle on $e in H$ is:
$
  overline(E)|_(H^A) &= (I - Tb)( I - Pi_c ) \
                 &= I - Tb quad (I - Pi_c = I "on" H)\
$

the composant along H are given by:
$
  (I- Pi_c)overline(E)|_(H^A) &= I - (I- Pi_c)Tb \
                          &= I - X
$

if $X approx I "on" H^A$, the smoother is effective and the error is almost entirely annihilated: \
$lambda_"min"(X) approx 1$ and $||E||_A^2 = 1 - lambda_"min"(X) approx 0$

if X has a small eigenvalue, the associated direction in H will escape both the effect of the smoother and the coarse space correction: converge will stagnate.

#align(center)[
  #cetz.canvas({
    import cetz.draw: *

    let circle_radius = 0.06
    let margin = 1
    let main-purple = rgb("#7570b3")
    let red-accent  = rgb("#e04040")

    // Define main points
    let orig = (0, 0)
    let e-pt = (4.5, 3.5)          // Initial error vector e
    let e-H = (0, 3.5)            // High-frequency component e_H
    let e-P = (4.5, 0)            // Coarse-space component e_P
    let e-smoothed = (3, 0.9)   // Error after smoother: (I - Tb) e
    let e-H-smoothed = (0, 0.9) // Reduced high-frequency component (I - X) e_H
    let e-P-smoothed = (3, 0) // Slightly reduced coarse component

    // Set coordinate system and axes
    line((-0.5, 0), (5.5, 0), stroke: 1pt, name: "x-axis")
    line((0, -0.5), (0, 4.5), stroke: 1pt, name: "y-axis")

    content((5.6, 0), $text(range)(P)$, anchor: "west")
    content((0, 4.6), $H^A$, anchor: "south")

    // --- INITIAL STATE e ---
    // Dashed projection lines for e
    line(e-pt, e-H, stroke: (dash: "dashed", paint: gray))
    line(e-pt, e-P, stroke: (dash: "dashed", paint: gray))

    // Projections of e on axes
    // circle(e-H, radius: 0.06, fill: blue, stroke: none)
    content(e-H, padding: 0.15, $e_H$, anchor: "east")
    // circle(e-P, radius: 0.06, fill: blue, stroke: none)
    content(e-P, padding: 0.15, $e_P$, anchor: "north")

    // --- STATE AFTER SMOOTHER (I - Tb)e ---
    // Dashed projection lines for smoothed error
    line(e-smoothed, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))

    // Vector (I - Tb)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - Tb)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(e-H-smoothed, padding: 0.15, $(I - Pi_c)(I - Tb)e_H\ =(I - X)e_H$, anchor: "east")
    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H^A component (reduction arrow)
    line((0, 3.5), (0, 0.9), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1pt))
    content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })
]

==== The operator Z

The inverse of $X$ on $W_P^(perp_A)$ can be explicitly written as:
$
  Z = (Tb)^(-1) (I - Q_c) quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A) 
$

using $(u,v)_A = (Tb u,v)_(Ri)$ we have indeed, $forall u in W, forall v_P in W_P$:

$
  (Z u, v)_A &= ( (I - Q_c) u, v )_(Ri) = 0
$
meaning $im(Z) subset W_P^(perp_A)$, and we have// don't we need the equality here ?

$
  X Z = (I - Pi_c)(I - Q_c) = I "on" W_P^(perp_A)
$

$"so" lambda_min (X) = 1 / (lambda_max (Z)) $

interpretation of $Z$:\
after you apply the smoother the remaining smooth error should be as close as possible to $range(P)$\
the distance is $ d_Ri (e,range(P)) = min_(e_P in range(P)) ||e - e_P||_Ri = ||(I - Q_c)v||^2_Ri$

#align(center)[
  #cetz.canvas({
    import cetz.draw: *

    let circle_radius = 0.06
    let main-purple = rgb("#7570b3")
    let red-accent  = rgb("#e04040")
    let blue-accent = rgb("#386cb0")

    let vp = 4.2
    let vh = 3.2

    // Coordinates
    let orig = (0, 0)
    let v-pt = (vp, vh)            // Vector v in V
    let Qc-v = (vp, 0)             // Projection Q_c v on range(P)
    let IQc-v = (0, vh)            // Projection (I - Q_c) v on W_P^\perp_A

    // Axes
    line((-0.5, 0), (5.5, 0), stroke: 1pt, name: "x-axis")
    line((0, -0.5), (0, 4.2), stroke: 1pt, name: "y-axis")

    content((5.6, 0), $text(range)(P)$, anchor: "west")
    content((0, 4.3), $W_P^perp_Ri$, anchor: "south")

    // --- PROJECTIONS OF v ---
    // Dashed lines forming the rectangular decomposition
    line(v-pt, IQc-v, stroke: (dash: "dashed", paint: gray))
    line(v-pt, Qc-v, stroke: (dash: "dashed", paint: gray))

    // Vector v
    line(
      orig, v-pt,
      mark: (end: "triangle", fill: main-purple),
      stroke: (paint: main-purple, thickness: 1.5pt)
    )
    circle(v-pt, radius: circle_radius, fill: black, stroke: none)
    content(v-pt, padding: 0.15, $v$, anchor: "south-west")

    // Coarse component Q_c v on range(P)
    content(Qc-v, padding: 0.15, $Q_c v$, anchor: "north")

    // Complement component (I - Q_c)v on y-axis
    content(IQc-v, padding: 0.15, $(I - Q_c)v$, anchor: "east")

    // --- DISTANCE ANNOTATION (Double-edged arrow along y-axis) ---
    line(
      (vp + 0.15, 0.1), (vp + 0.15, 3.1),
      mark: (start: "triangle", end: "triangle", fill: red-accent),
      stroke: (paint: red-accent, thickness: 0.6pt),
      name: "distance"
    )

    content(
      (rel: (1.5, 0), to: "distance"),
      text(size: 10pt)[
        $d_Ri (v, text(range)(P)) \ = ||(I - Q_c)v||_Ri$
      ],
      anchor: "center"
    )
  })
]

$Z$ allows to express this distance in the $A$-geometry:
$
  (Z v, v)_A = (Tb^(-1)(I - Q_c)v, v)_A &= (Tb Tb^(-1)(I - Q_c)v, v)_Ri \
                                    &= ||(I - Q_c)v||_Ri^2 
$
and so in the worst case scenario we have:
$
  lambda_max (Z) &= max_(v in W_P^(perp_A)) r_Z (v) =  max_(v in W_P^(perp_A))( (Z v, v)_A ) / (v,v)_A\

                 &= max_(v in W_P^(perp_A)) (||(I - Q_c)v||_Ri^2) / (||v||_A^2 )
$

which is also the maximum for $w in W,"decomposed into" w = v + v_P$, $v in W_P^(perp_A) "and" v_P in W_P$
$
  r_Z (w) = (||(I - Q_c)(v + v_P)||_Ri^2) / (||v||_A^2 + ||v_P||_A^2) lt.eq (||(I - Q_c)v||_Ri^2) / (||v||_A^2)
$
$ "so " lambda_max (Z)= max_(w in W) r_Z (w) = max_(v in W_P^(perp_A)) r_Z (v) = K(V_c) $

$ ||E||_A^2 = 1 - 1 / (K (V_c)) $

which can be understood as follow: the effectiveness of a two-grids methods depends on 
how well the coarse space can approximate the smooth error (quantified by $Z$) that escaped the smoother (quantified by $X$)

=== Optimal Coarse Space

We want to chose the optimal coarse space $V_c$ so that $range(V_c)$ is the best complement to the space of high frequency errors H.

This is achieved when $K(V_c)$ is minimal.
$
  1 / (K(V_c)) = min_(w in W) max_(w_P in W_P) (||w||_A^2)/(||w - w_P||_Ri^2)
$
If we take the minimum over the space of high frequency error $H^Ri = (W_P^(perp_Ri) inter W) subset W$,\ for which $min_(w in H^Ri) ||w - w_P||_Ri^2 = ||w||_Ri^2 $ // dessin ?
$
  1 / (K(V_c)) &lt.eq min_(w in H^Ri) max_(w_P in W_P) (||w||_A^2)/(||w - w_P||_Ri^2) \
  &= min_(w in H^Ri)(||w||_A^2)/(||w||_Ri^2) \
  &= min_(w in H^Ri)(||Rb A w||_Ri^2)/(||w||_Ri^2) \
  &= min_(w in H^Ri) r_(Rb A) (w)
$

and so
$
  max_(dim V_c = n_c) 1 / (K(V_c)) &= max_(dim V_c = n_c) quad min_(w in H^Ri) quad r_(Rb A) (w)\
  &= max_(dim H^Ri = n - (n_c + 1) + 1) quad min_(w in H^Ri) quad r_(Rb A) (w)\
  &= mu_(n_c + 1)
$
by the max-min theorem (Courant-Fisher), with ${mu_j , q_j}$ the eigenpairs of $Rb A$ in increasing order.

The optimal two-grids convergence is then $ ||E||_A = 1 - mu_(n_c + 1) $
which is achieved by choosing $V_c$ so that $range(P) = "vect"{q_1, ..., q_n_c}$, \ ie $H^Ri = "vect"{q_(n_c + 1), ..., q_n}$ since:
$
  min_(w in "vect"{q_(n_c + 1), ..., q_n})(||Rb A w||_Ri^2)/(||w||_Ri^2) = mu_(n_c + 1)
$

which can be obtained by setting $V_c = RR^(n_c)$ and $P = (q_1, q_2, ..., q_n)$
// un exemple pour rigoler ?

== Algebraic high and low frequencies

=== Smooth error

An error is smooth when it cannot be reduced by the smoother.

Given a smoother $R: quad V arrow V$ and $epsilon in ]0,1[$, we say that $v in V$ is $epsilon$-smooth if
$ ||v||_A^2 <= epsilon ||v||_Ri^2 $

The effect of the smoother on any error $v in V$ is given by $S v = (I - R A)v$\
we have:
$
  ||S v||_A^2 &= ((I - R A) v, (I - R A) v)_A \
  &= ((I - Rb A) v, v)_A \
  &= ||v||_A^2 - (Rb A v, v)_A
$

The norm reduction effect of the smoother is given by $(Rb A v, v)_A$, which will be small for $epsilon$-smooth errors:
$
  (Rb A v, v)_A <= epsilon ||v||_A arrow.double v "is" epsilon"-smooth"
$

$
  ||v||_A^2 = (Rb A v, Ri v) &<= (Rb A v, v)_A^(1/2)||v||_Ri\
  & <= sqrt(epsilon) ||v||_A||v||_Ri \
  arrow.double ||v||_A^2 <= epsilon ||v||_Ri
$

and equivalently because $||S v||_A^2 = ||v||_A^2 - (Rb A v, v)_A$, we have

$
  (Rb A v, v)_A <= epsilon ||v||_A arrow.l.r.double.long (||S v||_A^2)/(||v||_A^2) >= 1 - epsilon
$

which is another way to see that the (normalized) norm reduction effect of the smoother will be small for $epsilon$-smooth errors

== Numerical examples

For a general operator $E$:

$
  ||E||_A^2 = sup_(||x||_A =  1 )||E x||_A^2 = sup_(||x||_A =  1 )(E x, E x)_A
$
when E is not self adjoint with respect to A:
$
  ||E||_A^2 = sup_(||x||_A =  1 )(E^T A E x, x) = sup_(||x||_A =  1 )(A^(-1)E^T A E x, x)_A
$
with $M = A^(-1)E^T A E$ self adjoint with respect to A:
$
  ||E||_A^2 = lambda_max (A^(-1)E^T A E)
$
instead when E is self adjoint with respect to A, this simplifies to:
$
  ||E||_A^2 = sup_(||x||_A =  1 )(E x, E x)_A = sup_(||x||_A =  1 )(E^*E x, x)_A = lambda_max (E^*E)
$
which is why for $S = I - R A$ and $S^*S = overline(E) = I - Rb A$ we have:
$
  ||S||_A^2 = lambda_max (Sb) = 1 - lambda_min (Rb A)
$

in order to quantify the effect of the smoother, we will look at the eigenvalue of $Rb A$ instead of the more obvious $R A$

For a 1D Laplacian operator A = tridiag(-1, 2, -1) of size n, 
we have the eigenpairs for $ k in [| 1,n |]$:
$
  lk (A) = 4sin((k pi) / (2(n+1)))\
  q_k (A) = ( sin((j k pi) / (n+1)) )_j
$
=== Weighted Jacobi

$
  R = omega D^(-1) = omega /2 I 
$
for $omega = 2/3$,
$
  R = 1/3 I, quad S = I - 1/3 A \
  Rb = R^T + R - R^T A R = R(2 I - A R) = 1/3(2 I - 1/3 A)  
$
$A, R, S "and" Rb A$ share the same eigenvectors:
$
  lk (Rb A) = lk(Rb) lk(A) = 1/3lk(A)(2 - 1/3lk(A))
$
and we have:
$ A = U mat(
  lambda_1(A), , , ;
  , lambda_2(A), , ;
  , , dots.down, ;
  , , , lambda_n(A)
) U^T, quad
Tb = U mat(
  lambda_1(Tb), , , ;
  , lambda_2(Tb), , ;
  , , dots.down, ;
  , , , lambda_n(Tb)
) U^T $

to measure the action on $H^A$ we look at $ X = Tb (I - Pi_c)$

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
  , , , lambda_(n_c + 1)(Tb), , ;
  , , ,  dots.down, ,;
  , , , ,lambda_n (Tb)
) U^T quad X: V arrow V $

So 
$ X|_H = U_H mat(
   lambda_(n_c + 1)(Tb), , ;
   dots.down , ;
   quad quad quad quad lambda_n (Tb)
) U_H^T quad X|_H: H^A arrow H^A, quad U_H = U[ n_c + 1 : n] $

The two grid convergence is governed by
$
||E||_A^2 &= 1 - lmin(X)\ 
&= lambda_(n_c + 1)(Tb)\
&= 1 - 1/3lambda_(n_c + 1)(A)(2 - 1/3lambda_(n_c + 1)(A))
$

=== Gauss-Seidel and Additive Schwartz

Here $R$ is no longer symmetric, and we do not have a simple expression for $Rb "or" Tb$.

We can compute the eigenvalues of $S = I - R A$ and $X$ directly for small values of $n$.

// add figures later

The optimal coarse space is no longer composed of eigenvectors of A.\
To compute the orthogonal diagonalisation of $Tb$, we can note that the problem is equivalent to:
$
  Tb v = lambda v arrow.double.l.r & A Rb A v = lambda A v \
  arrow.double.l.r & M v = lambda A v \
$

$M$ is symmetric, we solve the generalized eigenproblem $M v = lambda A v$ for symmetric M and A, \ 
and we get an A-orthogonal basis:

$
  M v = lambda A v   arrow.double.l.r & (L^(-1) M L^(-T)) L^T v = lambda L^(-1) L L^T v quad, quad A = L L^T \
  arrow.double.l.r & tilde(M) w = lambda w \
$

Solving $tilde(M) w = lambda w$ yields an orthonormal basis $W = L^T V$,\
such that $W^T W = I arrow.double V^T L L^T V = I => V^T A V = I$

which we note is equivalent to getting an $Rb^(-1)$-orthonormal basis since:
$
  Tb v = lambda v & arrow.double.long Rb A v = lambda v \
  & arrow.double.long A v = lambda Rb^(-1) v
$
Solving $A v = lambda Rb^(-1) v$ yields a $Rb^(-1)$-orthonormal basis $V^T Rb^(-1) V = I$.

Furthermore, $Tb = Rb A$ is self-adjoint with respect to both the $A$-inner product and the $Rb^(-1)$-inner product:
$
  ( Tb x, y )_A &= x^T A (Rb A y) = x^T A Rb A y = ( x, Tb y )_A \
  ( Tb x, y )_(Rb^(-1)) &= x^T Rb^(-1) (Rb A y) = x^T A y = ( x, Tb y )_(Rb^(-1))
$
// both formulations generate the same optimal eigenspaces $V_(n_c) = "span"(v_1, ..., v_(n_c))$. ?

= Personal notes

== Matrix representation of X

$ X|_V = (I - Pi_c)Rb A : V arrow V$ is A-self adjoint\
Its matrix representation is $tilde(X)|_V in MM_(n times n)(RR)$

what is the matrix representation of $X : H arrow H$ ?

let $Q in RR^(n times n)$ be an A-orthonormal basis of V, $Q_H$ be an orthogonal basis of H 
$
 Q^T A Q = I\
 Q = mat(q_1 | dots | q_n) \
 Q_H = mat(q_(n_c+1) | dots | q_n) \
$

In operator form, 

$
 v in H subset V, quad v = sum_(i = n_c + 1)^(n)v_i bold(q_i) = Q_H tb(v) quad tb(v) in RR^(n-n_c)
$
$w = X v$ is well defined and $w in H$
$
  w = X v &= sum_(i = n_c + 1)^(n)(I - Pi_c)Rb A v_i bold(q_i) \
  & = Q_H tb(w) quad tb(w) in RR^(n-n_c)
$
$
  w_j &= ( bold(q_j) , X v)_A 
$
and we have:
$
  w = X v &arrow.double Q_H tb(w) = X Q_H tb(v) \
  &arrow.double tb(w) = Q_H^T A X Q_H tb(v)\
  &arrow.double tb(w) = tilde(X) tb(v)\
$
$ tilde(X) = Q_H^T A (I - Pi_c)Rb A Q_H in RR^((n-n_c) times (n-n_c))$ is the matrix representation of X (in the euclidean basis ? $Q_H$ ?)
noting that $tilde(X)$ is symmetric (X isn't)

In numerical application, if we have an euclidean orthonormal basis $U$ of eigenvectors of $A$ SPD:
$
  A = U D U^T, quad U^T U = I\
  A = H^2, quad H = U D^(1/2) U^T
$

we have
$
  U^T U = I &arrow.double U^T H^(-1) A H^(-1) U = I\
  &arrow.double (H^(-1) U)^T A (H^(-1) U) = I\
$
a valid $Q$ is $Q = H^(-1) U = U D^(-1/2)$
