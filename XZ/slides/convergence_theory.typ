// Import the Touying package and the Metropolis theme
#import "@preview/touying:0.5.3": *
#import themes.metropolis: *
#import "@preview/cetz:0.5.1"

// #set text(font: "Fira Sans")

// Initialize the theme
#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Multigrid Convergence],
    author: [Gabriel Clave],
    date: datetime.today().display(),
  ),
)

#let Ri = $overline(R)^(-1)$
#let range = math.op("range")
#let nul = math.op("null")
#let im = math.op("Im")
#let ker = math.op("ker")

// Generate the title slide
#title-slide()

== Multigrid recap

A multigrid solver relies on three core ingredients:

- A *Smoother $S$* eliminates error in a subspace $H$.
- A *Coarse Grid* & Interpolation ($P$)
- A *Coarse Solver* (can recursively be another multigrid cycle).

*Fundamental Principle:*\
The central objective when designing a multigrid method is to construct $P$ and $S$ such that $range(P)$ and $H$ are *complementary*:
$ V approx range(P) xor_A H $

Ensuring that error components missed by the smoother are effectively captured and eliminated by the coarse grid correction.

#pagebreak()

#align(center)[
  #image("../../pictures/MG Convergence(2).png", width: 80%)
  ]

= Classic Convergence Theory

== Smooth error

An error is smooth when it cannot be reduced effectively by the smoother
#v(1cm)
$ S e &approx e \
  (I - M^(-1)N)e &approx 0 \
  (M - N)e &approx 0 \
  A e = r &approx 0 $

== Oscillatory error

An error is oscillatory when has a "high energy":\
it has directions associated with large eigenvalues of A.
#v(1cm)

$ ||e||_A^2 = ( A e, e ) $

$ ( A e, e ) = sum_k lambda_k c_k^2 $

Huge if components are along the eigenvectors associated with large eigenvalues.

== Smoothing property
*The Smoothing Property*: The smoother $S$ reduces high-frequency errors.
$ ||S e||_A^2 <= ||e||_A^2 - alpha ||A e||_(D^(-1))^2 $

- $S$ is a contraction in the $A$-norm, and the relaxation converge.
- The negative term quantifies the drop: 
  - If the error is smooth, $A e approx 0$, and progress stalls.
  - If the error is oscillatory, $||r||_(D^(-1))^2 = ||A e||_(D^(-1))^2 = sum_k lambda_k^2 c_k^2$ is huge and progress is fast.

== The approximation property
*The Approximation Property*: The coarse grid can approximate smooth errors.
$ min_(u_c) ||e - P u_c||_D^2 <= beta ||e||_A^2 $

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
    content((0, 4.3), $H^A$, anchor: "south")

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
    content(v-pt, padding: 0.15, $e$, anchor: "south-west")

    // Coarse component Q_c v on range(P)
    content(Qc-v, padding: 0.15, $P u_c$, anchor: "north")

    // Complement component (I - Q_c)v on y-axis
    content(IQc-v, padding: 0.15, $e - P u_c$, anchor: "east")

    // --- DISTANCE ANNOTATION (Double-edged arrow along y-axis) ---
    line(
      (vp + 0.15, 0.1), (vp + 0.15, 3.1),
      mark: (start: "triangle", end: "triangle", fill: red-accent),
      stroke: (paint: red-accent, thickness: 0.6pt),
      name: "distance"
    )

    content(
      (rel: (2.5, 0), to: "distance"),
      text(size: 14pt)[
        $d_D (e, text(range)(P)) \ = min_(u_c)||e - P u_c||_D$
      ],
      anchor: "center"
    )
  })
]

== The approximation property
*The Approximation Property*: The coarse grid can approximate smooth errors.
$ min_(u_c) ||e - P u_c||_D^2 <= beta ||e||_A^2 $
- If the error is oscillatory, $||e||_A$ is huge, and the inequality is easily satisfied, regardless of the quality of the interpolation.
- If $e$ is a smooth error, the interpolation needs to accurately represent it: smooth errors must be near $op("range")(P)$.

== Convergence theorem
*Convergence theorem* \
If the smoothing and the approximation property is true, and if $alpha < beta$,\ the iteration converges, with:
$ ||S T G||_A <= sqrt(1 - alpha/beta) < 1 $

#pagebreak()

== The problemsW_P\^\(perp_A\)

- developed in a geometric settings
- In practice $alpha$ and $beta$ are out of reach
- sufficient but not necessary conditions


= XZ Theory

== The promise 

- "In this paper, we try to develop a *unified framework* and theory that can be used to derive and analyze different algebraic multigrid methods in a coherent manner."

- "Our theory applies to *most existing multigrid methods*, including the standard geometric multigrid method, the classic AMG, energy-minimization AMG, unsmoothed and smoothed aggregation AMG, and spectral AMGe."
// and they develop these method in their framework

// - they allegedly developed a unified theory ("subspace correction methods") applicable for domain decomposition and multigrid

== Some notations (smoothers)
The first component of a multigrid solver is a smoother.\
It is a linear operator $R : V arrow.r.long V$ that can be used to solve the system through the  iteration:
$
  v^((k+1)) = v^((k)) + R(f - A v^((k)))
$
with associated error propagation operator
$
  E = I - R A
$
If $R$ is not symmetric, it can be symmetrized through:

$
  overline(R) = R + R^T - R^T A R
$

with
$
  overline(E) = I - overline(R) A = (I - R^T A)(I - R A)\  
$

== Some notations (subspace decomposition)

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

    let brace-margin = 2.5*margin

    set-style(
      content: (padding: 0.1),
      text: (size: 20pt)
    )

    // Entire space V
    rect((x0, y0), (xend, yend), stroke: 1pt, name: "V_rect")

    // Complementary space H^A
    content(((xP + xend) / 2, (y0 + yend) / 2), $H^A = "range"(P)^(perp_A)$)

    // Outer dashed box for Range(P)
    rect(
      (x0, y0),
      (xP, yend),
      fill: light-purple
    )

    content(((xP / 2), (y0 + yend) / 2), $"range"(P)$)

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

== Some notations (projections)

We define the following $A$-orthogonal projections:

$
  Pi_c : V &arrow.r.long range(P)\

$

#v(1cm)

We define the following $Ri$-orthogonal projections:

$
 Q_c: V &arrow.r.long range(P)\
$

#pagebreak()

== Two-Grids Operator

A two-grid correction operator is an operator $B : V' -> V$ that acts on $f in V'$ with two steps:
- Apply the subspace correction:
  $ w = P A_c^(-1)P' f $
- Apply the smoother:
  $ B f = w + R(f - A w) $

The associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

#pagebreak()

== The Convergence Theorem

We have:
$ ||E||_A^2 = 1 - lambda_"min"(X) $
where $ X = (I - Pi_c) overline(T) quad : quad H^A -> H^A $

$X$ represents the effect of the (symmetrized) smoother on the non-smooth errors.

== The Operator X

$X$ represents the effect of the (symmetrized) smoother on the non-smooth errors.

#v(2cm)
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
    let e-smoothed = (3, 0.9)   // Error after smoother: (I - T_bar) e
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

    // --- STATE AFTER SMOOTHER (I - T_bar)e ---
    // Dashed projection lines for smoothed error
    line(e-smoothed, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))

    // Vector (I - T_bar)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - overline(T))e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(e-H-smoothed, padding: 0.15, $(I - X)e_H$, anchor: "east")
    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H^A component (reduction arrow)
    line((0, 3.5), (0, 0.9), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1pt))
    content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })
]

#pagebreak()

the effect of the two-grids cycle on $e in H$ is:
$ overline(E)|_H &= (I - overline(T))( I - Pi_c ) \
                 &= I - overline(T) quad (I - Pi_c = I "on" H) $

The components along $H$ are given by:
$ (I - Pi_c)overline(E)|_H &= I - (I - Pi_c)overline(T) \
                            &= I - X $

#pagebreak()

$ ||E||_A^2 = 1 - lambda_"min"(X) $

*Interpretation of $X$:* \
- If $X approx I$ on $H$, the smoother is effective and the error is almost entirely annihilated: \
  $lambda_"min"(X) approx 1$ and $||E||_A^2 = 1 - lambda_"min"(X) approx 0$

- If $X$ has a small eigenvalue, the associated direction in $H$ will escape both the effect of the smoother and the coarse space correction: convergence will stagnate.

==  The operator Z

The inverse of $X$ on $H^A$ can be explicitly written as:
$ Z = overline(T)^(-1) (I - Q_c) quad : quad H^A -> H^A $

So we have $ lambda_"min"(X) = 1 / (lambda_"max" (Z)) $

#pagebreak()

*Interpretation of $Z$:* \
After applying the smoother, the remaining smooth error should be as close as possible to $op("range")(P)$. \
The distance is $d_(Ri) (e, op("range")(P)) = min_(e_P in op("range")(P)) ||e - e_P||_(Ri)^2 = ||(I - Q_c)v||^2_(Ri)$

$Z$ allows us to express this distance in the $A$-geometry:
$ (Z v, v)_A =  ||(I - Q_c)v||_(Ri)^2 $

== Z

And so in the worst-case scenario we have:
$ lambda_"max" (Z) &= max_(v in H^A) r_Z (v) \
                  &= max_(v in H^A) (Z v, v)_A / (v,v)_A \
                  &= max_(v in H^A) (||(I - Q_c)v||_(Ri)^2) / (||v||_A^2) \
                  &= K(V_c) $

#pagebreak()

$K(V_c)$ is (more or less) the minimum K that verifies the approximation property:
$ ||(I - Q_c)v||_(Ri)^2 <= K||v||_A^2 $


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
    content((0, 4.3), $H^Ri$, anchor: "south")

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

#pagebreak()

$ ||E||_A^2 = 1 - 1 / K(V_c) $

This can be understood as follows: \
the effectiveness of a two-grids method depends on how well the coarse space can approximate the smooth error (quantified by $Z$) that escaped the smoother (quantified by $X$).

== Optimal coarse space

Given the eigenpairs of $overline(R)A$ $(q_i , mu_i)$

The optimal two-grids convergence is when $ ||E||_A = 1 - mu_(n_c + 1) $
This is achieved by choosing $V_c$ so that $op("range")(P) = op("span"){q_1, ..., q_(n_c)}$,\
 i.e., $H^(Ri) = op("span"){q_(n_c + 1), ..., q_n}$, since:

$ min_(w in op("span"){q_(n_c + 1), ..., q_n}) (||overline(R)A w||_(Ri)^2) / (||w||_(Ri)^2) = mu_(n_c + 1) $

This can be obtained by setting $V_c = RR^(n_c)$ and $P = (q_1, q_2, ..., q_(n_c))$.