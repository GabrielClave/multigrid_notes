#import "@preview/cetz:0.5.1": canvas, draw, decorations

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

= Subspace decomposition

#let fig-subspace-decomposition = figure(
  caption: none,
  canvas({
    import draw: *

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
    decorations.brace((xP - 0.25*margin, y0 - brace-margin), (x0, y0 - brace-margin), name: "b_low")
    content(
      (rel: (0, -margin), to: "b_low.center"),
      [#text("low frequency error")],
      anchor: "north"
    )

    // high freq brace
    decorations.brace((xend, y0 - brace-margin), (xP + 0.25*margin, y0 - brace-margin), name: "b_low")
    content(
      (rel: (0, -margin), to: "b_low.center"),
      [#text("High frequency error")],
      anchor: "north"
    )
  })
)

#fig-subspace-decomposition

= X

#let fig-X = figure(
  caption: none,
  canvas({
    import draw: *

    let circle_radius = 0.06
    let margin = 1
    let main-purple = rgb("#7570b3")
    let red-accent  = rgb("#e04040")

    // Define main points
    let orig = (0, 0)
    let eh = 5.5
    let ehs = 0.9
    let ep = 7.5
    let eps = 5
    let e-pt = (ep, eh)          // Initial error vector e
    let e-H = (0, eh)          // High-frequency component e_H
    let e-P = (ep, 0)            // Coarse-space component e_P
    let e-smoothed = (eps, ehs)   // Error after smoother: (I - T_bar) e
    let e-H-smoothed = (0, ehs) // Reduced high-frequency component (I - X) e_H
    let e-P-smoothed = (eps, 0) // Slightly reduced coarse component
    let midpoint = (0, (ehs+eh)/2)

    // Set coordinate system and axes
    line((-0.5, 0), (10, 0), stroke: 1pt, name: "x-axis")
    line((0, -0.5), (0, 7), stroke: 1pt, name: "y-axis")

    content((10.3, 0), $text(range)(P)$, anchor: "north-west")
    content((0, 7.3), $H^A$, anchor: "east")

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
    content(e-smoothed, padding: 0.15, $(I - Rb A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(midpoint, padding: 0.15, $(I - Pi_c)Rb A e_H\ = X e_H$, anchor: "east")
    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H^A component (reduction arrow)
    line((0, eh), (0, ehs), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1pt))
    content(midpoint, padding: 0.15, text(size: 14pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "west")

    // content(e-H-smoothed, padding: 0.15, $(I - Pi_c)(I - Rb A) e_H$, anchor: "east")
    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })
)

#fig-X

= Z

#let fig-Z = figure(
  canvas({
    import draw: *

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
)

#fig-Z

= unused 1


// #align(center)[
//   #canvas({
//     import draw: *

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