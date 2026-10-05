#import "@preview/cetz:0.5.1": canvas, draw, decorations

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

== AMG

#let fig-amg(step: 3) = canvas({
// #let fig-amg = canvas({
  import draw: *

  let circle_radius = 0.06
  let trajectory-purple = rgb("#7570b3")
  let trajectory-teal   = rgb("#2a9d8f")

  let eh = 5.5
  let ep = 7.5
  let eh-smoothed = 1.2
  let ep-smoothed = 5.2

  let e-pt = (ep, eh)
  let e-H = (0.0, eh)
  let e-P = (ep, 0.0)
  let e-smoothed = (ep-smoothed, eh-smoothed)
  let e-P-smoothed = (ep-smoothed, 0.0)
  let e-final = (0.0, eh-smoothed)

  // Axes & Initial state (Step 1+)
  line((-0.5, 0.0), (9.0, 0.0), stroke: 1pt)
  line((0.0, -0.5), (0.0, 7.0), stroke: 1pt)
  content((9.4, 0.0), $text(range)(P)$, anchor: "west")
  content((0.0, 7.4), $H^A$, anchor: "south")

  line(e-pt, e-H, stroke: (dash: "dashed", paint: gray))
  line(e-pt, e-P, stroke: (dash: "dashed", paint: gray))
  content(e-H, padding: 0.15, $e_H$, anchor: "east")
  content(e-P, padding: 0.15, $e_P$, anchor: "north")
  circle(e-pt, radius: circle_radius, fill: black, stroke: none)
  content(e-pt, $e$, padding: 0.15, anchor: "south-west")

  // Smoother action (Step 2+)
  if step == 2 {
    line(e-smoothed, (0,eh-smoothed), stroke: (dash: "dashed", paint: gray))
  }
  if step >= 2 {
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))
    line(
      e-pt, e-smoothed,
      mark: (end: "triangle", fill: trajectory-teal),
      stroke: (paint: trajectory-teal, thickness: 2pt)
    )
    circle(e-smoothed, radius: circle_radius, fill: black, stroke: none)
    content(e-smoothed, padding: 3, text(fill: trajectory-teal)[smoother], anchor: "south")
  }

  // Coarse grid correction (Step 3+)
  if step >= 3 {
    line(
      e-smoothed, e-final,
      mark: (end: "triangle", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 2pt)
    )
    circle(e-final, radius: circle_radius, fill: black, stroke: none)
    content(e-final, padding: 0.15, text(fill: trajectory-purple)[Coarse correction], anchor: "east")
  }
})

// #fig-amg()
#fig-amg(step: 1)
#fig-amg(step: 2)
#fig-amg(step: 3)

== subspace decomposition SPD

#let fig-subspace-decomposition-SPD = canvas({

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

== smoother error

#let fig-smoother-error = canvas({
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
    line((0, -0.5), (0, 8), stroke: 1pt, name: "y-axis")

    content((10 + 0.5, 0), $text(range)(P)$, anchor: "west")
    content((0, 8 + 0.5), $H^A$, anchor: "south")

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
    content(e-smoothed, padding: 0.15, $(I - R A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(midpoint, padding: 0.15, $(I - Pi_c)R A e_H\ = "important part"$, anchor: "east")
    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H^A component (reduction arrow)
    line((0, eh), (0, 0.9), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1pt))
    // content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })

  #fig-smoother-error

  == optimal snoother

#let fig-optimal-smoother = canvas({
    import draw: *

    let circle_radius = 0.06
    let margin = 1
    let main-purple = rgb("#7570b3")
    let red-accent  = rgb("#e04040")

    // Define main points
    let orig = (0, 0)
    let eh = 5.5
    let ehs = 0
    let ep = 7.5
    let eps = 7.5
    let e-pt = (ep, eh)          // Initial error vector e
    let e-H = (0, eh)          // High-frequency component e_H
    let e-P = (ep, 0)            // Coarse-space component e_P
    let e-smoothed = (eps, ehs)   // Error after smoother: (I - T_bar) e
    let e-H-smoothed = (0, ehs) // Reduced high-frequency component (I - X) e_H
    let e-P-smoothed = (eps, 0) // Slightly reduced coarse component
    let midpoint = (0, (ehs+eh)/2)

    // Set coordinate system and axes
    line((-0.5, 0), (10, 0), stroke: 1pt, name: "x-axis")
    line((0, -0.5), (0, 8), stroke: 1pt, name: "y-axis")

    content((10 + 0.5, 0), $text(range)(P)$, anchor: "west")
    content((0, 8 + 0.5), $H^A$, anchor: "south")

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

    // Vector (I - RA)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - R A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
        content(midpoint, padding: 0.15, $(I - Pi_c)R A e_H$, anchor: "east")

    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H^A component (reduction arrow)
    line((0, eh), (0, 0), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1.5pt))
    // content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    content(e-H-smoothed, padding: 0.5, $(I - Pi_c)(I - R A) e_H$, anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })