#import "@preview/cetz:0.5.1": canvas, draw, decorations
#import "@preview/touying:0.8.0": *
#import themes.metropolis: *

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
#let Mi = $overline(M)^(-1)$
#let Mb = $overline(M)$
#let Tb = $overline(T)$
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let lk = $lambda_k$
#let lmin = $lambda_min$

== Smoother spectrum

// #let fig-spectrum-split(step: 2) = canvas({
//   import draw: *

//   let trajectory-red  = rgb("#e76f51")
//   let trajectory-teal = rgb("#2a9d8f")

//   // 1. Render image inside draw.image with named anchor
//   image("../../pictures/plot_J_spectrum.png", width: 80%)

//   // 2. Coordinates tuned to the plot area relative to plot.south-west
//   let x-split = 0.52
//   let split-top    = (rel: (x-split * 14, 6.2), to: "plot.south-west")
//   let split-bottom = (rel: (x-split * 14, 1.1), to: "plot.south-west")

//   // 3. Conditional drawing within CeTZ canvas
//   if step >= 2 {
//     // Shaded region for H (effective smoothing domain)
//     rect(
//       split-top,
//       (rel: (5.8, -5.1), to: split-top),
//       fill: rgb("#e76f5115"),
//       stroke: none
//     )

//     // Vertical dashed demarcation line
//     line(
//       split-top,
//       split-bottom,
//       stroke: (dash: "dashed", paint: trajectory-red, thickness: 1.5pt)
//     )

//     // Label content definitions (wrapped properly to prevent join errors)
//     let label-L = [
//       #text(fill: gray.darken(30%), weight: "bold", size: 0.9em)[$L$ (Low Frequencies)] \
//       #text(fill: gray.darken(30%), size: 0.8em)[$lambda(S) approx 1 arrow.r$ Stagnation]
//     ]

//     let label-H = [
//       #text(fill: trajectory-red, weight: "bold", size: 0.9em)[$H$ (High Frequencies)] \
//       #text(fill: trajectory-red, size: 0.8em)[$lambda(S) thin \ll thin 1 arrow.r$ Damped]
//     ]

//     // Subspace L
//     content(
//       (rel: (-2.5, -0.4), to: split-top),
//       label-L,
//       anchor: "south"
//     )

//     // Subspace H
//     content(
//       (rel: (2.5, -0.4), to: split-top),
//       label-H,
//       anchor: "south"
//     )
//   }
// })

== Smoother effect

#let fig-smoother(step: 2) = canvas({

  import draw: *

  let circle_radius = 0.06
  let trajectory-purple = rgb("#7570b3")
  let trajectory-teal   = rgb("#2a9d8f")
  let trajectory-red    = rgb("#e76f51")

  let xmax = 15
  let ymax = 11

  let eh = 10
  let ep = 14
  let eh-smoothed = 1.2
  let ep-smoothed = 10.2

  let e-pt = (ep, eh)
  let e-H = (0.0, eh)
  let e-P = (ep, 0.0)
  let e-smoothed = (ep-smoothed, eh-smoothed)
  let e-P-smoothed = (ep-smoothed, 0.0)
  let e-H-smoothed = (0.0, eh-smoothed)

  let padding = - 0.3

  // Axes & Initial state (Step 1+)
  line((-0.5, 0.0), (xmax, 0.0), stroke: 1pt)
  line((0.0, -0.5), (0.0, ymax), stroke: 1pt)
  content((xmax - padding, 0.0), $L$, anchor: "west")
  content((0.0, ymax - padding), $H$, anchor: "south")

  line(e-pt, e-H, stroke: (dash: "dashed", paint: gray))
  line(e-pt, e-P, stroke: (dash: "dashed", paint: gray))
  content(e-H, padding: 0.15, $e_H$, anchor: "east")
  content(e-P, padding: 0.15, $e_L$, anchor: "north")
  circle(e-pt, radius: circle_radius, fill: black, stroke: none)
  content(e-pt, $e$, padding: 0.15, anchor: "south-west")

  // Smoother action (Step 2+)
  if step == 2 {
    line(e-smoothed, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
  }
  if step >= 2 {
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))
    line(
      e-pt, e-smoothed,
      mark: (end: "triangle", fill: trajectory-teal),
      stroke: (paint: trajectory-teal, thickness: 2pt)
    )

    // Projections to axes
    line(e-smoothed, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))

    content(e-H-smoothed, padding: 0.15, text()[$(S e)_H$], anchor: "east")
    content(e-P-smoothed, padding: 0.15, text()[$(S e)_L$], anchor: "north")

    circle(e-smoothed, radius: circle_radius, fill: black, stroke: none)
    content(e-smoothed, padding: 3, text(fill: trajectory-teal)[smoother \ effect], anchor: "south")
    content(e-smoothed, padding: 0.2, text()[$S e$], anchor: "north-west")
  }

  // Component decomposition & reduction highlighting (Step 3)
  if step >= 3 {


    // Damping annotations along the axes
    // Strong reduction in H (large drop)
    line(
      (0, eh), (0, eh-smoothed),
      mark: ( end: "straight", fill: trajectory-red),
      stroke: (paint: trajectory-red, thickness: 1.5pt)
    )
    content(
      (padding, (eh + eh-smoothed) / 2),
      text(fill: trajectory-red, size: 0.85em)[#align(center)[strong \ damping]],
      anchor: "east"
    )

    // Weak reduction in L (small drop)
    line(
      (ep, 0), (ep-smoothed, 0),
      mark: (end: "straight", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 1.5pt)
    )
    content(
      ((ep + ep-smoothed) / 2, padding),
      text(fill: trajectory-purple, size: 0.85em)[#align(center)[small \ reduction]],
      anchor: "north"
    )
  }

})

#fig-smoother(step: 1)
#fig-smoother(step: 2)
#fig-smoother(step: 3)

== coarse space correction (projection)

#let fig-coarse-correction(step: 2) = canvas({

  import draw: *

  let circle_radius = 0.06
  let trajectory-purple = rgb("#7570b3")
  let trajectory-teal   = rgb("#2a9d8f")
  let trajectory-orange = rgb("#e76f51")

  let xmax = 15
  let ymax = 11

  // Start from the exact same initial error point e
  let eh = 8
  let ep = 11
  let e-pt = (ep, eh)

  // Coarse space range axis along x, orthogonal complement along y
  let proj-y-pt = (0.0, eh)
  let proj-x-pt = (ep, 0.0)

  let padding = -0.3

  // Axes: x-axis represents range(P)
  line((-0.5, 0.0), (xmax, 0.0), stroke: 1pt)
  line((0.0, -0.5), (0.0, ymax), stroke: 1pt)
  content((xmax - padding, 0.0), $range(P)$, anchor: "west")

    // Dashed line dropping down to range(P)
    line(
      e-pt, proj-x-pt,
      stroke: (dash: "dashed", paint: gray)
    )
    line(
      e-pt, proj-y-pt,
      stroke: (dash: "dashed", paint: gray)
    )
    if step == 1 {
      content(
        proj-x-pt,
        padding: 0.15,
        $e_L$,
        anchor: "north"
      )
      content(
        proj-y-pt,
        padding: 0.15,
        $e_H$,
        anchor: "east"
      )
    }
  // Initial error point e
  circle(e-pt, radius: circle_radius, fill: black, stroke: none)
  content(e-pt, $e$, padding: 0.15, anchor: "south-west")

  // Step 1: Projection vector & dashed alignment to range(P)
  if step > 1 {

    // Purple vector showing coarse grid correction to (I - Pi) e
    line(
      e-pt, proj-y-pt,
      mark: (end: "triangle", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 2pt)
    )

    // Might angle indicator at the y-axis (0, eh)
    let ra-size = 0.4
    line(
      (0.0, eh - ra-size),
      (ra-size, eh - ra-size),
      (ra-size, eh),
      stroke: (paint: gray, thickness: 0.8pt)
    )

    // Projected point on the y-axis: remaining error (I - Pi) e
    circle(proj-y-pt, radius: circle_radius, fill: trajectory-purple, stroke: none)
    content(
      proj-y-pt,
      padding: 0.15,
      text()[$(I - Pi_c) e$],
      anchor: "east"
    )
    content(
      proj-x-pt,
      padding: 0.15,
      text()[$Pi_c e$],
      anchor: "north"
    )
  }

})

#fig-coarse-correction(step: 1)
#fig-coarse-correction(step: 2)

== two-grid iteration

#let fig-twogrid(step: 4) = canvas({

  import draw: *

  let circle_radius = 0.06
  let trajectory-purple = rgb("#7570b3")
  let trajectory-teal   = rgb("#2a9d8f")
  let trajectory-red    = rgb("#e76f51")

  let xmax = 15
  let ymax = 11

  // Tilt angle for the L/H basis relative to range(P) and its orthogonal complement
  let theta = 8deg

  // Helper projection function onto tilted axes (L rotated by theta)
  // Maps coordinates (u_L, u_H) in the tilted basis to canvas (x, y)
  let to-canvas(u_L, u_H) = (
    u_L * calc.cos(theta) - u_H * calc.sin(theta),
    u_L * calc.sin(theta) + u_H * calc.cos(theta)
  )

  let padding = -0.3

  // -------------------------------------------------------------
  // Axes Layout
  // -------------------------------------------------------------

  // 2. Dashed Tilted Axes: L and H (extended past 0 by 0.5 units)
  let axis-ext = 1
  let L-start = (-axis-ext * calc.cos(theta), -axis-ext * calc.sin(theta))
  let L-end   = (xmax * calc.cos(theta), xmax * calc.sin(theta))

  let H-start = (axis-ext * calc.sin(theta), -axis-ext * calc.cos(theta))
  let H-end   = (-ymax * calc.sin(theta), ymax * calc.cos(theta))
  
  line(L-start, L-end, stroke: (dash: "dashed", paint: gray, thickness: 0.8pt))
  line(H-start, H-end, stroke: (dash: "dashed", paint: gray, thickness: 0.8pt))
  
  content(L-end, padding: 0.15, $L$, anchor: "south-west")
  content(H-end, padding: 0.15, $H$, anchor: "south-east")

  // 1. Solid Axes: range(P) along x-axis, and its orthogonal complement along y-axis
  line((-axis-ext, 0.0), (xmax, 0.0), stroke: 1pt)
  line((0.0, -axis-ext), (0.0, ymax), stroke: 1pt)
  content((xmax - padding, 0.0), $range(P)$, anchor: "west")

  // -------------------------------------------------------------
  // Trajectory State Coordinates
  // -------------------------------------------------------------
  // Step 1: Initial Error e (large L component, large H component)
  let e-pt = (12.0, 8.5)
  
  // Step 2: Coarse Correction (I - Pi_c) e -> Orthogonal drop onto y-axis (0, y_e)
  let e1-pt = (0.0, e-pt.at(1))

  // Step 3: Smoother S e1 -> Damps H component heavily, small change in L in tilted basis
  // Express e1 in (L, H) coordinates, shrink H component, project back to canvas
  // e1 in canvas is (0, 8.5). In tilted basis: u_L = 8.5*sin(theta), u_H = 8.5*cos(theta)
  // Damping u_H heavily gives e2-pt
  let e2-pt = (3, 1.8)

  // Step 4: Second Coarse Correction (I - Pi_c) e2 -> Orthogonal drop onto y-axis
  let e3-pt = (0.0, e2-pt.at(1))

  // -------------------------------------------------------------
  // Step 1: Initial Error
  // -------------------------------------------------------------
  circle(e-pt, radius: circle_radius, fill: black, stroke: none)
  content(e-pt, $e^(0)$, padding: 0.15, anchor: "south-west")

  // -------------------------------------------------------------
  // Step 2: Coarse Grid Correction
  // -------------------------------------------------------------
  if step >= 2 {
    // Drop orthogonal to range(P) onto y-axis
    line(
      e-pt, e1-pt,
      mark: (end: "triangle", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 2pt)
    )

    circle(e1-pt, radius: circle_radius, fill: trajectory-purple, stroke: none)
    content(
      e1-pt,
      padding: 0.3,
      text()[$e^(1) = (I - Pi_c) e^(0)$],
      anchor: "east"
    )
  }

  // -------------------------------------------------------------
  // Step 3: Smoother Sweep
  // -------------------------------------------------------------
  if step >= 3 {
    // Smoother action: damps H along the tilted space
    line(
      e1-pt, e2-pt,
      mark: (end: "triangle", fill: trajectory-teal),
      stroke: (paint: trajectory-teal, thickness: 2pt)
    )

    circle(e2-pt, radius: circle_radius, fill: trajectory-teal, stroke: none)
    content(
      e2-pt,
      padding: 0.15,
      text()[$e^(2) = S e^(1)$],
      anchor: "west"
    )
  }

  // -------------------------------------------------------------
  // Step 4: Second Coarse Grid Correction
  // -------------------------------------------------------------
  if step >= 4 {
    // Second orthogonal drop onto y-axis
    line(
      e2-pt, e3-pt,
      mark: (end: "triangle", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 2pt)
    )

    circle(e3-pt, radius: circle_radius, fill: trajectory-purple, stroke: none)
    content(
      e3-pt,
      padding: 0.4,
      text()[$e^(3) = (I - Pi_c) e^(2)$],
      anchor: "east"
    )
  }
})

#fig-twogrid(step: 1)
#fig-twogrid(step: 2)
#fig-twogrid(step: 3)
#fig-twogrid(step: 4)


== Smoother convergence (informal)

// add 2 steps
// 1: two-grid cycle
// 2: make the "red arrow" smoother effect appear

#let fig-smoother-convergence(step: 2) = canvas({

  import draw: *

  let circle_radius = 0.06
  let main-purple = rgb("#7570b3")
  let red-accent  = rgb("#e04040")

  let xmax = 15
  let ymax = 11
  let padding = -0.3

  // 16:9 Scaled Coordinates
  let eh = 8.5
  let ehs = 1.8
  let ep = 12.0
  let eps = 3.5

  let e-pt = (ep, eh)            // Initial error e
  let e-H = (0.0, eh)            // e_H = (I - \Pi_c) e
  let e-P = (ep, 0.0)            // e_P = \Pi_c e

  let e-H-smoothed = (0.0, ehs)  // Remaining error on H
  let e-final = (eps, ehs)       // Final error state
  let e-P-final = (eps, 0.0)

  let midpoint-H = (0.0, (ehs + eh) / 2)

  // Axes
  line((-0.5, 0.0), (xmax, 0.0), stroke: 1pt)
  line((0.0, -0.5), (0.0, ymax), stroke: 1pt)

  content((xmax - padding, 0.0), $range(P)$, anchor: "west")
  content((0.0, ymax - padding), $H^A$, anchor: "south")

  // --- INITIAL STATE e & DECOMPOSITION ---
  line(e-pt, e-H, stroke: (dash: "dashed", paint: gray))
  line(e-pt, e-P, stroke: (dash: "dashed", paint: gray))

  circle(e-pt, radius: circle_radius, fill: black, stroke: none)
  content(e-pt, $e$, padding: 0.15, anchor: "south-west")

  content(e-H, padding: 0.15, $e_H = (I - Pi_c) e$, anchor: "east")

  circle(e-P, radius: circle_radius, fill: black, stroke: none)
  content(e-P, padding: 0.15, $e_P = Pi_c e$, anchor: "north")

  // --- STEP 1: Two-grid cycle trajectory ---
  if step >= 1 {
    // 1st Purple Arrow: Coarse projection e -> e_H
    line(
      e-pt, e-H,
      mark: (fill: main-purple),
      stroke: (paint: main-purple, thickness: 1.5pt)
    )

    // 2nd Purple Arrow: Smoother action from e_H to final error
    line(
      e-H, e-final,
      mark: (end: "triangle", fill: main-purple),
      stroke: (paint: main-purple, thickness: 1.5pt)
    )

    // Dashed projections for final error
    line(e-final, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
    line(e-final, e-P-final, stroke: (dash: "dashed", paint: gray))

    circle(e-final, radius: circle_radius, fill: black, stroke: none)
    content(e-final, padding: 0.15, $(I - M A)e_H$, anchor: "west")
  }

  if step == 1{
    content(e-H-smoothed, padding: 0.15, $(I - Pi_c)(I - M A) e_H$, anchor: "east")
  }

  // --- STEP 2: Med smoother effect along H ---
  if step >= 2 {
    line(
      e-H, e-H-smoothed,
      mark: (end: "triangle", fill: red-accent),
      stroke: (paint: red-accent, thickness: 1.5pt)
    )

    content(e-H-smoothed, padding: 0.15, $(I - X) e_H$, anchor: "east")

    content(
      midpoint-H,
      padding: 0.2,
      text(size: 18pt, fill: red-accent)[
        Smoother effect on $H^A:$ \
        $ (I - Pi_c) M A e_H $
      ],
      anchor: "east"
    )
  }
})

#fig-smoother-convergence(step: 1)
#fig-smoother-convergence(step: 2)

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

    // Complementary space H
    content(((xP + xend) / 2, (y0 + yend) / 2), $H = "range"(P)^(perp_A)$)

    // Outer dashed box for range(P)
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
    content((0, 8 + 0.5), $H$, anchor: "south")

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
    content(e-smoothed, padding: 0.15, $(I - M A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(midpoint, padding: 0.15, $(I - Pi_c)M A e\ = "important part"$, anchor: "east")
    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H component (reduction arrow)
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
    content((0, 8 + 0.5), $H$, anchor: "south")

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

    // Vector (I - MA)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - M A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
        content(midpoint, padding: 0.15, $(I - Pi_c)M A e_H$, anchor: "east")

    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H component (reduction arrow)
    line((0, eh), (0, 0), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1.5pt))
    // content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    content(e-H-smoothed, padding: 0.5, $(I - Pi_c)(I - M A) e_H$, anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })

  #fig-optimal-smoother

== V-cycle

#let fig-v-cycle = canvas({
  import draw: *

  let dx = 1.9
  let dy = 3.5
  let r-node = 0.1
  let p_padding = 0.7

  let c-primary = rgb("#7570b3")
  let c-coarse = rgb("#ad1b1bd7")
  let c-guide = luma(210)

  let spaces = ($V$, $V_1$, $V_2$, $V_3$)
  let eqs = ($A e = r$, $A_1 e_1 = r_1$, $A_2 e_2 = r_2$, $e_3 = A_3^(-1) r_3$)
  let rests = ($P_1^T$, $P_2^T$, $P_3^T$)
  let prols = ($P_1$, $P_2$, $P_3$)

  // Subtle horizontal level guides and space labels on the left
  for k in std.range(4) {
      let y = (3 - k) * dy
      let y-eq = if k == 3 { y - 1 } else { y + 0.3 }
      line((-1.4, y), (6 * dx + 1.4, y), stroke: (paint: c-guide, dash: "dashed", thickness: 0.6pt))
      content((-1.7, y), anchor: "east", spaces.at(k))
      content((3 * dx, y-eq), anchor: "south", eqs.at(k))
  }

  // Descent nodes (pre-smoothing levels: V, V_1, V_2)
  for k in std.range(3) {
    let x = k * dx
    let y = (3 - k) * dy
    let name = "down-" + str(k)
    circle((x, y), radius: r-node, fill: white, stroke: 1.2pt + c-primary, name: name)
  }

  // Coarsest grid node (V_3)
  let x-bot = 3 * dx
  circle((x-bot, 0), radius: r-node, fill: c-coarse, stroke: 1.2pt + c-coarse, name: "bottom")

  // Ascent nodes (post-smoothing levels: V_2, V_1, V)
  for k in std.range(3) {
    let lvl = 2 - k
    let x = (4 + k) * dx
    let y = (k + 1) * dy
    let name = "up-" + str(lvl)
    circle((x, y), radius: r-node, fill: white, stroke: 1.2pt + c-primary, name: name)
  }

  // Restriction arrows (downward)
  let down-names = ("down-0", "down-1", "down-2", "bottom")
  for k in std.range(3) {
    let edge = "rest-" + str(k)
    line(
      down-names.at(k),
      down-names.at(k + 1),
      stroke: 1.1pt + black,
      mark: (end: ">", fill: black, scale: 0.8),
      name: edge,
    )
    content(edge + ".mid", anchor: "east", padding: p_padding, rests.at(k))
  }

  // Prolongation arrows (upward)
  let up-names = ("bottom", "up-2", "up-1", "up-0")
  for k in std.range(3) {
    let prol-idx = 2 - k
    let edge = "prol-" + str(prol-idx)
    line(
      up-names.at(k),
      up-names.at(k + 1),
      stroke: 1.1pt + black,
      mark: (end: ">", fill: black, scale: 0.8),
      name: edge,
    )
    content(edge + ".mid", anchor: "west", padding: p_padding, prols.at(prol-idx))
  }
})

#fig-v-cycle 