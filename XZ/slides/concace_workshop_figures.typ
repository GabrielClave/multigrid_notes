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
#let Bi = $overline(B)^(-1)$
#let Bb = $overline(B)$
#let Tb = $overline(T)$
#let Sb = $overline(S)$
#let lm = $lambda_max$
#let lk = $lambda_k$
#let lmin = $lambda_min$

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
  let ep-smoothed = 8.2

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
    content(e-smoothed, padding: 2, text(fill: trajectory-teal)[smoother \ effect], anchor: "south")
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
      text(fill: trajectory-red, size: 0.85em)[strong \ damping],
      anchor: "east"
    )

    // Weak reduction in L (small drop)
    line(
      (ep, 0), (ep-smoothed, 0),
      mark: (end: "straight", fill: trajectory-purple),
      stroke: (paint: trajectory-purple, thickness: 1.5pt)
    )
    // content(
    //   ((ep + ep-smoothed) / 2, padding),
    //   text(fill: trajectory-purple, size: 0.85em)[small reduction],
    //   anchor: "north"
    // )
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

    // Bight angle indicator at the y-axis (0, eh)
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

  let e-H-smoothed = (0.0, ehs)  // Bemaining error on H
  let e-final = (eps, ehs)       // Final error state
  let e-P-final = (eps, 0.0)

  let midpoint-H = (0.0, (ehs + eh) / 2)

  // Axes
  line((-0.5, 0.0), (xmax, 0.0), stroke: 1pt)
  line((0.0, -0.5), (0.0, ymax), stroke: 1pt)

  content((xmax - padding, 0.0), $range(P)$, anchor: "west")
  content((0.0, ymax - padding), $range(P)^(perp_A)$, anchor: "south")

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
    content(e-final, padding: 0.15, $(I - B A)e_H$, anchor: "west")
  }

  if step == 1{
    content(e-H-smoothed, padding: 0.15, $(I - Pi_c)(I - B A) e_H$, anchor: "east")
  }

  // --- STEP 2: Bed smoother effect along H ---
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
      text(size: 15pt, fill: red-accent)[
        Smoother effect on $range(P)^(perp_A):$ \
        $ X = (I - Pi_c) B A e_H $
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
    let e-H-smoothed = (0, ehs) // Beduced high-frequency component (I - X) e_H
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

    // --- STATE AFTEB SMOOTHEB (I - T_bar)e ---
    // Dashed projection lines for smoothed error
    line(e-smoothed, e-H-smoothed, stroke: (dash: "dashed", paint: gray))
    line(e-smoothed, e-P-smoothed, stroke: (dash: "dashed", paint: gray))

    // Vector (I - T_bar)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - B A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
    content(midpoint, padding: 0.15, $(I - Pi_c)B A e\ = "important part"$, anchor: "east")
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
    let e-H-smoothed = (0, ehs) // Beduced high-frequency component (I - X) e_H
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

    // --- STATE AFTEB SMOOTHEB (I - T_bar)e ---

    // Vector (I - BA)e
    line(e-pt, e-smoothed, mark: (end: "triangle", fill: main-purple), stroke: (paint: main-purple, thickness: 1.5pt)) // effect of smoother
    circle(e-smoothed, radius: 0.06, fill: black, stroke: none)
    content(e-smoothed, padding: 0.15, $(I - B A)e$, anchor: "south-east")

    // Projections of smoothed error on axes
    // circle(e-H-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)
        content(midpoint, padding: 0.15, $(I - Pi_c)B A e_H$, anchor: "east")

    // circle(e-P-smoothed, radius: 0.06, fill: rgb("d9534f"), stroke: none)

    // Action of smoother on H component (reduction arrow)
    line((0, eh), (0, 0), mark: (end: "triangle", fill: red-accent), stroke: (paint: red-accent, thickness: 1.5pt))
    // content((-0.1, 2.2), text(size: 8pt, fill: red-accent)[$text("Smoother action on \n high frequency error")$], anchor: "east")

    content(e-H-smoothed, padding: 0.5, $(I - Pi_c)(I - B A) e_H$, anchor: "east")

    // Vector e
    circle(e-pt, radius: 0.06, fill: black, stroke: none)
    content(e-pt, $e$, padding: 0.15, anchor: "south-west")
  })

  #fig-optimal-smoother