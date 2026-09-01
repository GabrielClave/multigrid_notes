# brought to you by llm

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_text,
    facet_grid,
    geom_hline,
    geom_line,
    ggplot,
    labs,
    scale_color_manual,
    theme,
    theme_minimal,
    geom_point,
    facet_wrap,
    scale_x_continuous,
    scale_y_log10,
    theme_minimal,
    element_blank,
    element_rect,
    ggsave
)
import scipy.sparse as sp
import scipy.sparse.linalg as spla

bg_color = "#fafafa"

# ==============================================================================
# Smoother effect 
# ==============================================================================

N = 100
h = 1.0 / (N + 1)
x = np.linspace(h, 1.0 - h, N)

main_diag = 2.0 * np.ones(N) / (h**2)
off_diag = -1.0 * np.ones(N - 1) / (h**2)
A = sp.diags([off_diag, main_diag, off_diag], [-1, 0, 1], format="csr")

omega = 2 / 3
D_inv = 1.0 / (2.0 / (h**2))


def damped_jacobi(A, b, x0, iterations=10, omega=2 / 3):
  x = x0.copy()
  for _ in range(iterations):
    r = b - A @ x
    x = x + omega * D_inv * r
  return x


# 2. Modes Setup
b = np.zeros(N)
v_high = np.sin(20 * np.pi * x)
v_low = np.sin(2 * np.pi * x)
v_comp = v_high + v_low

num_iters = 10
v_high_smooth = damped_jacobi(A, b, v_high, iterations=num_iters, omega=omega)
v_low_smooth = damped_jacobi(A, b, v_low, iterations=num_iters, omega=omega)
v_comp_smooth = damped_jacobi(A, b, v_comp, iterations=num_iters, omega=omega)

# 3. Assemble Tidy Dataframe for plotnine
df_list = []
signals = [
    ("1. High Freq", v_high, v_high_smooth),
    ("2. Low Freq", v_low, v_low_smooth),
    # ("3. Composite", v_comp, v_comp_smooth),
]

for signal_name, initial, smoothed in signals:
  df_list.append(
      pd.DataFrame({
          "x": x,
          "value": initial,
          "state": "Initial",
          "signal": signal_name,
      })
  )
  df_list.append(
      pd.DataFrame({
          "x": x,
          "value": smoothed,
          "state": f"After {num_iters} iterations",
          "signal": signal_name,
      })
  )

df = pd.concat(df_list, ignore_index=True)

# Preserve order in facet panels and legends
df["signal"] = pd.Categorical(
    df["signal"], categories=[s[0] for s in signals], ordered=True
)
df["state"] = pd.Categorical(
    df["state"],
    categories=["Initial", f"After {num_iters} iterations"],
    ordered=True,
)

# 4. Plot using Plotnine Syntax
plot = (
    ggplot(df, aes(x="x", y="value", color="state"))
    + geom_hline(yintercept=0, linetype="dashed", color="#cccccc")
    + geom_line(size=0.9)
    + facet_grid("signal ~ state", scales="free_y")
    + scale_color_manual(values=["#d95f02", "#7570b3"], name="Smoother State")
    + labs(
        title=f"Damped Jacobi Smoothing Effect",
        x="Grid coordinate x",
        y="Amplitude",
    )
    + theme_minimal()
    + theme(
        figure_size=(10, 7),
        legend_position="none",
        plot_background=element_rect(fill=bg_color, color=bg_color),
        strip_text_x=element_text(size=10, weight="bold"),
        strip_text_y=element_text(size=10, weight="bold"),
        title=element_text(size=13, weight="bold"),
    )
)

plot

ggsave(plot, filename="damped_jacobi.png", dpi=300, width=10, height=7)

# ==============================================================================
###### Jacobi used as a solver
# ==============================================================================

# Right-hand side b and exact solution x_star
b = np.ones(N)
x_star = spla.spsolve(A, b)

# 3. Solver Iterations tracking error ||e_k||_2
np.random.seed(42)
x_k = np.random.randn(N)  # Random initial guess with high + low frequencies

max_iters = 150
iterations = []
error_norms = []
rel_error_norms = []

for k in range(max_iters + 1):
  e_k = x_star - x_k
  err_norm = np.linalg.norm(e_k, 2)

  iterations.append(k)
  error_norms.append(err_norm)

  if k < max_iters:
    r_k = b - A @ x_k
    x_k = x_k + omega * D_inv * r_k

df_conv = pd.DataFrame({
    "iteration": iterations,
    "error_norm": error_norms,
    "rel_error_norm": np.array(error_norms) / error_norms[0],
})

# 4. Plot Convergence Curve (Semi-log scale)
conv_plot = (
    ggplot(df_conv, aes(x="iteration", y="rel_error_norm"))
    + geom_line(color="#7570b3", size=1.0)
    + geom_point(
        data=df_conv[df_conv["iteration"] % 10 == 0],
        color="#7570b3",
        size=2.0,
    )
    + scale_y_log10()
    + labs(
        title=f"Jacobi Solver Convergence",
        x="Iteration",
        y="Relative Error Norm (log scale)",
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        plot_background=element_rect(fill=bg_color, color=bg_color),
        strip_text_x=element_text(size=10, weight="bold"),
        title=element_text(size=12, weight="bold"),
    )
)

conv_plot

ggsave(conv_plot, filename="jacobi_convergence.png", dpi=300, width=10, height=7)

# ==============================================================================
###### coarse grid effect
# ==============================================================================

# 1. Setup Grids
N_h = 31  # Fine grid interior points
x_h = np.linspace(0, 1, N_h + 2)  # Include boundary points 0 and 1
x_H = x_h[::2]  # Coarse grid: 1 point every 2 fine points

# High-resolution grid for smooth background sine wave
x_dense = np.linspace(0, 1, 500)

# Mode definition: Oscillatory mode k = 24 (high frequency on fine grid)
k = 8
v_dense = np.sin(k * np.pi * x_dense)
v_h = np.sin(k * np.pi * x_h)
v_H = np.sin(k * np.pi * x_H)

# Interpolation back to fine grid
v_interp = np.interp(x_h, x_H, v_H)

# ==============================================================================
# Plot 1: Sampling aligned with continuous sine in the background
# ==============================================================================
# Data for discrete grid points
df_discrete = pd.concat([
    pd.DataFrame({
        "x": x_h,
        "y": v_h,
        "grid": f"1. Fine Grid (h): {len(x_h)} points",
    }),
    pd.DataFrame({
        "x": x_H,
        "y": v_H,
        "grid": f"2. Coarse Grid (H): {len(x_H)} points",
    }),
])

df_discrete["grid"] = pd.Categorical(
    df_discrete["grid"],
    categories=[
        f"1. Fine Grid (h): {len(x_h)} points",
        f"2. Coarse Grid (H): {len(x_H)} points",
    ],
    ordered=True,
)

# Data for continuous background sine wave across both facets
df_background = pd.concat([
    pd.DataFrame({
        "x": x_dense,
        "y": v_dense,
        "grid": f"1. Fine Grid (h): {len(x_h)} points",
    }),
    pd.DataFrame({
        "x": x_dense,
        "y": v_dense,
        "grid": f"2. Coarse Grid (H): {len(x_H)} points",
    }),
])

df_background["grid"] = pd.Categorical(
    df_background["grid"],
    categories=[
        f"1. Fine Grid (h): {len(x_h)} points",
        f"2. Coarse Grid (H): {len(x_H)} points",
    ],
    ordered=True,
)

plot_sampling = (
    ggplot()
    + geom_hline(yintercept=0, linetype="dashed", color="#aaaaaa")
    # Continuous background sine wave
    + geom_line(
        df_background,
        aes(x="x", y="y"),
        color="#666464ac",
        size=0.8,
    )
    # Grid points
    + geom_point(df_discrete, aes(x="x", y="y", color="grid"), size=2.0)
    + facet_wrap("~grid", ncol=1)
    + scale_color_manual(values=["#d95f02", "#7570b3"], guide=None)
    + scale_x_continuous(breaks=[0, 0.25, 0.5, 0.75, 1.0], expand=(0, 0))
    + labs(
        title=f"Coarse Grid Sampling",
        x="Spatial Coordinate x",
        y="Amplitude",
    )
    + theme_minimal()
    + theme(
        figure_size=(10, 7),
        plot_background=element_rect(fill=bg_color, color=bg_color),
        strip_text=element_text(size=11, weight="bold"),
        title=element_text(size=13, weight="bold"),
        panel_grid_major_x=element_text(color="#eeeeee"),
        panel_spacing_y=0.08,
    )
)

plot_sampling

ggsave(plot_sampling, filename="signal_sampling.png", dpi=300, width=10, height=7)


# ==============================================================================
# Plot 2: Superposition showing Aliasing / Frequency Misalignment
# ==============================================================================
# Mapping dense x coordinates [0, 1] to index space for background plotting
idx_dense_h = x_dense * (len(x_h) - 1)
idx_dense_H = x_dense * (len(x_H) - 1)

df_super_bg = pd.concat([
    pd.DataFrame({
        "point_index": idx_dense_h,
        "y": v_dense,
        "type": "Fine Grid Signal (h)",
    }),
    pd.DataFrame({
        "point_index": idx_dense_H,
        "y": v_dense,
        "type": "Coarse Grid Signal (H - stops at index 16)",
    }),
])

df_super_points = pd.concat([
    pd.DataFrame({
        "point_index": np.arange(len(x_h)),
        "y": v_h,
        "type": "Fine Grid Signal (h)",
    }),
    pd.DataFrame({
        "point_index": np.arange(len(x_H)),
        "y": v_H,
        "type": "Coarse Grid Signal (H - stops at index 16)",
    }),
])

color_fine_line = "#daa67e"  # Background wave color for fine grid (e.g., orange)
color_coarse_line = "#7570b3"  # Background wave color for coarse grid (e.g., purple)

color_fine_point = "#d95f02"  # Discrete point color for fine grid (red)
color_coarse_point = "#7570b3"  # Discrete point color for coarse grid (blue)

plot_superposition_index = (
    ggplot()
    + geom_hline(yintercept=0, linetype="dashed", color="#aaaaaa")
    # Smooth continuous sine waves with independent line colors
    + geom_line(
        df_super_bg[df_super_bg["type"] == "Fine Grid Signal (h)"],
        aes(x="point_index", y="y"),
        color=color_fine_line,
        size=0.9,
        alpha = 0.5
    )
    + geom_line(
        df_super_bg[
            df_super_bg["type"]
            == "Coarse Grid Signal (H - stops at index 16)"
        ],
        aes(x="point_index", y="y"),
        color=color_coarse_line,
        size=0.9,
        alpha = 0.5
    )
    # Discrete grid points with their own mapped colors for the legend
    + geom_point(
        df_super_points,
        aes(x="point_index", y="y", color="type"),
        size=2.2,
    )
    + scale_color_manual(
        values=[color_coarse_point, color_fine_point],
        name="Grid Points",
    )
    + labs(
        title=f"Superposition by Index: Frequency Doubling Effect (k = {k})",
        x="Grid Point Index (i)",
        y="Amplitude",
    )
    + theme_minimal()
    + theme(
        figure_size=(10, 5),
        plot_background=element_rect(fill=bg_color, color=bg_color),
        title=element_text(size=13, weight="bold"),
        legend_position="bottom",
        legend_title=element_text(size=10, weight="bold"),
        panel_grid_major_x=element_text(color="#eeeeee"),
    )
)

plot_superposition_index

ggsave(plot_superposition_index, filename="coarse_restriction.png", dpi=300, width=10, height=5)

# ==============================================================================
###### Two-level Grid operator
# ==============================================================================

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_text,
    facet_grid,
    geom_hline,
    geom_line,
    ggplot,
    labs,
    scale_color_manual,
    theme,
    theme_minimal,
)

# 1. Discretization Setup (1D Poisson using dense arrays)
N_h = 63  # Fine grid interior points (2^k - 1)
h = 1.0 / (N_h + 1)
x_h = np.linspace(h, 1.0 - h, N_h)

# 1D Laplace Matrix A_h (N_h x N_h)
main_diag = 2.0 * np.ones(N_h) / (h**2)
off_diag = -1.0 * np.ones(N_h - 1) / (h**2)
A_h = (
    np.diag(main_diag) + np.diag(off_diag, k=-1) + np.diag(off_diag, k=1)
)

# Damped Jacobi Smoother Setup (S = I - omega * D^-1 * A)
omega = 2 / 3
D_inv = h**2 / 2.0


def apply_jacobi(x, q):
  x_out = x.copy()
  for _ in range(q):
    r = -A_h @ x_out  # Homogeneous system b = 0
    x_out += omega * D_inv * r
  return x_out


# 2. Dense Grid Transfer Operators
N_H = (N_h - 1) // 2  # Coarse grid interior points (31)

# Full-weighting Restriction R (N_H x N_h)
R = np.zeros((N_H, N_h))
for i in range(N_H):
  fine_idx = 2 * i + 1
  R[i, fine_idx - 1] = 0.25
  R[i, fine_idx] = 0.50
  R[i, fine_idx + 1] = 0.25

# Linear Prolongation P (N_h x N_H)
P = 2.0 * R.T

# Galerkin Coarse Grid Operator A_H (N_H x N_H)
A_H = R @ A_h @ P


# 3. Two-Grid Error Propagation Function T_q(e)
def apply_two_grid(e_init, q_presmoothing):
  # Step 1: Pre-smoothing (q iterations)
  e_smoothed = apply_jacobi(e_init, q_presmoothing)

  # Step 2: Residual computation (b = 0, so r = -A*e)
  r_h = -A_h @ e_smoothed

  # Step 3: Restriction to coarse grid
  r_H = R @ r_h

  # Step 4: Direct solve on coarse grid using standard dense solver
  e_H = np.linalg.solve(A_H, r_H)

  # Step 5: Prolongation & Coarse Grid Correction
  e_corrected = e_smoothed + P @ e_H

  return e_corrected


# 4. Define Modes & Apply T_q
k_low, k_high = 2, 10
v_low = np.sin(k_low * np.pi * x_h)
v_high = np.sin(k_high * np.pi * x_h)
v_comp = v_low + v_high

signals = [
    ("1. Low Freq (k=2)", v_low),
    ("2. High Freq (k=10)", v_high),
    ("3. Composite Signal", v_comp),
]

q_values = [0, 1, 20]

df_list = []
for sig_name, v_init in signals:
  # Initial state
  df_list.append(
      pd.DataFrame({
          "x": x_h,
          "error": v_init,
          "state": "Initial Error",
          "signal": sig_name,
      })
  )

  # Apply T_q for different q
  for q in q_values:
    e_tg = apply_two_grid(v_init, q_presmoothing=q)
    df_list.append(
        pd.DataFrame({
            "x": x_h,
            "error": e_tg,
            "state": f"TG (q={q})",
            "signal": sig_name,
        })
    )

df = pd.concat(df_list, ignore_index=True)

# Preserve facet ordering
df["signal"] = pd.Categorical(
    df["signal"], categories=[s[0] for s in signals], ordered=True
)
df["state"] = pd.Categorical(
    df["state"],
    categories=["Initial Error", "TG (q=0)", "TG (q=1)", "TG (q=3)"],
    ordered=True,
)

# 5. Plot using plotnine
plot_tg = (
    ggplot(df, aes(x="x", y="error", color="state"))
    + geom_hline(yintercept=0, linetype="dashed", color="#cccccc")
    + geom_line(size=0.8)
    + facet_grid("signal ~ state", scales="free_y")
    + scale_color_manual(
        values=["#000000", "#e41a1c", "#377eb8", "#4daf4a"],
        name="Operator Stage",
    )
    + labs(
        title="Effect of Dense Two-Grid Operator TG_q = (I - P A_H^-1 R A) S^q",
        x="Spatial Coordinate x",
        y="Error Amplitude",
    )
    + theme_minimal()
    + theme(
        figure_size=(11, 7),
        strip_text_x=element_text(size=10, weight="bold"),
        strip_text_y=element_text(size=10, weight="bold"),
        title=element_text(size=12, weight="bold"),
        legend_position="bottom",
    )
)

plot_tg

