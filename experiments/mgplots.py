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
    theme_bw,
    geom_point,
    facet_wrap,
    scale_x_continuous
)
import scipy.sparse as sp

###### Smoother effect 

# 1. Discretization Setup
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
    ("1. High Freq (k=20)", v_high, v_high_smooth),
    ("2. Low Freq (k=2)", v_low, v_low_smooth),
    ("3. Composite (k=2 + k=20)", v_comp, v_comp_smooth),
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
        title=f"Damped Jacobi Smoothing Effect (omega = {omega:.2f})",
        x="Grid coordinate x",
        y="Amplitude",
    )
    + theme_bw()
    + theme(
        figure_size=(10, 7),
        strip_text_x=element_text(size=10, weight="bold"),
        strip_text_y=element_text(size=10, weight="bold"),
        title=element_text(size=13, weight="bold"),
    )
)

plot

###### coarse grid effect

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
    + scale_color_manual(values=["#e41a1c", "#377eb8"], guide=None)
    + scale_x_continuous(breaks=[0, 0.25, 0.5, 0.75, 1.0], expand=(0, 0))
    + labs(
        title=f"Coarse Grid Sampling & Aliasing Effect (k = {k})",
        x="Spatial Coordinate x",
        y="Amplitude",
    )
    + theme_bw()
    + theme(
        figure_size=(10, 5),
        strip_text=element_text(size=11, weight="bold"),
        title=element_text(size=13, weight="bold"),
        panel_grid_major_x=element_text(color="#eeeeee"),
        panel_spacing_y=0.08,
    )
)

plot_sampling


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

color_fine_line = "#cf4d4dff"  # Background wave color for fine grid (e.g., orange)
color_coarse_line = "#6481a77b"  # Background wave color for coarse grid (e.g., purple)

color_fine_point = "#e41a1c"  # Discrete point color for fine grid (red)
color_coarse_point = "#377eb8"  # Discrete point color for coarse grid (blue)

plot_superposition_index = (
    ggplot()
    + geom_hline(yintercept=0, linetype="dashed", color="#aaaaaa")
    # Smooth continuous sine waves with independent line colors
    + geom_line(
        df_super_bg[df_super_bg["type"] == "Fine Grid Signal (h)"],
        aes(x="point_index", y="y"),
        color=color_fine_line,
        size=0.9,
    )
    + geom_line(
        df_super_bg[
            df_super_bg["type"]
            == "Coarse Grid Signal (H - stops at index 16)"
        ],
        aes(x="point_index", y="y"),
        color=color_coarse_line,
        size=0.9,
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
    + theme_bw()
    + theme(
        figure_size=(10, 5),
        title=element_text(size=13, weight="bold"),
        legend_position="bottom",
        legend_title=element_text(size=10, weight="bold"),
        panel_grid_major_x=element_text(color="#eeeeee"),
    )
)

plot_superposition_index