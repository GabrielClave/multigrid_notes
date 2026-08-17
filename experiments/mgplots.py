import matplotlib.pyplot as plt

# Force crisp, high-DPI inline rendering in Jupyter interactive window
plt.rcParams["figure.dpi"] = 150

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
)
import scipy.sparse as sp

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