import pandas as pd
import plotnine as p9
import numpy as np

# ==============================================================================
# LOAD DATA
# ==============================================================================
df_components = pd.read_csv("data/df_components.csv")
df_2d = pd.read_csv("data/df_2d.csv")
df_lines = pd.read_csv("data/df_lines.csv")
# 2. Take absolute value of components
df_components["Component"] = np.abs(df_components["Component"])

# 3. Ensure correct categorical ordering
stage_order = ["1. Initial (e0)", "2. After Smoother (e1)", "3. After CGC (e2)"]
df_components["Stage"] = pd.Categorical(
    df_components["Stage"], categories=stage_order, ordered=True
)

basis_order = ["Euclidean", "Eigenvector"]
df_components["Basis"] = pd.Categorical(
    df_components["Basis"], categories=basis_order, ordered=True
)

# ==============================================================================
# PLOT 1: error coordinates
# ==============================================================================

# 4. Determine uniform Y-axis maximum (with ~5% padding)
y_max = df_components["Component"].max() * 1.05

# 5. Build Plot
plot1 = (
    p9.ggplot(df_components, p9.aes(x="Index", y="Component", fill="Stage"))
    + p9.geom_col(show_legend=False, width=0.7)
    + p9.facet_grid("Stage ~ Basis")  # Shared Y scale across all facets
    + p9.scale_y_continuous(limits=(0, y_max))
    + p9.scale_fill_manual(values=["#2b5c8f", "#d95f02", "#7570b3"])
    + p9.theme_minimal()
    + p9.theme(
        # Text sizes
        title=p9.element_text(size=16, weight="bold"),
        axis_title=p9.element_text(size=14, weight="bold"),
        axis_text=p9.element_text(size=12),
        strip_text=p9.element_text(size=13, weight="bold"),
        # Clean up layout
        panel_grid_minor=p9.element_blank(),
        figure_size=(10, 6),
    )
    + p9.labs(
        x="Component Index (1 to n)",
        y="|Error Component|",
        title="Error Components: Standard Euclidean vs. Eigenbasis",
    )
)

# Save
# plot1.save("../pictures/smoother_error_component.png", dpi=300)

# ==============================================================================
# PLOT 2: 2D Subspace Norm Trajectory
# ==============================================================================

# 2. Axis limits: Extend from slightly negative to ~1.5x max
x_max = 1.5*df_2d["x"].max()
x_min = -0.1*x_max

y_max = 1.5*df_2d["y"].max()
y_min = -0.1*y_max

# 3. Build Plot
plot2 = (
    p9.ggplot()
    # Explicit X and Y axis reference lines at 0
    + p9.geom_hline(yintercept=0, color="black", size=0.8, linetype="solid")
    + p9.geom_vline(xintercept=0, color="black", size=0.8, linetype="solid")
    + p9.scale_x_continuous(limits=(x_min, x_max))
    + p9.scale_y_continuous(limits=(y_min, y_max))
    + p9.scale_color_manual(
        values=["#2b5c8f", "#d95f02", "#1b9e77", "#7570b3"]
    )
    + p9.theme_minimal()
    + p9.theme(
        # Text sizes
        title=p9.element_text(size=16, weight="bold"),
        axis_title=p9.element_text(size=14, weight="bold"),
        axis_text=p9.element_text(size=12),
        legend_title=p9.element_text(size=13, weight="bold"),
        legend_text=p9.element_text(size=12),
        # Clean up layout
        panel_grid_minor=p9.element_blank(),
        figure_size=(10, 6),
    )
    # Connecting trajectory segments
    + p9.geom_segment(
        df_lines, # type: ignore
        p9.aes(x="x", y="y", xend="xend", yend="yend"), # type: ignore
        linetype="dotted",
        color="gray",
        size=1.0,
    )
    # Stage points
    + p9.geom_point(
        df_2d, # type: ignore
        p9.aes(x="x", y="y", color="Stage"), # type: ignore
        size=5,
    )
    + p9.labs(
        x="High-Frequency Subspace Norm",
        y="Low-Frequency Subspace Norm",
        title="Error Trajectory in 2D Norm Space",
        color="State",
    )
)

# Save Plot 2
# plot2.save("../pictures/plot2_subspace_trajectory.png", dpi=300)

# ==============================================================================
# PLOT 3: 2D Subspace Norm Trajectory with Jacobi and Gauss-Seidel
# ==============================================================================

# 1. Load Data
df_2d = pd.read_csv("data/df_2d_real.csv")
df_lines = pd.read_csv("data/df_lines_real.csv")

# Updated stage category list (without e_cgc)
stage_order = [
    "e0 (Initial)",
    "e1 (After Smoother 1)",
    "e2 (After CGC 1)",
    "e3 (After Smoother 2)",
    "e4 (After CGC 2)",
]
df_2d["State"] = pd.Categorical(
    df_2d["State"], categories=stage_order, ordered=True
)

# ------------------------------------------------------------------------------
# PLOT 2 COMBINED (Cycle 1 Only)
# ------------------------------------------------------------------------------

# 1. Filter Data & Lines for Cycle 1
cycle1_states = ["e0 (Initial)", "e1 (After Smoother 1)", "e2 (After CGC 1)"]
df_2d_c1 = df_2d[df_2d["State"].isin(cycle1_states)]

# Keep only the first two segments per smoother (e0 -> e1 and e1 -> e2)
df_lines_c1 = df_lines.groupby("Smoother").head(2)

# Re-calculate limits based on Cycle 1 data
x_max = 1.5 * df_2d_c1["x"].max()
x_min = -0.1 * x_max
y_max = 1.5 * df_2d_c1["y"].max()
y_min = -0.1 * y_max

# 2. Build Plot
plot2_combined = (
    p9.ggplot()
    + p9.geom_hline(yintercept=0, color="black", size=0.8)
    + p9.geom_vline(xintercept=0, color="black", size=0.8)
    + p9.geom_segment(
        df_lines_c1, # type: ignore
        p9.aes(x="x", y="y", xend="xend", yend="yend", linetype="Smoother"), # type: ignore
        color="gray",
        size=0.9,
    )
    + p9.geom_point(
        df_2d_c1, # type: ignore
        p9.aes(x="x", y="y", color="State", shape="Smoother"), # type: ignore
        size=5,
    )
    + p9.scale_x_continuous(limits=(x_min, x_max))
    + p9.scale_y_continuous(limits=(y_min, y_max))
    + p9.scale_color_manual(values=["#2b5c8f", "#d95f02", "#7570b3"])
    + p9.theme_minimal()
    + p9.theme(
        title=p9.element_text(size=16, weight="bold"),
        axis_title=p9.element_text(size=14, weight="bold"),
        axis_text=p9.element_text(size=12),
        legend_title=p9.element_text(size=13, weight="bold"),
        legend_text=p9.element_text(size=12),
        panel_grid_minor=p9.element_blank(),
        figure_size=(9, 7),
    )
    + p9.labs(
        x="High-Frequency Subspace Norm ||P_H e||",
        y="Low-Frequency Subspace Norm ||P_L e||",
        title="2D Discrete Laplacian: Damped Jacobi vs. Gauss-Seidel (Cycle 1)",
        color="State",
        shape="Smoother",
        linetype="Smoother",
    )
)

plot2_combined.save("../pictures/plot2_real_smoothers.png", dpi=300)

# 2. Build Plot 3
plot3 = (
    p9.ggplot()
    + p9.geom_hline(yintercept=0, color="black", size=0.8)
    + p9.geom_vline(xintercept=0, color="black", size=0.8)
    # Trajectory lines
    + p9.geom_segment(
        df_lines, # type: ignore
        p9.aes(x="x", y="y", xend="xend", yend="yend", linetype="Smoother"), # type: ignore
        color="gray",
        size=0.9,
    )
    # Trajectory points
    + p9.geom_point(
        df_2d, # type: ignore
        p9.aes(x="x", y="y", color="State", shape="Smoother"), # type: ignore
        size=4.5,
    )
    + p9.scale_x_log10()
    + p9.scale_y_log10()
    + p9.scale_color_manual(
        values=[
            "#2b5c8f",
            "#d95f02",
            "#1b9e77",
            "#7570b3",
            "#e7298a",
            "#66a61e",
        ]
    )
    + p9.theme_minimal()
    + p9.theme(
        title=p9.element_text(size=16, weight="bold"),
        axis_title=p9.element_text(size=14, weight="bold"),
        axis_text=p9.element_text(size=12),
        legend_title=p9.element_text(size=13, weight="bold"),
        legend_text=p9.element_text(size=11),
        panel_grid_minor=p9.element_blank(),
        figure_size=(9, 7),
    )
    + p9.labs(
        x="High-Frequency Subspace Norm ||P_H e||",
        y="Low-Frequency Subspace Norm ||P_L e||",
        title="2-Grid Iteration Trajectory Across 2 Cycles",
        color="State",
        shape="Smoother",
        linetype="Smoother",
    )
)

# 3. Save
plot3.save("plot3_two_cycles.png", dpi=300)
print("Saved plot3_two_cycles.png")