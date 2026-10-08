import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_vline, geom_hline,
    scale_color_manual, labs, theme_minimal, theme, geom_point
)
import os


# -------------------------------------------------------------
# Eigenvalues of S = I - RA for Jacobi, GS, Schwartz
# -------------------------------------------------------------

# 1. Load the spectrum data
df_spec_raw = pd.read_csv('../data/smoother_spectrum.csv')

# Normalize the index axis to [0, 1] for scale invariance
if 'normalized_index' not in df_spec_raw.columns:
    df_spec_raw['normalized_index'] = df_spec_raw['k'] / df_spec_raw['k'].max()

# 2. Reshape from wide to long format
df_spec_plot = df_spec_raw.melt(
    id_vars=['k', 'normalized_index'],
    # value_vars=['Jacobi', 'WeightedJacobi', 'GaussSeidel', 'AdditiveSchwarz'],
    value_vars=['Jacobi','Jacobi_true', 'WeightedJacobi','WeightedJacobi_true', 'GaussSeidel', 'AdditiveSchwarz'],

    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
)

# Clean up display names
smoother_labels = {
    'Jacobi': 'Jacobi (ω=1.0)',
    'WeightedJacobi': 'Weighted Jacobi (ω=2/3)',
    'GaussSeidel': 'Gauss-Seidel',
    'AdditiveSchwarz': 'Additive Schwarz',
    'Jacobi_true': 'Jacobi (predicted)',
    'WeightedJacobi_true': 'Weighted Jacobi (ω=2/3, predicted)',
}
df_spec_plot['Smoother'] = df_spec_plot['Smoother'].map(smoother_labels)

# 3. Custom color palette matching the previous plot
colors = {
    'Jacobi (ω=1.0)': '#e76f51',              # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',       # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f',                 # Teal
    'Additive Schwarz': '#264653',   # Dark slate
    'Jacobi (predicted)': "#99A2AC",
    'Weighted Jacobi (ω=2/3, predicted)': "#C298BC",
}

# 4. Generate plot
n = len(df_spec_raw)

plot_spec = (
    ggplot(df_spec_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8) # Convergence boundary |λ| = 1
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(S)|',
        title='Spectrum of S = ( I - RA) (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        # legend_title=theme(text='Smoother Type'),
    )
)

# Display or save plot
plot_spec

# -------------------------------------------------------------
# Attenuation factor ||Se|| of S = I - RA for Jacobi, GS, Schwartz
# -------------------------------------------------------------

df_raw = pd.read_csv('../data/fourier_attenuation.csv')

# 2. Reshape from wide to long format for plotnine
df_plot = df_raw.melt(
    id_vars=['k', 'normalized_frequency'],
    value_vars=[
        'Jacobi_attenuation', 
        'WeightedJacobi_attenuation', 
        'GaussSeidel_attenuation', 
        'AdditiveSchwarz_attenuation'
    ],
    var_name='Smoother',
    value_name='Attenuation'
)

# Clean up label names
smoother_labels = {
    'Jacobi_attenuation': 'Jacobi (ω=1.0)',
    'WeightedJacobi_attenuation': 'Jacobi (ω=2/3)',
    'GaussSeidel_attenuation': 'Gauss-Seidel',
    'AdditiveSchwarz_attenuation': 'Additive Schwarz'
}
df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)

# 3. Custom color palette
colors = {
    'Jacobi (ω=1.0)': '#e76f51',              # Terracotta
    'Jacobi (ω=2/3)': '#f4a261',       # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f',                 # Teal
    'Additive Schwarz': '#264653'   # Dark slate
}

# 4. Generate plotnine figure
n = 499
high_freq_threshold = 0.5  # Boundary separating low and high frequencies (n/2)

plot = (
    ggplot(df_plot, aes(x='normalized_frequency', y='Attenuation', color='Smoother'))
    + geom_line(size=0.7, alpha=0.85)
    + geom_vline(xintercept=high_freq_threshold, linetype='dashed', color='#8d99ae', size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Frequency Mode',
        y='Attenuation Factor ||S q_k||₂',
        title='Smoother Mode Attenuation across Fourier Spectrum',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        # legend_title=theme(text='Smoother Type'),
    )
)

# Display or save plot
# plot.save("smoother_attenuation.png", dpi=300)
plot

# -------------------------------------------------------------
# Eigenvalues of S_bar = I - R_bar A for Jacobi, GS, Schwartz
# -------------------------------------------------------------

# 1. Load the spectrum data
df_raw = pd.read_csv('../data/smoother_S_spectrum.csv')

# Normalize the index axis to [0, 1] for scale invariance
if 'normalized_index' not in df_raw.columns:
    df_raw['normalized_index'] = df_raw['k'] / df_raw['k'].max()

# 2. Reshape from wide to long format
df_plot = df_raw.melt(
    id_vars=['k', 'normalized_index'],
    value_vars=['Jacobi', 'WeightedJacobi', 'GaussSeidel', 'AdditiveSchwarz'],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
)

# Clean up display names
smoother_labels = {
    'Jacobi': 'Jacobi (ω=1.0)',
    'WeightedJacobi': 'Weighted Jacobi (ω=2/3)',
    'GaussSeidel': 'Gauss-Seidel',
    'AdditiveSchwarz': 'Additive Schwarz'
}
df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)

# 3. Custom color palette matching the previous plot
colors = {
    'Jacobi (ω=1.0)': '#e76f51',              # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',       # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f',                 # Teal
    'Additive Schwarz': '#264653'   # Dark slate
}

# 4. Generate plot
n = len(df_raw)

plot = (
    ggplot(df_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8) # Convergence boundary |λ| = 1
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(S)|',
        title='Spectrum of S_bar = ( I - R_bar A) (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        # legend_title=theme(text='Smoother Type'),
    )
)

# Display or save plot
plot

# -------------------------------------------------------------
# Eigenvalues of X = (I - Pi_c) R_bar A for Jacobi, GS, Schwartz
# -------------------------------------------------------------

# 1. Load the spectrum data
df_raw = pd.read_csv('../data/smoother_X_spectrum.csv')

# Normalize the index axis to [0, 1] for scale invariance
if 'normalized_index' not in df_raw.columns:
    df_raw['normalized_index'] = df_raw['k'] / df_raw['k'].max()

# 2. Reshape from wide to long format
df_plot = df_raw.melt(
    id_vars=['k', 'normalized_index'],
    value_vars=['Jacobi', 'WeightedJacobi', 'GaussSeidel', 'AdditiveSchwarz'],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
)

# Clean up display names
smoother_labels = {
    'Jacobi': 'Jacobi (ω=1.0)',
    'WeightedJacobi': 'Weighted Jacobi (ω=2/3)',
    'GaussSeidel': 'Gauss-Seidel',
    'AdditiveSchwarz': 'Additive Schwarz'
}
df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)

# 3. Custom color palette matching the previous plot
colors = {
    'Jacobi (ω=1.0)': '#e76f51',              # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',       # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f',                 # Teal
    'Additive Schwarz': '#264653'   # Dark slate
}

# 4. Generate plot
n = len(df_raw)

plot = (
    ggplot(df_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8) # Convergence boundary |λ| = 1
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(S)|',
        title='Spectrum of X with linear interpolation P (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | dim(H) = n - n_c = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        # legend_title=theme(text='Smoother Type'),
    )
)

# Display or save plot
plot

# -------------------------------------------------------------
# Eigenvalues of optimal X = (I - Pi_c_optimal) R_bar A for Jacobi, GS, Schwartz
# -------------------------------------------------------------

# 1. Load the spectrum data
df_raw = pd.read_csv('../data/smoother_X_opt_spectrum.csv')

# Normalize the index axis to [0, 1] for scale invariance
if 'normalized_index' not in df_raw.columns:
    df_raw['normalized_index'] = df_raw['k'] / df_raw['k'].max()

# 2. Reshape from wide to long format
df_plot = df_raw.melt(
    id_vars=['k', 'normalized_index'],
    # value_vars=['Jacobi', 'WeightedJacobi', 'GaussSeidel', 'AdditiveSchwarz'],
    value_vars=['Jacobi','Jacobi_true', 'WeightedJacobi','WeightedJacobi_true', 'GaussSeidel', 'AdditiveSchwarz'],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
)

# Clean up display names
smoother_labels = {
    'Jacobi': 'Jacobi (ω=1.0)',
    'WeightedJacobi': 'Weighted Jacobi (ω=2/3)',
    'GaussSeidel': 'Gauss-Seidel',
    'AdditiveSchwarz': 'Additive Schwarz',
    'Jacobi_true': 'Jacobi (predicted)',
    'WeightedJacobi_true': 'Weighted Jacobi (ω=2/3, predicted)',
}

drawing_order = [
    'Jacobi (predicted)',
    'Jacobi (ω=1.0)',
    'Weighted Jacobi (ω=2/3, predicted)',
    'Weighted Jacobi (ω=2/3)',
    'Gauss-Seidel',
    'Additive Schwarz'  # Plotted last -> rendered on top
]

# 2. Map labels
df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)

# 3. Convert to Categorical and sort the DataFrame rows
df_plot['Smoother'] = pd.Categorical(
    df_plot['Smoother'], 
    categories=drawing_order, 
    ordered=True
)
df_plot = df_plot.sort_values('Smoother')

# 3. Custom color palette matching the previous plot
colors = {
    'Jacobi (ω=1.0)': '#e76f51',              # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',       # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f',                 # Teal
    'Additive Schwarz': '#264653',   # Dark slate
    'Jacobi (predicted)': "#99A2AC",
    'Weighted Jacobi (ω=2/3, predicted)': "#71777E",
}

# 4. Generate plot
n = len(df_raw)

plot = (
    ggplot(df_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8) # Convergence boundary |λ| = 1
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(S)|',
        title='Spectrum of X with optimal interpolation P (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | dim(H) = n - n_c = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        # legend_title=theme(text='Smoother Type'),
    )
)

# Display or save plot
plot


# -------------------------------------------------------------
# Eigenvalue of X for Jacobi
# -------------------------------------------------------------

# Parameters
n = 49                  # Grid size / total modes
n_c = 24                 # Number of coarse modes kept
omega = 2.0 / 3.0        # Optimal weighted Jacobi parameter
k_c = (n_c + 1) / (n + 1) # Normalized cutoff frequency

# Continuous mode parameter k_norm = k / (n + 1) in (0, 1]
k_norm = np.linspace(1 / (n + 1), 1.0, 1000)

# 1. Eigenvalues of A: lambda_A(k) = 4 * sin^2(k_norm * pi / 2)
lambda_A = 4.0 * (np.sin(k_norm * np.pi / 2.0) ** 2)

# 2. Eigenvalues of T_bar = R_bar * A where R_bar = R(2I - AR), R = (omega/2)I
# lambda_Tbar = (omega/2) * lambda_A * (2 - (omega/2) * lambda_A)
lambda_Tb = (omega / 2.0) * lambda_A * (2.0 - (omega / 2.0) * lambda_A)

# 3. Eigenvalues of X = T_bar * (I - Pi_c)
# Modes k <= n_c are projected out to 0; modes k > n_c retain lambda_Tb
lambda_X = np.where(k_norm <= k_c, 0.0, lambda_Tb)

# Prepare DataFrame for plotnine
df = pd.DataFrame({
    'k_norm': np.tile(k_norm, 2),
    'Eigenvalue': np.concatenate([lambda_Tb, lambda_X]),
    'Operator': (
        ['T_bar (Smoother)'] * len(k_norm) +
        ['X = T_bar(I - Pi_c)'] * len(k_norm)
    )
})

# Color palette matching operator roles
colors = {
    'T_bar (Smoother)': '#e07a5f',     # Terracotta
    'X = T_bar(I - Pi_c)': '#3d405b'   # Dark slate
}

# # -------------------------------------------------------------
# # Eigenvalue of X for GS
# # -------------------------------------------------------------

# # Load spectrum data generated by Julia
# df_raw = pd.read_csv('../data/spectrum_X_GS.csv')

# n = df_raw.shape[0]      # Grid size / total modes
# n_c = 250                 # Number of coarse modes kept
# k_c = (n_c + 1) / (n + 1)

# # Melt DataFrame for plotnine grouped representation
# df_plot = pd.melt(
#     df_raw,
#     id_vars=['k_mode', 'k_norm'],
#     value_vars=['Modal', 'GMG', 'XZ_Optimal'],
#     var_name='Coarse_Space',
#     value_name='Eigenvalue'
# )

# # Color palette for coarse spaces
# colors = {
#     'Modal': '#e76f51',       # Coral / Terracotta
#     'GMG': '#2a9d8f',         # Teal
#     'XZ_Optimal': '#264653'   # Dark slate
# }

# # Plot restricted spectrum on H^A
# plot = (
#     ggplot(df_plot, aes(x='k_norm', y='Eigenvalue', color='Coarse_Space'))
#     # + geom_line(size=1.1)
#     + geom_point(size = 0.7, alpha=0.8)
#     + geom_vline(xintercept=k_c, linetype='dashed', color='#8d99ae', size=0.8)
#     + scale_color_manual(values=colors)
#     + labs(
#         x='Normalized Frequency Mode',
#         y='Eigenvalue Magnitude λ(X)',
#         title='Restricted Spectrum of Two-Grid Operator X on H',
#         subtitle=f'Gauss-Seidel Smoother | n = {n}, n_c = {n_c}'
#     )
#     + theme_minimal()
#     + theme(
#         figure_size=(9, 5),
#         legend_position='right',
#     )
# )

# # Save or show plot
# # plot.save('spectrum_X_comparison.png', dpi=300)
# plot

# # -------------------------------------------------------------
# # Eigenvalue of X for Additive Schwartz
# # -------------------------------------------------------------

df_raw = pd.read_csv('../data/spectrum_X_Schwartz.csv')

# 2. Reshape from wide to long format
df_plot = df_raw.melt(
    id_vars=['k_mode', 'k_norm'],
    value_vars=['Modal', 'GMG', 'XZ_Optimal'],
    var_name='Interpolation',
    value_name='Eigenvalue_Magnitude'
)

# Clean up display labels
interp_labels = {
    'Modal': 'Modal P',
    'GMG': 'Standard GMG P',
    'XZ_Optimal': 'XZ-Optimal P'
}

drawing_order = [
    'Standard GMG P',
    'Modal P',
    'XZ-Optimal P'  # Plotted last -> rendered on top
]

df_plot['Interpolation'] = df_plot['Interpolation'].map(interp_labels)
df_plot['Interpolation'] = pd.Categorical(
    df_plot['Interpolation'], 
    categories=drawing_order, 
    ordered=True
)
df_plot = df_plot.sort_values('Interpolation')

# 3. Custom palette & line styles
colors = {
    'Standard GMG P': '#e76f51',    # Terracotta
    'Modal P': '#2a9d8f',           # Teal
    'XZ-Optimal P': '#264653',      # Dark Slate
}

# 4. Generate plot
n = len(df_raw)

plot = (
    ggplot(df_plot, aes(
        x='k_norm', 
        y='Eigenvalue_Magnitude', 
        color='Interpolation'
    ))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color='#4e4774', size=0.8) # Convergence boundary |λ| = 1
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Mode Index (k / n)',
        y='Eigenvalue Magnitude |λ(X)|',
        title='Additive Schwarz Spectrum Across Interpolation Operators P',
        subtitle=f'1D Laplacian | Problem Dimension n = {n}',
        color='Interpolation P',
        linetype='Interpolation P'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
    )
)

plot
