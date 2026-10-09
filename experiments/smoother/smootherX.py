import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_vline, geom_hline,
    scale_color_manual, labs, theme_minimal, theme, geom_point, scale_y_log10, scale_linetype_manual
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
n = df_raw.shape[0]
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
# Attenuation factor ||Se||_A of S = I - RA for Jacobi, GS, Schwartz
# -------------------------------------------------------------

df_raw = pd.read_csv('../data/fourier_attenuation_anorm.csv')

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
n = df_raw.shape[0]
high_freq_threshold = 0.5  # Boundary separating low and high frequencies (n/2)

plot = (
    ggplot(df_plot, aes(x='normalized_frequency', y='Attenuation', color='Smoother'))
    + geom_line(size=0.7, alpha=0.85)
    + geom_vline(xintercept=high_freq_threshold, linetype='dashed', color='#8d99ae', size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Frequency Mode',
        y='Attenuation Factor ||S q_k||_A / ||q_k||_A',
        title='Smoother Mode Attenuation across Fourier Spectrum in A norm',
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

# # -------------------------------------------------------------
# # Observed vs estimated convergence rate 
# # -------------------------------------------------------------

# 1. Load data exported from Julia
df_raw = pd.read_csv('../data/two_grid_convergence.csv')

# 2. Reshape from wide to long format
df_plot = df_raw.melt(
    id_vars=['k'],
    value_vars=['res_J', 'bound_J', 'res_wJ', 'bound_wJ', 'res_gs', 'bound_gs'],
    var_name='Series',
    value_name='Relative_Residual'
)

# Filter out NaNs to keep the plot clean after convergence
df_plot = df_plot.dropna(subset=['Relative_Residual'])

# 3. Map series names to Smoother type and Line Type (Observed vs Bound)
series_mapping = {
    'res_J': ('Jacobi (ω=1.0)', 'Observed'),
    'bound_J': ('Jacobi (ω=1.0)', 'XZ Bound'),
    'res_wJ': ('Weighted Jacobi (ω=2/3)', 'Observed'),
    'bound_wJ': ('Weighted Jacobi (ω=2/3)', 'XZ Bound'),
    'res_gs': ('Gauss-Seidel', 'Observed'),
    'bound_gs': ('Gauss-Seidel', 'XZ Bound'),
}

df_plot[['Smoother', 'Type']] = pd.DataFrame(
    df_plot['Series'].map(series_mapping).tolist(), 
    index=df_plot.index
)

# 4. Set drawing order and categorical factors
drawing_order = [
    'Jacobi (ω=1.0)',
    'Weighted Jacobi (ω=2/3)',
    'Gauss-Seidel'
]

df_plot['Smoother'] = pd.Categorical(
    df_plot['Smoother'], 
    categories=drawing_order, 
    ordered=True
)

df_plot['Type'] = pd.Categorical(
    df_plot['Type'], 
    categories=['Observed', 'XZ Bound'], 
    ordered=True
)

# Sort so observed curves render cleanly on top of theoretical bounds
df_plot = df_plot.sort_values(by=['Smoother', 'Type'], ascending=[True, False])

# 5. Define color scheme and linetypes matching your spectrum plot
colors = {
    'Jacobi (ω=1.0)': '#e76f51',          # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',  # Amber / Sandy
    'Gauss-Seidel': '#2a9d8f'              # Teal
}

linetypes = {
    'Observed': 'solid',
    'XZ Bound': 'dashed'
}

# 6. Generate plot
plot = (
    ggplot(df_plot, aes(x='k', y='Relative_Residual', color='Smoother', linetype='Type'))
    + geom_line(size=1.1, alpha=0.85)
    + scale_y_log10()
    + scale_color_manual(values=colors)
    + scale_linetype_manual(values=linetypes)
    + labs(
        x='Iteration (k)',
        y='Normalized Residual Norm ||r_k|| / ||r_0||',
        title='Two-Grid Convergence vs Theoretical XZ Bound',
        subtitle=r'Observed residual decay compared against $\sqrt{1 - \mu_{n_c+1}}^k$',
        color='Smoother',
        linetype='Line Type'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right'
    )
)

# Display or save plot
plot