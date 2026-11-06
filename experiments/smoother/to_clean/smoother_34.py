import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_vline, geom_hline, element_rect,
    scale_color_manual, labs, theme_minimal, theme
)
import os

bg_color = "#FAFAFA"

# -------------------------------------------------------------
# Color palette and labels
# -------------------------------------------------------------

smoother_labels = {
    'Jacobi': 'Jacobi (ω=1.0)',
    'Jacobi_true': 'Jacobi (predicted)',
    'WeightedJacobi': 'Weighted Jacobi (ω=2/3)',
    'WeightedJacobi_23': 'Weighted Jacobi (ω=2/3)',
    'WeightedJacobi_true': 'Weighted Jacobi (ω=2/3, predicted)',
    'WeightedJacobi_06': 'Weighted Jacobi (ω=0.601)',
    'WeightedJacobi_06_attenuation': 'Weighted Jacobi (ω=0.601)',
    'GaussSeidel': 'Gauss-Seidel',
    'AdditiveSchwarz': 'Additive Schwarz',
    'Jacobi_attenuation': 'Jacobi (ω=1.0)',
    'WeightedJacobi_attenuation': 'Weighted Jacobi (ω=2/3)',
    'GaussSeidel_attenuation': 'Gauss-Seidel',
    'AdditiveSchwarz_attenuation': 'Additive Schwarz'
}

colors = {
    'Jacobi (ω=1.0)': '#e76f51',                      # Terracotta
    'Weighted Jacobi (ω=2/3)': '#f4a261',               # Amber / Sandy
    'Weighted Jacobi (ω=0.601)': '#e9c46a',             # Warm Yellow (3:4 optimal)
    'Gauss-Seidel': '#2a9d8f',                         # Teal
    'Additive Schwarz': '#264653',                       # Dark slate
    'Jacobi (predicted)': '#99A2AC',
    'Weighted Jacobi (ω=2/3, predicted)': '#C298BC'
}

# -------------------------------------------------------------
# 1. Eigenvalues of S = I - RA
# -------------------------------------------------------------

df_spec_raw = pd.read_csv('../data/df_34_smoother_spectrum.csv')

if 'normalized_index' not in df_spec_raw.columns:
    df_spec_raw['normalized_index'] = df_spec_raw['k'] / df_spec_raw['k'].max()

df_spec_plot = df_spec_raw.melt(
    id_vars=['k', 'normalized_index'],
    value_vars=[
        'Jacobi', 'Jacobi_true', 
        'WeightedJacobi', 'WeightedJacobi_true', 
        'WeightedJacobi_06', 
        'GaussSeidel', 'AdditiveSchwarz'
    ],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
).dropna()

df_spec_plot['Smoother'] = df_spec_plot['Smoother'].map(smoother_labels)
n = len(df_spec_raw)

plot_spec = (
    ggplot(df_spec_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(S)|',
        title='Spectrum of S = I - RA (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        plot_background=element_rect(fill=bg_color, color=bg_color)
    )
)

# -------------------------------------------------------------
# 2. Attenuation factor ||S q_k||_2 across Fourier Spectrum
# -------------------------------------------------------------

df_raw = pd.read_csv('../data/df_34_fourier_attenuation.csv')

df_plot = df_raw.melt(
    id_vars=['k', 'normalized_frequency'],
    value_vars=[
        'Jacobi_attenuation', 
        'WeightedJacobi_attenuation', 
        'WeightedJacobi_06_attenuation',
        'GaussSeidel_attenuation', 
        'AdditiveSchwarz_attenuation'
    ],
    var_name='Smoother',
    value_name='Attenuation'
).dropna()

df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)
n = df_raw.shape[0]

# Updated threshold for 3:4 coarsening (3m/8 = 0.375)
cutoff_34 = 0.375

plot_fourier = (
    ggplot(df_plot, aes(x='normalized_frequency', y='Attenuation', color='Smoother'))
    + geom_line(size=0.7, alpha=0.85)
    + geom_vline(xintercept=cutoff_34, linetype='dashed', color='#8d99ae', size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Frequency Mode',
        y='Attenuation Factor ||S q_k||₂',
        title='Smoother Mode Attenuation across Fourier Spectrum',
        subtitle=f'1D Laplacian | n = {n} (Vertical line at 3:4 cutoff = 0.375)'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        plot_background=element_rect(fill=bg_color, color=bg_color)
    )
)

# -------------------------------------------------------------
# 3. Eigenvalues of T_bar = R_bar * A
# -------------------------------------------------------------

df_raw = pd.read_csv('../data/df_34_smoother_T_spectrum.csv')

if 'normalized_index' not in df_raw.columns:
    df_raw['normalized_index'] = df_raw['k'] / df_raw['k'].max()

df_plot = df_raw.melt(
    id_vars=['k', 'normalized_index'],
    value_vars=['Jacobi', 'WeightedJacobi', 'WeightedJacobi_06', 'GaussSeidel', 'AdditiveSchwarz'],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
).dropna()

df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)
n = len(df_raw)

plot_T = (
    ggplot(df_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude λ(T_bar)',
        title='Spectrum of T_bar = R_bar * A (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        plot_background=element_rect(fill=bg_color, color=bg_color)
    )
)

# -------------------------------------------------------------
# 4. Eigenvalues of X = (I - Pi_c) (I - R_bar A) for 3:4 transfer (P_34)
# -------------------------------------------------------------

df_raw = pd.read_csv('../data/df_34_smoother_X_spectrum.csv')

if 'normalized_index' not in df_raw.columns:
    df_raw['normalized_index'] = df_raw['k'] / df_raw['k'].max()

df_plot = df_raw.melt(
    id_vars=['k', 'normalized_index'],
    value_vars=['Jacobi', 'WeightedJacobi_23', 'WeightedJacobi_06', 'GaussSeidel', 'AdditiveSchwarz'],
    var_name='Smoother',
    value_name='Eigenvalue_Magnitude'
).dropna()

df_plot['Smoother'] = df_plot['Smoother'].map(smoother_labels)
n = len(df_raw)

plot_X_34 = (
    ggplot(df_plot, aes(x='normalized_index', y='Eigenvalue_Magnitude', color='Smoother'))
    + geom_line(size=1.1, alpha=0.85)
    + geom_hline(yintercept=1.0, linetype='dashed', color="#4e4774", size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Eigenvalue Index',
        y='Eigenvalue Magnitude |λ(X)|',
        title='Spectrum of X with 3:4 Interpolation P_3/4 (Sorted Eigenvalue Magnitudes)',
        subtitle=f'1D Laplacian | n = {n}'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
        plot_background=element_rect(fill=bg_color, color=bg_color)
    )
)