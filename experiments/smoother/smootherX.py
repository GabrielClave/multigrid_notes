import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_vline, 
    scale_color_manual, labs, theme_minimal, theme
)

# Parameters
n = 100                  # Grid size / total modes
n_c = 40                 # Number of coarse modes kept
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

# Plot using plotnine
plot = (
    ggplot(df, aes(x='k_norm', y='Eigenvalue', color='Operator'))
    + geom_line(size=1.2)
    + geom_vline(xintercept=k_c, linetype='dashed', color='#8d99ae', size=0.8)
    + scale_color_manual(values=colors)
    + labs(
        x='Normalized Frequency Mode (k / (n + 1))',
        y='Eigenvalue Magnitude',
        title='Spectrum of T_bar, and Two-Grid Operator X',
        subtitle=f'Coarse space dimension n_c = {n_c}, Jacobi weight omega = 2/3'
    )
    + theme_minimal()
    + theme(
        figure_size=(9, 5),
        legend_position='right',
    )
)

# Display / save plot
plot