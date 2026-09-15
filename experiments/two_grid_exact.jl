using LinearAlgebra
using Random
using DataFrames
using AlgebraOfGraphics
using CairoMakie

Random.seed!(42)

n = 20
k = 10

# 1. Generate a random orthogonal matrix Q
M = randn(n, n)
F = qr(M)
# Normalize signs of R's diagonal to ensure strict uniform Haar distribution
d = sign.(diag(F.R))
d[d .== 0] .= 1.0
Q = F.Q * Diagonal(d)

# 2. Generate SPD matrix A with spectrum log-spaced from 10^-6 to 10^3
# (Using log-spacing ensures an even distribution of orders of magnitude)
λ = range(100, 0.1, length=n)
D = Diagonal(λ)
A = Q * D * Q'

# Ensure exact symmetry numerically (mitigates float precision noise)
A = Symmetric(A)

# 3. Generate injective P such that range(P) = range(Q[:, k+1:n])
dim_subspace = n - k

# Non-orthogonal invertible matrix B
C = randn(dim_subspace, dim_subspace) 

P = Q[:, k+1:n] * C

λ_B = zeros(n)
λ_B[1:k] = 1 ./ λ[1:k]

# B = Q * Diagonal(λ_B) * Q'
B = Q * Diagonal(λ_B) * Q'

# Smoother error propagation matrix
E_s = I - B * A

# 5. Coarse grid correction operators
A_c = P' * A * P
# Correction update map: v <- v + P * (A_c \ (P' * r))
# Error operator: E_c = I - P * (A_c \ (P' * A))
E_c = I - P * (A_c \ (P' * A))

# 6. Setup problem: Au = f
f = ones(n)
u = A \ f
v0 = randn(n)

e0 = u - v0
r0 = f - A * v0

# --- STEP 1: Apply the Smoother ---
v1 = v0 + B * r0
e1 = u - v1

# in the Q basis we should see the first k components vanish:
e_L = Q'*e1

# and the opposite for the coarse grid correction
e_H = Q'*E_c*e0

# Verify theoretical smoother error propagation: e1 == (I - BA)e0
println("Smoother error residual norm: ", norm(e1 - E_s * e0)) # e-16

# --- STEP 2: Apply the Coarse Grid Correction ---
r1 = f - A * v1
v2 = v1 + P * (A_c \ (P' * r1))
e2 = u - v2

# Verify theoretical combined error propagation: e2 == E_c * e1
println("CGC error residual norm vs theoretical: ", norm(e2 - E_c * e1)) # e-11:

# --- VERIFICATION ---
println("Residual norm after 1 full 2-grid cycle ||f - Av2||: ", norm(f - A * v2)) #e-12
println("Error norm after 1 full 2-grid cycle ||e2||: ", norm(e2)) #e-12

# States
e1 = E_s * e0              # After Smoother
e_cgc_only = E_c * e0      # After CGC alone
e2 = E_c * E_s * e0        # After full 2-grid cycle

# ==============================================================================
# VISUALIZATION 1: Eigenbasis Component Bar Plots
# ==============================================================================

# Transform errors into the eigenbasis Q
e0_hat = Q' * e0
e1_hat = Q' * e1
e_cgc_hat = Q' * e_cgc_only

df_components = DataFrame(
    Index = repeat(1:n, 3),
    Component = vcat(e0_hat, e1_hat, e_cgc_hat),
    Stage = vcat(
        fill("1. Initial (e0)", n),
        fill("2. After Smoother (e1)", n),
        fill("3. After CGC alone", n)
    )
)

plot_components = ggplot(df_components, aes(x = :Index, y = :Component, fill = :Stage)) +
    geom_col(position = "dodge") +
    geom_vline(xintercept = k + 0.5, linetype = "dashed") +
    labs(
        title = "Error Spectrum Decomposed in Eigenbasis Q",
        x = "Eigenvector Index (1:k = High Freq, k+1:n = Low Freq)",
        y = "Component Value"
    )

display(plot_components)

# ==============================================================================
# VISUALIZATION 2: 2D Subspace Norm Projection Space
# ==============================================================================

function get_subspace_coords(e, Q, k)
    e_hat = Q' * e
    norm_H = norm(e_hat[1:k])
    norm_L = norm(e_hat[k+1:end])
    return norm_H, norm_L
end

x0, y0 = get_subspace_coords(e0, Q, k)
x1, y1 = get_subspace_coords(e1, Q, k)
xc, yc = get_subspace_coords(e_cgc_only, Q, k)
x2, y2 = get_subspace_coords(e2, Q, k)

df_2d = DataFrame(
    x = [x0, x1, xc, x2],
    y = [y0, y1, yc, y2],
    Label = ["e0 (Initial)", "e1 (After Smoother)", "e_cgc (After CGC)", "e2 (Full Cycle)"]
)

plot_subspace = ggplot(df_2d, aes(x = :x, y = :y, color = :Label)) +
    geom_point(size = 4) +
    geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), color = "gray", linetype = "dotted") +
    geom_segment(aes(x = x0, y = y0, xend = xc, yend = yc), color = "gray", linetype = "dotted") +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray", linetype = "dotted") +
    labs(
        title = "Error Trajectory in High vs Low Frequency Norm Space",
        x = "High-Frequency Subspace Norm ||P_H e||",
        y = "Low-Frequency Subspace Norm ||P_L e||"
    )

display(plot_subspace)