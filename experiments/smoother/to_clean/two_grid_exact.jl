using LinearAlgebra
using Random
using DataFrames
using CSV

Random.seed!(1111)

# ==============================================================================
# "ideal" smoother
# ==============================================================================

n = 20
k = 10

# Generate a random orthogonal matrix Q
M = randn(n, n)
F = qr(M)
# Normalize signs of R's diagonal to ensure strict uniform Haar distribution (whatever pleases you)
d = sign.(diag(F.R))
d[d .== 0] .= 1.0
Q = F.Q * Diagonal(d)

# Generate SPD matrix A with spectrum spaced from 10 to 1
λ = range(10, 1, length=n)
D = Diagonal(λ)
A = Q * D * Q'

A = Symmetric(A)

# Generate injective P such that range(P) = range(Q[:, k+1:n])
dim_subspace = n - k

# Non-orthogonal invertible matrix B
C = randn(dim_subspace, dim_subspace) 

P = Q[:, k+1:n] * C

λ_B = zeros(n)
λ_B[1:k] = 1 ./ λ[1:k]

B = Q * Diagonal(λ_B) * Q'
# ideal smoother: perfectly invert A on the first k eigenvectors

# Smoother error propagation matrix
E_s = I - B * A

# 5. Coarse grid correction operators
A_c = P' * A * P
# Correction update map: v <- v + P * (A_c \ (P' * r))
# Error operator: E_c = I - P * (A_c \ (P' * A))
E_c = I - P * (A_c \ (P' * A))

# 6. Setup problem: Au = f
f = range(10, 0, length=n)
u = A \ f
v0 = zeros(n)

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

e_cgc_only = E_c * e0      # After CGC alone

# ==============================================================================
# DATASET 1: Vector Components (Euclidean vs Eigenbasis)
# ==============================================================================

# Calculate basis components
e0_hat = Q' * e0
e1_hat = Q' * e1
e2_hat = Q' * e2

# total series: 3 stages x 2 bases
df_components = DataFrame(
    Index = repeat(1:n, 6),
    Component = vcat(
        # Euclidean basis
        e0, e1, e2,
        # Eigenvector basis Q
        e0_hat, e1_hat, e2_hat
    ),
    Stage = vcat(
        fill("1. Initial (e0)", n),
        fill("2. After Smoother (e1)", n),
        fill("3. After CGC (e2)", n),
        fill("1. Initial (e0)", n),
        fill("2. After Smoother (e1)", n),
        fill("3. After CGC (e2)", n)
    ),
    Basis = vcat(
        fill("Euclidean", 3 * n),
        fill("Eigenvector", 3 * n)
    )
)

# ==============================================================================
# DATASET 2: 2D Subspace Norms & Trajectories
# ==============================================================================

function get_subspace_coords(e, Q, k)
    e_hat = Q' * e
    return norm(e_hat[1:k]), norm(e_hat[k+1:end])
end

x0, y0 = get_subspace_coords(e0, Q, k)
x1, y1 = get_subspace_coords(e1, Q, k)
xc, yc = get_subspace_coords(e_cgc_only, Q, k)
x2, y2 = get_subspace_coords(e2, Q, k)

df_2d = DataFrame(
    x = [x0, x1, xc, x2],
    y = [y0, y1, yc, y2],
    Stage = ["e0 (Initial)", "e1 (After Smoother)", "e_cgc (After CGC)", "e2 (Full Cycle)"]
)

df_lines = DataFrame(
    x = [x0, x1, x0],
    y = [y0, y1, y0],
    xend = [x1, x2, xc],
    yend = [y1, y2, yc]
)

# ==============================================================================
# EXPORT TO CSV
# ==============================================================================

CSV.write("experiments/data/df_components.csv", df_components)
CSV.write("experiments/data/df_2d.csv", df_2d)
CSV.write("experiments/data/df_lines.csv", df_lines)