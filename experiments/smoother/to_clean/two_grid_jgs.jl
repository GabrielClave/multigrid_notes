using LinearAlgebra
using SparseArrays
using Random
using DataFrames
using CSV

Random.seed!(42)

# ==============================================================================
# 1. GENERATE 2D LAPLACIAN MATRIX
# ==============================================================================
nx, ny = 50, 50               # Grid dimensions
n = nx * ny                   # Total degrees of freedom
k = div(n, 2)                 # Cutoff: first k are high frequencies

# 1D 3-point stencil Laplacian
function laplacian_1d(N)
    T = tridiagonal(N)
    return T
end

function tridiagonal(N)
    d = fill(2.0, N)
    dl = fill(-1.0, N-1)
    return SpDiagMat(dl, d, dl, N)
end

function SpDiagMat(dl, d, du, N)
    return SparseMatrixCSC(SymTridiagonal(d, dl))
end

# 2D Laplacian via Kronecker sum: A = I_ny ⊗ T_x + T_y ⊗ I_nx
Tx = Array(tridiagonal(nx))
Ty = Array(tridiagonal(ny))
A_sparse = kron(I(ny), Tx) + kron(Ty, I(nx))
A = Symmetric(Matrix(A_sparse))

# Spectral decomposition (sort descending: λ_1 >= λ_2 >= ... >= λ_n)
F = eigen(A)
idx = sortperm(F.values, rev=true)
λ = F.values[idx]
Q = F.vectors[:, idx]

# ==============================================================================
# 2. DEFINE SMOOTHERS & COARSE GRID CORRECTION
# ==============================================================================
D = Diagonal(diag(A))
L = LowerTriangular(A)

# Damped Jacobi smoother (ω = 2/3)
ω = 2/3
B_jacobi = ω * inv(D)

# Gauss-Seidel smoother
B_gs = inv(Matrix(D + L))

# Coarse grid space (spanned by low frequencies q_{k+1} ... q_n)
dim_subspace = n - k
C = randn(dim_subspace, dim_subspace)
P = Q[:, k+1:n] * C
A_c = P' * A * P
E_c = I - P * (A_c \ (P' * A))

# ==============================================================================
# 3. COMPUTE ERROR TRAJECTORIES
# ==============================================================================
function get_subspace_coords(e, Q, k)
    e_hat = Q' * e
    return norm(e_hat[1:k]), norm(e_hat[k+1:end])
end

# Initial random error e0
e0 = randn(n)
x0, y0 = get_subspace_coords(e0, Q, k)

results_2d = DataFrame()
results_lines = DataFrame()

smoothers = [("Damped Jacobi (ω=2/3)", B_jacobi), ("Gauss-Seidel", B_gs)]
for (name, B) in smoothers
    E_s = I - B * A
    
    # --- CYCLE 1 ---
    e1 = E_s * e0
    e2 = E_c * e1

    # --- CYCLE 2 ---
    e3 = E_s * e2
    e4 = E_c * e3

    x1, y1 = get_subspace_coords(e1, Q, k)
    x2, y2 = get_subspace_coords(e2, Q, k)
    x3, y3 = get_subspace_coords(e3, Q, k)
    x4, y4 = get_subspace_coords(e4, Q, k)

    # 2D Coordinates (e0 -> e1 -> e2 -> e3 -> e4)
    df_pts = DataFrame(
        x = [x0, x1, x2, x3, x4],
        y = [y0, y1, y2, y3, y4],
        State = [
            "e0 (Initial)", 
            "e1 (After Smoother 1)", 
            "e2 (After CGC 1)", 
            "e3 (After Smoother 2)", 
            "e4 (After CGC 2)"
        ],
        Smoother = fill(name, 5)
    )
    append!(results_2d, df_pts)

    # Sequential Trajectory segments
    df_lns = DataFrame(
        x = [x0, x1, x2, x3],
        y = [y0, y1, y2, y3],
        xend = [x1, x2, x3, x4],
        yend = [y1, y2, y3, y4],
        Smoother = fill(name, 4)
    )
    append!(results_lines, df_lns)
end

# Export exclusively the 2D plot data
CSV.write("experiments/data/df_2d_real.csv", results_2d)
CSV.write("experiments/data/df_lines_real.csv", results_lines)