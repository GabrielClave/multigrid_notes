using LinearAlgebra
using CSV
using DataFrames

n = 501
n_c = 250

# 1D Laplacian
A = SymTridiagonal(2.0 * ones(n), -1.0 * ones(n - 1))

# Exact Fourier basis Q (Euclidean orthonormal)
k_indices = 1:n
λ_A = [4.0 * sin(k * π / (2 * (n + 1)))^2 for k in k_indices]
Q = [sqrt(2 / (n + 1)) * sin(i * k * π / (n + 1)) for i in 1:n, k in 1:n]

# Additive Schwartz
function build_additive_schwarz(A, d::Int, overlap::Int)
    n = size(A, 1)
    @assert d >= 1 "Number of subdomains d must be at least 1."
    @assert d <= n "Number of subdomains d cannot exceed n."

    # Base size of each non-overlapping partition
    block_size = div(n, d)
    
    # Store subdomains as index ranges
    subdomains = Vector{UnitRange{Int}}(undef, d)
    
    for i in 1:d
        # Core non-overlapping bounds
        start_idx = (i - 1) * block_size + 1
        end_idx = (i == d) ? n : i * block_size  # last domain absorbs remainder
        
        # Expand bounds with overlap
        ov_start = max(1, start_idx - overlap)
        ov_end = min(n, end_idx + overlap)
        
        subdomains[i] = ov_start:ov_end
    end

    # Build R_as matrix
    R_as = zeros(eltype(A), n, n)
    A_mat = Matrix(A)

    for idx in subdomains
        # Subdomain solver: R_i^T * (A_i)^(-1) * R_i
        R_as[idx, idx] .+= inv(A_mat[idx, idx])
    end

    return R_as, subdomains
end

ndomain = 10
overlap = 10

R_as, _ = build_additive_schwarz(A, ndomain, overlap)

S_as = Matrix(I, n, n) - R_as * A_mat

R_bar = R_as' + R_as - R_as' * A * R_as
# R_bar should always be symmetric
issymmetric(R_bar)  #no
isapprox(R_bar', R_bar, atol=1e-11) # true
isapprox(0.5 * (R_bar + R_bar'), R_bar, atol=1e-11) # true

# probably the condition Number: 
cond(A)
cond(R_as)

norm( Symmetric(0.5 * (R_bar + R_bar') ) - R_bar, Inf) #e-14
norm( Symmetric(R_bar) - R_bar, Inf) #e-14

T_bar = R_bar*A
issymmetric(A*T_bar)  #no
isapprox(A*T_bar, (A*T_bar)', atol=1e-11) # true
norm( Symmetric(A*T_bar) - A*T_bar, Inf) #e-13

# -------------------------------------------------------------
# 1. MODAL COARSE SPACE (First n_c Fourier modes)
# -------------------------------------------------------------
P_modal = Q[:, 1:n_c] #range(P) = low frequency modes
Q_high = Q[:, (n_c + 1):n] # H = high frequency modes
Λ_high = Diagonal(λ_A[(n_c + 1):n])

Π_c_modal = P_modal * inv(P_modal' * A * P_modal) * P_modal' * A # A-orthogonal projection onto range(P)
X_modal_full = (I - Π_c_modal) * T_bar # n x n matrix here, acts on the full V, not just H

# Check A-self-adjointness on V: A*X == X'*A
println("X_modal is A-self-adjoint on V: ", isapprox(A * X_modal_full, X_modal_full' * A, atol=1e-11))
# not A-self-adjoint -> are we missing something or is it floating point cachotteries ? 

# Restricted operator on H^A using A-orthonormal basis Q_high_A
Q_high_A = Q_high * inv(sqrt(Λ_high))
X_modal_H = Q_high_A' * A * X_modal_full * Q_high_A # Euclidean symmetric -> TO DO
println("X_modal_H is symmetric: ", isapprox(X_modal_H, X_modal_H', atol=1e-11))

λ_X_modal = sort(real(eigen(X_modal_H).values), rev=true)

λ_X_modal = sort(real(eigen(X_modal_full).values), rev=true)

# -------------------------------------------------------------
# 2. GMG LINEAR INTERPOLATION COARSE SPACE
# -------------------------------------------------------------
P_gmg = zeros(n, n_c)
for j in 1:n_c
    i = 2 * j
    P_gmg[i - 1, j] = 0.5
    P_gmg[i,     j] = 1.0
    P_gmg[i + 1, j] = 0.5
end

Π_c_gmg = P_gmg * inv(P_gmg' * A * P_gmg) * P_gmg' * A
X_gmg_full = (I - Π_c_gmg) * T_bar

λ_X_gmg = sort(real(eigen(X_gmg_full).values), rev=true)

# -------------------------------------------------------------
# 3. XZ-OPTIMAL COARSE SPACE (Smallest eigenmodes of T_bar)
# -------------------------------------------------------------

R_bar = Symmetric(0.5 * (R_bar + R_bar')) # enforce symmetry
M = A * R_bar * A

issymmetric(M)  # no
isapprox(M, M', atol=1e-11) # true
norm( Symmetric(M) - M, Inf) #e-14

M = Symmetric(M)

# we want to solve R_bar * A * v = λv
# which is equivalent to A * R_bar * A * v = λA * v
F_Tbar = eigen(M, A)
# returns an A-orthonormal basis V, such that V^T * A * V = I

idx_opt = sortperm(F_Tbar.values) # Sort ascending μ_1 <= μ_2 ...
μ_Tbar = F_Tbar.values[idx_opt]
Q_Tbar = F_Tbar.vectors[:, idx_opt] # A-orthogonal

# in XZ they use a R_bar^-1 orthogonal basis
# we want to solve R_bar * A * v = λv
# which is equivalent to A * v = λR_bar^-1 * v
# so we could use F_Tbar = eigen(A, inv(R_bar))

rank(Q_Tbar) # n
println("||QT*A*Q - I||: ", norm(Q_Tbar'*A*Q_Tbar - I, Inf)) #e-13

P_opt = Q_Tbar[:, 1:n_c] # optimal range(P)

Π_c_opt = P_opt * inv(P_opt' * A * P_opt) * P_opt' * A
X_opt_full = (I - Π_c_opt) * T_bar

λ_X_opt = sort(real(eigen(X_opt_full).values), rev=true)
sum(λ_X_opt .> 1e-12) #235 > n_c

# -------------------------------------------------------------
# EXPORT TO CSV
# -------------------------------------------------------------
k_high_norm = [(n_c + k) / (n + 1) for k in 1:n]

df = DataFrame(
    k_mode = 1:n,
    k_norm = k_high_norm,
    Modal = λ_X_modal,
    GMG = λ_X_gmg,
    XZ_Optimal = λ_X_opt
)

CSV.write("experiments/data/spectrum_X_Schwartz.csv", df)