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

overlap = 10
mid = div(n, 2)
idx1 = 1:(mid + overlap)
idx2 = (mid - overlap + 1):n

R_as = zeros(n, n)
A_mat = Matrix(A)

# Additive Schwarz preconditioner: R_AS = R_1^T (A_1)^(-1) R_1 + R_2^T (A_2)^(-1) R_2
R_as[idx1, idx1] .+= inv(A_mat[idx1, idx1])
R_as[idx2, idx2] .+= inv(A_mat[idx2, idx2])

S_as = Matrix(I, n, n) - R_as * A_mat

R_bar = R_as' + R_as - R_as' * A * R_as
issymmetric(R_bar)  #no
isapprox(R_bar', R_bar, atol=1e-11) #I said no

T_bar = R_bar*A
issymmetric(A*T_bar)  #no
isapprox(A*T_bar, (A*T_bar)', atol=1e-11) #I said no

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

# -------------------------------------------------------------
# 2. XZ-OPTIMAL COARSE SPACE (Smallest eigenmodes of T_bar)
# -------------------------------------------------------------
F_Tbar = eigen(T_bar)
idx_opt = sortperm(real(F_Tbar.values)) # Sort ascending μ_1 <= μ_2 ...
μ_Tbar = F_Tbar.values[idx_opt]
Q_Tbar = F_Tbar.vectors[:, idx_opt] # is this one orthogonal ? Ri-orthogonal ? does it matter ?

rank(Q_Tbar)

P_opt = Q_Tbar[:, 1:n_c]

# Find exact zero columns
zero_cols = findall(j -> count(!iszero, view(P_opt, :, j)) == 0, axes(P_opt,2))

# Find columns with negligible norm
col_norms = [norm(view(P_opt, :, j)) for j in axes(P_opt,2)]
near_zero_cols = findall(col_norms .< 1e-12)

println("Zero columns: ", zero_cols)
println("Near-zero columns: ", near_zero_cols)

# Numerical rank using SVD tolerance
r = rank(Matrix(P_opt))
println("Size of P: $(size(P_opt)), Numerical Rank: $r")

S = svdvals(Matrix(P_opt))
println("Smallest 5 singular values of P: ", S[end-4:end])

Π_c_opt = P_opt * inv(P_opt' * A * P_opt) * P_opt' * A
X_opt_full = (I - Π_c_opt) * T_bar

λ_X_opt_nonzero = sort(filter(x -> x > 1e-8, real(eigen(X_opt_full).values)), rev=true)

# -------------------------------------------------------------
# 3. GMG LINEAR INTERPOLATION COARSE SPACE
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

λ_X_gmg_nonzero = sort(filter(x -> x > 1e-8, real(eigen(X_gmg_full).values)), rev=true)

# -------------------------------------------------------------
# EXPORT TO CSV
# -------------------------------------------------------------
num_high_modes = n - n_c
k_high_norm = [(n_c + k) / (n + 1) for k in 1:num_high_modes]

df = DataFrame(
    k_mode = 1:num_high_modes,
    k_norm = k_high_norm,
    Modal = λ_X_modal,
    GMG = λ_X_gmg_nonzero,
    XZ_Optimal = λ_X_opt_nonzero
)

CSV.write("experiments/data/spectrum_X_GS.csv", df)