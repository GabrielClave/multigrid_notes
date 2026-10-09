# computing the eigenvalues of the symmetrized smoother T_bar = R_bar A
# for Jacobi, Gauss Seidel, and Additive Schwartz in 1D Laplacian problem

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

# Jacobi

# R = 1/2*I
# S = I - 1/2*A
# Jacobi (unweighted, ω = 1)
R_jacobi_mat = 0.5 * Matrix(I, n, n)
S_jacobi_mat = Matrix(I, n, n) - 0.5 * Matrix(A)
λ_S = 1 .- 1/2 .* λ_A
λ_S = sort(abs.(λ_S))
# R_bar*A = 2/3*A(2I - 2/3A) 
λ_X = 1/2 * ( 2 .- 1/2 .* λ_A ) .* λ_A
λ_X = sort(λ_X)
# λ_X = sort(abs.(λ_X))

# Weighted Jacobi (ω = 2/3)
R_wjacobi_mat = (1.0/3) * Matrix(I, n, n)
S_wjacobi_mat = Matrix(I, n, n) - (1.0/3) * A
λ_Sw = 1 .- 1/3 .* λ_A
λ_Sw = sort(abs.(λ_Sw))
# R_bar*A = 1/3*A(2I - 1/3A) 
λ_Xw = 1/3 * ( 2 .- 1/3 .* λ_A ) .* λ_A
λ_Xw = sort(λ_Xw)
# λ_Xw = sort(abs.(λ_Xw))

# Gauss-Seidel
D_L = LowerTriangular(A)
R_gs = inv(D_L)
S_gs = I - R_gs*A

# Additive Schwarz
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
S_as = I - R_as*A

# functions

function get_sorted_abs_eigenvalues(S::AbstractMatrix)
    evs = eigen(S).values
    return sort(abs.(evs))
end

function compute_mode_attenuation(S::AbstractMatrix, Q::AbstractMatrix, anorm = false)
    if anorm
        return [sqrt(dot(S*Q[:, k], A, S*Q[:, k])/dot(Q[:, k],A,Q[:, k])) for k in axes(Q, 2)] # ||S * q_k||_A for each column k
    end
    # Q columns are q_k; compute ||S * q_k||_2 for each column k
    return [norm(S * Q[:, k]) for k in axes(Q, 2)]
end

function get_sorted_S_eigenvalues(R::AbstractMatrix)

    # R_bar = Symmetric(R' + R - R' * A * R)
    # T_bar = Symmetric(R_bar*A)
    R_bar = R' + R - R' * A * R
    T_bar = R_bar*A
    # evs = eigen(I - T_bar).values
    evs = real.(eigen(I - T_bar).values)
    return sort(abs.(evs))
end

Λ_high = Diagonal(λ_A[(n_c + 1):n])
Q_high = Q[:, (n_c + 1):n]
Q_high_A = Q_high * inv(sqrt(Λ_high))

P_gmg = zeros(n, n_c)
for j in 1:n_c
    i = 2 * j
    P_gmg[i - 1, j] = 0.5
    P_gmg[i,     j] = 1.0
    P_gmg[i + 1, j] = 0.5
end

function get_sorted_X_eigenvalues(R::AbstractMatrix, P::AbstractMatrix)

    # R_bar = Symmetric(R' + R - R' * A * R)
    # T_bar = Symmetric(R_bar*A)
    R_bar = R' + R - R' * A * R
    T_bar = R_bar*A
    Π_c = P * inv(P' * A * P) * P' * A # A-orthogonal projection onto range(P)
    X_full = (I - Π_c) * T_bar # n x n matrix here, acts on the full V, not just H
    # X_H = Symmetric(Q_high_A' * A * X_full * Q_high_A) # Euclidean symmetric -> is that what we want ?

    # evs = eigen(X_H).values
    # return sort(evs)
    evs = real(eigen(X_full).values) #no guarantee to be real
end

function get_sorted_X_eigenvalues(R::AbstractMatrix)

    # R_bar = Symmetric(R' + R - R' * A * R)
    R_bar = R' + R - R' * A * R
    # we want to solve R_bar * A * v = λv
    # which is equivalent to A * R_bar * A * v = λA * v
    M = Symmetric(A * R_bar * A)
    F_Tbar = eigen(M, A)

    idx_opt = sortperm(F_Tbar.values) # Sort ascending μ_1 <= μ_2 ...
    Q_Tbar = F_Tbar.vectors[:, idx_opt] # A-orthogonal
    P = Q_Tbar[:, 1:n_c]

    Π_c = P * inv(P' * A * P) * P' * A # A-orthogonal projection onto range(P)
    X_full = (I - Π_c) * R_bar*A # n x n matrix here, acts on the full V, not just H
    # X_H = Symmetric(Q_high_A' * A * X_full * Q_high_A) # Euclidean symmetric -> is that what we want ?

    # evs = eigen(X_H).values
    # return sort(evs)
    evs = real(eigen(X_full).values) #no guarantee to be real

end

# # -------------------------------------------------------------
# experiments

# Absolute values of eigenvalues, sorted in increasing order
df_sorted_abs_ev = DataFrame(
    k = 1:n,
    Jacobi = get_sorted_abs_eigenvalues(S_jacobi_mat),
    Jacobi_true = λ_S,
    WeightedJacobi = get_sorted_abs_eigenvalues(S_wjacobi_mat),
    WeightedJacobi_true = λ_Sw,
    GaussSeidel = get_sorted_abs_eigenvalues(S_gs),
    AdditiveSchwarz = get_sorted_abs_eigenvalues(S_as)
)

# Unsorted Fourier mode amplification factors a(k) = ||S q_k||_2
df_fourier_attenuation = DataFrame(
    k = 1:n,
    normalized_frequency = (1:n) ./ (n + 1),
    Jacobi_attenuation = compute_mode_attenuation(S_jacobi_mat, Q),
    WeightedJacobi_attenuation = compute_mode_attenuation(S_wjacobi_mat, Q),
    GaussSeidel_attenuation = compute_mode_attenuation(S_gs, Q),
    AdditiveSchwarz_attenuation = compute_mode_attenuation(S_as, Q)
)

# Unsorted Fourier mode amplification factors a(k) = ||S q_k||_A in A-norm
df_fourier_attenuation_anorm = DataFrame(
    k = 1:n,
    normalized_frequency = (1:n) ./ (n + 1),
    Jacobi_attenuation = compute_mode_attenuation(S_jacobi_mat, Q, true),
    WeightedJacobi_attenuation = compute_mode_attenuation(S_wjacobi_mat, Q, true),
    GaussSeidel_attenuation = compute_mode_attenuation(S_gs, Q, true),
    AdditiveSchwarz_attenuation = compute_mode_attenuation(S_as, Q, true)
)

# abs eigenvalues of S_bar, sorted in increasing order
df_sorted_S_ev = DataFrame(
    k = 1:n,
    Jacobi = get_sorted_S_eigenvalues(R_jacobi_mat),
    WeightedJacobi = get_sorted_S_eigenvalues(R_wjacobi_mat),
    GaussSeidel = get_sorted_S_eigenvalues(R_gs),
    AdditiveSchwarz = get_sorted_S_eigenvalues(R_as)
)

# eigenvalues of X_H, sorted in increasing order
df_sorted_X_ev = DataFrame(
    # k = 1:(n-n_c),
    k = 1:n,
    Jacobi = get_sorted_X_eigenvalues(R_jacobi_mat, P_gmg),
    WeightedJacobi = get_sorted_X_eigenvalues(R_wjacobi_mat, P_gmg),
    GaussSeidel = get_sorted_X_eigenvalues(R_gs, P_gmg),
    AdditiveSchwarz = get_sorted_X_eigenvalues(R_as, P_gmg)
)

# eigenvalues of optimal X_H, sorted in increasing order
df_sorted_X_opt_ev = DataFrame(
    # k = 1:(n-n_c),
    k = 1:n,
    Jacobi = get_sorted_X_eigenvalues(R_jacobi_mat),
    Jacobi_true = λ_X,
    WeightedJacobi = get_sorted_X_eigenvalues(R_wjacobi_mat),
    WeightedJacobi_true = λ_Xw,
    GaussSeidel = get_sorted_X_eigenvalues(R_gs),
    AdditiveSchwarz = get_sorted_X_eigenvalues(R_as)
)

# # -------------------------------------------------------------
# save data

CSV.write("experiments/data/smoother_spectrum.csv", df_sorted_abs_ev)
CSV.write("experiments/data/fourier_attenuation.csv", df_fourier_attenuation)
CSV.write("experiments/data/fourier_attenuation_anorm.csv", df_fourier_attenuation_anorm)
CSV.write("experiments/data/smoother_S_spectrum.csv", df_sorted_S_ev)
CSV.write("experiments/data/smoother_X_spectrum.csv", df_sorted_X_ev)
CSV.write("experiments/data/smoother_X_opt_spectrum.csv", df_sorted_X_opt_ev)

# # -------------------------------------------------------------
# two-grid convergence


function two_grid_step!(A, f, R, A_c_inv, P, r, v)
    # coarse grid correction

    # r .= f .- A*v # assumed to have been computed
    # r_c = P'*r
    # e_c = A_c_inv*e_c
    # v = v + P * e_c
    v .+= P* A_c_inv * P' * r

    # post-smooth one time
    r .= f .- A*v

    v .+= R*r    
end

function two_grid_cycle(A, f, R, P, r, v, residuals; maxiter = 100, tol = 1e-9)

    A_c_inv = inv(P' * A * P)

    k = 1
    r = f .- A*v

    nr = norm(r)
    residuals[k] = nr

    while k < maxiter

        if nr < tol
            print("converged in $k steps")
            return k, residuals
        end

        two_grid_step!(A, f, R, A_c_inv, P, r, v)
        k += 1
        r .= f .- A * v
        nr = norm(r)
        residuals[k] = nr
    end

    print("did not converge")

    return maxiter, residuals
end

maxiter = 25
f = rand(n)
r = zeros(n)
v = zeros(n)

k_J, res_hist_J = two_grid_cycle(A, f, R_jacobi_mat, P_gmg, r, copy(v), fill(NaN, maxiter); maxiter)
k_wJ, res_hist_wJ = two_grid_cycle(A, f, R_wjacobi_mat, P_gmg, r, copy(v), fill(NaN, maxiter); maxiter)
k_gs, res_hist_gs = two_grid_cycle(A, f, R_gs, P_gmg, r, copy(v), fill(NaN, maxiter); maxiter)
# k_as, res_hist_as = two_grid_cycle(A, f, R_as, P_gmg, r, copy(v), fill(NaN, maxiter); maxiter)

# theoretical XZ bound μ_nc+1
tol_zero = 1e-12
mu_J  = sort(filter(x -> abs(x) > tol_zero, df_sorted_X_ev.Jacobi))[1]
mu_wJ = sort(filter(x -> abs(x) > tol_zero, df_sorted_X_ev.WeightedJacobi))[1]
mu_gs = sort(filter(x -> abs(x) > tol_zero, df_sorted_X_ev.GaussSeidel))[1]

cvr_J  = sqrt(1 - mu_J)
cvr_wJ = sqrt(1 - mu_wJ)
cvr_gs = sqrt(1 - mu_gs)

bound_J  = [cvr_J^k  for k in 0:maxiter-1]
bound_wJ = [cvr_wJ^k for k in 0:maxiter-1]
bound_gs = [cvr_gs^k for k in 0:maxiter-1]

# save dataframe
df_convergence = DataFrame(
    k = 1:maxiter,
    res_J = res_hist_J ./ res_hist_J[1],
    bound_J = bound_J,
    res_wJ = res_hist_wJ ./ res_hist_wJ[1],
    bound_wJ = bound_wJ,
    res_gs = res_hist_gs ./ res_hist_gs[1],
    bound_gs = bound_gs,
    # res_as = res_hist_as ./ res_hist_as[1],
)

CSV.write("experiments/data/two_grid_convergence.csv", df_convergence)
