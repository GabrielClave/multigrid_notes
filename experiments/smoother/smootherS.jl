# computing the eigenvalues of the Smoother S = I - RA
# for Jacobi, Gauss Seidel, and Additive Schwartz in 1D Laplacian problem

using LinearAlgebra
using CSV
using DataFrames

n = 499
n_c = 249

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
S_jacobi_mat = Matrix(I, n, n) - 0.5 * Matrix(A)
λ_J = 1 .- 1/2 .* λ_A
# abs.(λ_J)

# Weighted Jacobi (ω = 2/3)
S_wjacobi_mat = Matrix(I, n, n) - (1/3) * Matrix(A)
λ_J = 1 .- 1/3 .* λ_A

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

function compute_mode_attenuation(S::AbstractMatrix, Q::AbstractMatrix)
    # Q columns are q_k; compute ||S * q_k||_2 for each column k
    return [norm(S * Q[:, k]) for k in axes(Q, 2)]
end

# Absolute values of eigenvalues, sorted in increasing order
df_sorted_abs_ev = DataFrame(
    k = 1:n,
    Jacobi = get_sorted_abs_eigenvalues(S_jacobi_mat),
    WeightedJacobi = get_sorted_abs_eigenvalues(S_wjacobi_mat),
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

CSV.write("experiments/data/smoother_spectrum.csv", df_sorted_abs_ev)
CSV.write("experiments/data/fourier_attenuation.csv", df_fourier_attenuation)