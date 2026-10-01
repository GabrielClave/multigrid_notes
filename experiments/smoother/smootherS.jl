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

# Additive Schwarz (2 overlapping blocks, domain split in two with overlap)
# Define subdomains covering 1:n_sub1 and n_sub2:n
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