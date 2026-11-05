using LinearAlgebra
using CSV
using DataFrames

# Load two-grid abstraction module
include("twoGrid.jl")

# -----------------------------------------------------------------------------
# 1. Domain & Operator Construction Helpers
# -----------------------------------------------------------------------------

function build_1d_laplacian(n::Int)
    return SymTridiagonal(2.0 * ones(n), -1.0 * ones(n - 1))
end

function build_fourier_basis(n::Int)
    k_indices = 1:n
    λ_A = [4.0 * sin(k * π / (2 * (n + 1)))^2 for k in k_indices]
    Q = [sqrt(2 / (n + 1)) * sin(i * k * π / (n + 1)) for i in 1:n, k in 1:n]
    return λ_A, Q
end

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

function build_P_gmg(n::Int, n_c::Int)
    P = zeros(n, n_c)
    for j in 1:n_c
        i = 2 * j
        P[i - 1, j] = 0.5
        P[i,     j] = 1.0
        P[i + 1, j] = 0.5
    end
    return P
end

function build_P_34(n::Int, n_c_34::Int)
    P = zeros(n, n_c_34)
    for i in 1:n
        j_real = (3.0 * i) / 4.0
        j_floor = floor(Int, j_real)
        t = j_real - j_floor

        if 1 <= j_floor <= n_c_34
            P[i, j_floor] += 1.0 - t
        end
        if 1 <= j_floor + 1 <= n_c_34
            P[i, j_floor + 1] += t
        end
    end
    return P
end

function compute_fourier_attenuation(tgs::AbstractDict{String, <:TwoGridOperator}, Q::AbstractMatrix, A::AbstractMatrix; anorm::Bool=false)
    n = size(Q, 1)
    df = DataFrame(k = 1:n, normalized_frequency = (1:n) ./ (n + 1))
    for (name, tg) in tgs
        S = tg.S
        if anorm
            df[!, name] = [sqrt(dot(S * Q[:, k], A, S * Q[:, k]) / dot(Q[:, k], A, Q[:, k])) for k in 1:n]
        else
            df[!, name] = [norm(S * Q[:, k]) for k in 1:n]
        end
    end
    return df
end

function compute_XZ_bound(tg::TwoGridOperator, max_k::Int; tol::Float64=1e-12)
    evs = tg_spectrum(tg, :X; absolute=true, sorted=true)
    non_zeros = filter(x -> x > tol, evs)
    mu = isempty(non_zeros) ? 0.0 : first(non_zeros)
    cvr = sqrt(1.0 - mu)
    return [cvr^k for k in 0:max_k-1]
end

# -----------------------------------------------------------------------------
# 2. Problem & Operator Setup
# -----------------------------------------------------------------------------

n = 501
n_c = 250
ndomain = 10
overlap = 10
n_c_34 = div(3 * n, 4)

A = build_1d_laplacian(n)
λ_A, Q = build_fourier_basis(n)

# Smoothers (R)
R_j         = 0.5 * Matrix(I, n, n)
R_wJ        = (1.0 / 3.0) * Matrix(I, n, n)
R_gs        = inv(Matrix(LowerTriangular(A)))
R_as, _        = build_additive_schwarz(A, ndomain, overlap)
R_bJ, _        = build_additive_schwarz(A, ndomain, 0) # block Jacobi
ω_opt_34    = 0.6
R_wJ06      = (ω_opt_34 / 2.0) * Matrix(I, n, n)

# Transfer operators (P)
P_gmg = build_P_gmg(n, n_c)
P_34  = build_P_34(n, n_c_34)

# Theoretical analytical spectra
λ_S_true  = sort(abs.(1.0 .- 0.5 .* λ_A))
λ_Sw_true = sort(abs.(1.0 .- (1.0 / 3.0) .* λ_A))
λ_gs_true = sort(sort([cos( i * pi / (n+1))^2 for i in 1:n]))
λ_X_true  = sort(0.5 .* (2.0 .- 0.5 .* λ_A) .* λ_A)
λ_Xw_true = sort((1.0 / 3.0) .* (2.0 .- (1.0 / 3.0) .* λ_A) .* λ_A)

# -----------------------------------------------------------------------------
# 3. TwoGridOperator Collections
# -----------------------------------------------------------------------------

# Standard GMG transfer
tgs_gmg = Dict(
    "Jacobi"          => TwoGridOperator(A, R_j, P_gmg),
    "WeightedJacobi"  => TwoGridOperator(A, R_wJ, P_gmg),
    "GaussSeidel"     => TwoGridOperator(A, R_gs, P_gmg),
    "AdditiveSchwarz" => TwoGridOperator(A, R_as, P_gmg),
    "BlockJacobi" => TwoGridOperator(A, R_bJ, P_gmg)
)

tgs_gmg["GaussSeidel"].S

# Optimal transfer
tgs_opt = Dict(
    "Jacobi"         => TwoGridOperator(A, R_j, n_c),
    "WeightedJacobi" => TwoGridOperator(A, R_wJ, n_c),
    "GaussSeidel"    => TwoGridOperator(A, R_gs, n_c),
    "AdditiveSchwarz"=> TwoGridOperator(A, R_as, n_c),
    "BlockJacobi" => TwoGridOperator(A, R_bJ, n_c)
)

# 3:4 Coarsening transfer
tgs_34 = Dict(
    "Jacobi"            => TwoGridOperator(A, R_j, P_34),
    "WeightedJacobi_23" => TwoGridOperator(A, R_wJ, P_34),
    "WeightedJacobi_06" => TwoGridOperator(A, R_wJ06, P_34),
    "GaussSeidel"       => TwoGridOperator(A, R_gs, P_34),
    "AdditiveSchwarz"   => TwoGridOperator(A, R_as, P_34)
)

# -----------------------------------------------------------------------------
# 4. Spectral Experiments & DataFrames
# -----------------------------------------------------------------------------

# Smoother spectrum
df_sorted_abs_ev = tg_spectrum_dataframe(tgs_gmg, :S)
df_sorted_abs_ev[!, :Jacobi_true]         = λ_S_true
df_sorted_abs_ev[!, :WeightedJacobi_true] = λ_Sw_true

# Fourier attenuation
df_fourier_attenuation       = compute_fourier_attenuation(tgs_gmg, Q, A; anorm=false)
df_fourier_attenuation_anorm = compute_fourier_attenuation(tgs_gmg, Q, A; anorm=true)

# Symmetrized smoother T_bar spectrum
df_sorted_T_ev = tg_spectrum_dataframe(tgs_gmg, :T_bar)

# Two-grid error matrix X spectrum (Standard GMG vs. Optimal)
df_sorted_X_ev = tg_spectrum_dataframe(tgs_gmg, :X)

df_sorted_X_opt_ev = tg_spectrum_dataframe(tgs_opt, :X)
df_sorted_X_opt_ev[!, :Jacobi_true]         = λ_X_true
df_sorted_X_opt_ev[!, :WeightedJacobi_true] = λ_Xw_true

# -----------------------------------------------------------------------------
# 5. Two-Grid Convergence Studies & XZ Bounds
# -----------------------------------------------------------------------------

maxiter = 25
f = rand(n)
v0 = zeros(n)

# Solves with Standard GMG (Euclidean & Energy Norms)
df_convergence       = solve_dataframe(tgs_gmg, f, v0; maxiter=maxiter, normA=false)
df_convergence_normA = solve_dataframe(tgs_gmg, f, v0; maxiter=maxiter, normA=true)

# Append XZ theoretical bounds to convergence DataFrames
df_convergence[!, :bound_J]  = compute_XZ_bound(tgs_gmg["Jacobi"], maxiter)
df_convergence[!, :bound_wJ] = compute_XZ_bound(tgs_gmg["WeightedJacobi"], maxiter)
df_convergence[!, :bound_gs] = compute_XZ_bound(tgs_gmg["GaussSeidel"], maxiter)

df_convergence_normA[!, :bound_J]  = df_convergence[!, :bound_J]
df_convergence_normA[!, :bound_wJ] = df_convergence[!, :bound_wJ]
df_convergence_normA[!, :bound_gs] = df_convergence[!, :bound_gs]

# Solves with Optimal P
df_convergence_opt = solve_dataframe(tgs_opt, f, v0; maxiter=maxiter, normA=true)
df_convergence_opt[!, :bound_J]  = compute_XZ_bound(tgs_opt["Jacobi"], maxiter)
df_convergence_opt[!, :bound_wJ] = compute_XZ_bound(tgs_opt["WeightedJacobi"], maxiter)
df_convergence_opt[!, :bound_gs] = compute_XZ_bound(tgs_opt["GaussSeidel"], maxiter)

# -----------------------------------------------------------------------------
# 6. 3:4 Coarsening Spectral Analysis
# -----------------------------------------------------------------------------

df_34_smoother_spectrum      = tg_spectrum_dataframe(tgs_34, :S)
df_34_fourier_attenuation    = compute_fourier_attenuation(tgs_34, Q, A; anorm=false)
df_34_fourier_attenuation_a  = compute_fourier_attenuation(tgs_34, Q, A; anorm=true)
df_34_smoother_T_spectrum    = tg_spectrum_dataframe(tgs_34, :T_bar)
df_34_smoother_X_spectrum    = tg_spectrum_dataframe(tgs_34, :X)

# -----------------------------------------------------------------------------
# 7. CSV Export
# -----------------------------------------------------------------------------

CSV.write("experiments/data/smoother_spectrum.csv", df_sorted_abs_ev)
CSV.write("experiments/data/fourier_attenuation.csv", df_fourier_attenuation)
CSV.write("experiments/data/fourier_attenuation_anorm.csv", df_fourier_attenuation_anorm)
CSV.write("experiments/data/smoother_T_spectrum.csv", df_sorted_T_ev)
CSV.write("experiments/data/smoother_X_spectrum.csv", df_sorted_X_ev)
CSV.write("experiments/data/smoother_X_opt_spectrum.csv", df_sorted_X_opt_ev)
CSV.write("experiments/data/two_grid_convergence.csv", df_convergence)
CSV.write("experiments/data/two_grid_convergence_normA.csv", df_convergence_normA)
CSV.write("experiments/data/two_grid_convergence_opt.csv", df_convergence_opt)
CSV.write("experiments/data/df_34_smoother_spectrum.csv", df_34_smoother_spectrum)
CSV.write("experiments/data/df_34_fourier_attenuation.csv", df_34_fourier_attenuation)
CSV.write("experiments/data/df_34_fourier_attenuation_anorm.csv", df_34_fourier_attenuation_a)
CSV.write("experiments/data/df_34_smoother_T_spectrum.csv", df_34_smoother_T_spectrum)
CSV.write("experiments/data/df_34_smoother_X_spectrum.csv", df_34_smoother_X_spectrum)