using LinearAlgebra
using CSV
using DataFrames

include("twoGrid.jl")
include("twoGrid_helpers.jl")

# -----------------------------------------------------------------------------
# Problem & Operator Setup
# -----------------------------------------------------------------------------

n = 501
n_c = 250
ndomain = 10
overlap = 10

A = build_1d_laplacian(n)
λ_A, Q = build_fourier_basis(n)

# Smoothers (R)
R_j         = 0.5 * Matrix(I, n, n)
R_wJ        = (1.0 / 3.0) * Matrix(I, n, n)
R_gs        = inv(Matrix(LowerTriangular(A)))
R_as, _        = build_additive_schwarz(A, ndomain, overlap)
R_bJ, _        = build_additive_schwarz(A, 20, 0) # block Jacobi

# Transfer operators (P)
P_gmg = build_P_gmg(n, n_c)

# Theoretical analytical spectrum
λ_S_true  = sort(abs.(1.0 .- 0.5 .* λ_A))
λ_Sw_true = sort(abs.(1.0 .- (1.0 / 3.0) .* λ_A))
λ_gs_true = sort(sort([cos( i * pi / (n+1))^2 for i in 1:n]))
λ_X_true  = sort(0.5 .* (2.0 .- 0.5 .* λ_A) .* λ_A)
λ_Xw_true = sort((1.0 / 3.0) .* (2.0 .- (1.0 / 3.0) .* λ_A) .* λ_A)

# -----------------------------------------------------------------------------
# TwoGridOperator Collections
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

# -----------------------------------------------------------------------------
# Spectral Experiments & DataFrames
# -----------------------------------------------------------------------------

# Smoother spectrum
df_sorted_abs_ev = tg_spectrum(tgs_gmg, :S)
df_sorted_abs_ev[!, :Jacobi_true]         = λ_S_true
df_sorted_abs_ev[!, :WeightedJacobi_true] = λ_Sw_true
df_sorted_abs_ev[!, :GaussSeidel_true] = λ_gs_true

# Fourier attenuation
df_fourier_attenuation       = compute_fourier_attenuation(tgs_gmg, Q, A; anorm=false)
df_fourier_attenuation_anorm = compute_fourier_attenuation(tgs_gmg, Q, A; anorm=true)

# Symmetrized smoother T_bar spectrum
df_sorted_T_ev = tg_spectrum(tgs_gmg, :T_bar)

# Two-grid error matrix X spectrum (Standard GMG vs. Optimal)
df_sorted_X_ev = tg_spectrum(tgs_gmg, :X)

df_sorted_X_opt_ev = tg_spectrum(tgs_opt, :X)
df_sorted_X_opt_ev[!, :Jacobi_true]         = λ_X_true
df_sorted_X_opt_ev[!, :WeightedJacobi_true] = λ_Xw_true

# -----------------------------------------------------------------------------
# Two-Grid Convergence Studies & XZ Bounds
# -----------------------------------------------------------------------------

maxiter = 25
f = rand(n)
v0 = zeros(n)

# Solves with Standard GMG (Euclidean & Energy Norms)
df_convergence       = solve(tgs_gmg, f, v0; maxiter=maxiter, normA=false)
df_convergence_normA = solve(tgs_gmg, f, v0; maxiter=maxiter, normA=true)

# Append XZ theoretical bounds to convergence DataFrames
df_convergence[!, :bound_J]  = compute_XZ_bound(tgs_gmg["Jacobi"], maxiter)
df_convergence[!, :bound_wJ] = compute_XZ_bound(tgs_gmg["WeightedJacobi"], maxiter)
df_convergence[!, :bound_gs] = compute_XZ_bound(tgs_gmg["GaussSeidel"], maxiter)

df_convergence_normA[!, :bound_J]  = df_convergence[!, :bound_J]
df_convergence_normA[!, :bound_wJ] = df_convergence[!, :bound_wJ]
df_convergence_normA[!, :bound_gs] = df_convergence[!, :bound_gs]

# Solves with Optimal P
df_convergence_opt = solve(tgs_opt, f, v0; maxiter=maxiter, normA=true)
df_convergence_opt[!, :bound_J]  = compute_XZ_bound(tgs_opt["Jacobi"], maxiter)
df_convergence_opt[!, :bound_wJ] = compute_XZ_bound(tgs_opt["WeightedJacobi"], maxiter)
df_convergence_opt[!, :bound_gs] = compute_XZ_bound(tgs_opt["GaussSeidel"], maxiter)

# -----------------------------------------------------------------------------
# CSV Export
# -----------------------------------------------------------------------------

# CSV.write("experiments/data/smoother_spectrum.csv", df_sorted_abs_ev)
# CSV.write("experiments/data/fourier_attenuation.csv", df_fourier_attenuation)
# CSV.write("experiments/data/fourier_attenuation_anorm.csv", df_fourier_attenuation_anorm)
# CSV.write("experiments/data/smoother_T_spectrum.csv", df_sorted_T_ev)
# CSV.write("experiments/data/smoother_X_spectrum.csv", df_sorted_X_ev)
# CSV.write("experiments/data/smoother_X_opt_spectrum.csv", df_sorted_X_opt_ev)
# CSV.write("experiments/data/two_grid_convergence.csv", df_convergence)
# CSV.write("experiments/data/two_grid_convergence_normA.csv", df_convergence_normA)
# CSV.write("experiments/data/two_grid_convergence_opt.csv", df_convergence_opt)