using LinearAlgebra

include("P_genEO.jl")
include("../twoGrid.jl")

# setup

Nx, Ny = 10, 4       # N = 40 subdomains
mx, my = 4, 4       # 5x5 nodes per subdomain
A, A_locals, R_list, D_list = build_subdomains(Nx, Ny, mx, my; c=1e-3, contrast=1e4)

# keep the p lowest eigenvalues
p = 5
P_0, V_0, alpha = build_coarse_basis(A, A_locals, R_list, D_list, p)
n = size(A)[1]
size(P_0)[2]

# One-level Additive Schwarz preconditioner M_AS
M_AS = sum(R_list[i]' * inv(R_list[i] * A * R_list[i]') * R_list[i] for i in eachindex(A_locals))
println("cond(A)         = ", cond(A))
println("cond(M_AS * A)  = ", cond(M_AS * A))

# Jacobi, Gauss-Seidel
R_j         = inv( Diagonal(A) )
# R_wJ        = (1.0 / 4.0) * inv( Diagonal(A) )
R_gs        = inv(Matrix(LowerTriangular(A)))

n_c = div(n,2)

# Optimal transfer
tgs_opt = Dict(
    "Jacobi"         => TwoGridOperator(A, R_j, n_c),
    # "WeightedJacobi" => TwoGridOperator(A, R_wJ, n_c),
    "GaussSeidel"    => TwoGridOperator(A, R_gs, n_c),
    "AdditiveSchwarz"=> TwoGridOperator(A, M_AS, n_c),
)

plot_spectrum(tgs_opt, :R)
plot_spectrum(tgs_opt, :T_bar)
plot_spectrum(tgs_opt, :X)

plot_history(tgs_opt, rand(n), rand(n))
df_hist = solve(tgs_opt, rand(n), rand(n))

tgs_AS = Dict(
    "optimal" => TwoGridOperator(A, M_AS, n_c),
    "genEO" => TwoGridOperator(A, M_AS, P_0),
)

plot_spectrum(tgs_AS, :R)
plot_spectrum(tgs_AS, :T_bar)
plot_spectrum(tgs_AS, :X)

plot_history(tgs_AS, rand(n), rand(n))

tg_genEO = TwoGridOperator(A, M_AS, P_0)
plot_spectrum(tg_genEO, :X)