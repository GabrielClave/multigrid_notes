using LinearAlgebra

include("twoGrid.jl")
include("twoGrid_helpers.jl")

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

S_as = I - R_as*A
S_wJ = I - R_wJ*A

# Transfer operators (P)
P_gmg = build_P_gmg(n, n_c)

# "randomized" P experiment

# build n x n_c random matrix V
V = rand(n, n_c)

# apply the smoother k times: V <- S^k * V
# V should be a basis of the low frequency space
P_rand = V
k = 1
for _ in 1:k P_rand = S_as*P_rand end
cond(P_rand) #1e16

tg_rand = TwoGridOperator(A, R_as, P_rand)
plot_spectrum(tg_rand, :X)
solve(tg_rand, rand(n), zeros(n))

k = 12
P_rand = V
for _ in 1:k P_rand = S_wJ*P_rand end
cond(P_rand) # k = 5, 1e4 / k = 10 1e6 / k = 100 1e9 and up
tg_rand = TwoGridOperator(A, R_wJ, P_rand)
plot_spectrum(tg_rand, :X)
# not even an orthogonal projection !

tgs_wJ = Dict(
    "gmg"  => TwoGridOperator(A, R_wJ, P_gmg),
    "rand" => TwoGridOperator(A, R_wJ, P_rand),
    "opt" => TwoGridOperator(A, R_wJ, n_c),
)

plot_spectrum(tgs_wJ, :X)
# not great
# probably not an actual basis: "too colinear" maybe ?

plot_history(tg_rand, rand(n), zeros(n))
plot_history(tgs_wJ, rand(n), zeros(n))

function build_smoothed_coarse_space(S, n_c, k)
    n = size(S, 1)
    # Start with a random orthogonal matrix
    Q, _ = qr(rand(n, n_c))
    Q = Matrix(Q)

    for _ in 1:k
        # Apply smoother
        V = S * Q
        # Re-orthonormalize
        Q, _ = qr(V)
        Q = Matrix(Q)
    end
    return Q
end

k = 20
P_smooth_wJ = build_smoothed_coarse_space(S_wJ, n_c, k)
cond(P_smooth_wJ) # logic 1

tgs_wJ = Dict(
    "gmg"  => TwoGridOperator(A, R_wJ, P_gmg),
    "rand" => TwoGridOperator(A, R_wJ, P_smooth_wJ),
    "opt" => TwoGridOperator(A, R_wJ, n_c),
)

plot_spectrum(tgs_wJ, :X)
plot_history(tgs_wJ, rand(n), zeros(n))
