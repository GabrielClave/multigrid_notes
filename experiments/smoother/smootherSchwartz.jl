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
issymmetric(R_as) #false
isapprox(R_as, R_as', atol=1e-12) #true

S_as = I - R_as * A
isapprox(A*S_as, (A*S_as)', atol=1e-12) #true
# S is A-self adjoint (almost)

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
issymmetric(X_modal_H) # false
println("X_modal_H is symmetric: ", isapprox(X_modal_H, X_modal_H', atol=1e-12)) #true

# λ_X_modal = sort(real(eigen(X_modal_full).values), rev=true)
λ_X_modal = eigen(X_modal_full).values
maximum(imag.(λ_X_modal)) # e-14

real(λ_X_modal)
λ_X_modal = sort(real(λ_X_modal), rev=true)

λ_X_modal_H = sort(real(eigen(X_modal_H).values), rev=true)

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
k_high_norm = [k/n for k in 1:n]

df = DataFrame(
    k_mode = 1:n,
    k_norm = k_high_norm,
    Modal = λ_X_modal,
    GMG = λ_X_gmg,
    XZ_Optimal = λ_X_opt
)

CSV.write("experiments/data/spectrum_X_Schwartz.csv", df)

# -------------------------------------------------------------
# Actual Multigrid algorithm
# -------------------------------------------------------------

function two_grid_iteration(CGC, e, maxiter = 500, tol = 1e-12)
    
    k = 0
    while k <= maxiter
        # if norm(e) > tol ||e||_2
        if sqrt(dot(e, A, e)) > tol # ||e||_A
            e = S_as*CGC*e
        else
            return k
        end
        k += 1
    end
    return maxiter
end

e = rand(n)
k_opt = two_grid_iteration(I-Π_c_opt, e) # 1 iterations
k_gmg = two_grid_iteration(I-Π_c_gmg, e) # 500 iteration
k_modal = two_grid_iteration(I-Π_c_modal, e) # 500 iteration

e_opt = S_as*(I-Π_c_opt)*e
norm(e_opt)
sqrt(dot(e_opt, A, e_opt))

e_gmg = S_as*(I-Π_c_gmg)*e
norm(e_gmg)
sqrt(dot(e_gmg, A, e_gmg))

e_gmg = S_as*(I-Π_c_gmg)*e_gmg
norm(e_gmg)
sqrt(dot(e_gmg, A, e_gmg))

function two_grid_iteration_reversed(CGC, e, maxiter = 500, tol = 1e-12)
    
    k = 0
    while k <= maxiter
        # if norm(e) > tol ||e||_2
        if sqrt(dot(e, A, e)) > tol # ||e||_A
            e = CGC*S_as*e
        else
            return k
        end
        k += 1
    end
    return maxiter
end

k_opt = two_grid_iteration_reversed(I-Π_c_opt, e) # 1 iterations
k_gmg = two_grid_iteration_reversed(I-Π_c_gmg, e) # 500 iteration
k_modal = two_grid_iteration_reversed(I-Π_c_modal, e) # 500 iteration

# using an actual multigrid algorithm

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

function two_grid_cycle(A, f, R, P, r, v; maxiter = 100, tol = 1e-9)

    A_c_inv = inv(P' * A * P)

    k = 0
    r = f .- A*v

    residuals = Float64[]
    nr = norm(r)
    push!(residuals, nr)

    while k <= maxiter

        if nr < tol
            print("converged in $k steps")
            return k, residuals
        end

        two_grid_step!(A, f, R, A_c_inv, P, r, v)
        k += 1
        r .= f .- A * v
        nr = norm(r)
        push!(residuals, nr)
    end

    print("did not converge")

    return maxiter, residuals
end

f = rand(n)
r = zeros(n)
v = zeros(n)

# optimal
k_opt, res_hist_opt = two_grid_cycle(A, f, R_as, P_opt, r, v)
#  1 step

# GMG
r = zeros(n)
v = zeros(n)
k_gmg, res_hist_gmg = two_grid_cycle(A, f, R_as, P_gmg, r, v)
# no convergence

# modal
r = zeros(n)
v = zeros(n)
k_gmg, res_hist_gmg = two_grid_cycle(A, f, R_as, P_gmg, r, v)
# no convergence


# -------------------------------------------------------------
# Fourier attenuation
# -------------------------------------------------------------

F_as = eigen(S_as) # complex
evalues = F_as.values
maximum(imag.(evalues)) # e-15
evect = F_as.vectors
maximum(imag.(evect)) # 0.6
cond(evect) #1e7