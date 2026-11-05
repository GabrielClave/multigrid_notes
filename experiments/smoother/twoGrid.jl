using LinearAlgebra
using DataFrames

# -----------------------------------------------------------------------------
# 1. Structure Definition & Constructors
# -----------------------------------------------------------------------------

struct TwoGridOperator{TA, TR, TP, TS, TTbar, TPic, TX}
    A::TA
    R::TR
    P::TP
    S::TS
    T_bar::TTbar
    Π_c::TPic
    X::TX
end

"""
    TwoGridOperator(A, R, P)

Computes and stores derived operators:
  - S    = I - R*A
  - T_bar = (R' + R - R'*A*R) * A
  - Π_c  = P * (P' * A * P)^(-1) * P' * A
  - X    = (I - Π_c) * T_bar
"""
function TwoGridOperator(A::AbstractMatrix, R::AbstractMatrix, P::AbstractMatrix)
    n = size(A, 1)
    I_n = Matrix{eltype(A)}(I, n, n)

    # Smoother iteration matrix
    S = I_n - R * A

    # Symmetrized smoother operator T_bar = R_bar * A
    R_bar = R' + R - R' * A * R
    T_bar = R_bar * A

    # A-orthogonal projection onto Range(P)
    A_c = P' * A * P
    Π_c = P * (A_c \ (P' * A))

    # Two-grid operator on error space
    X = (I_n - Π_c) * T_bar

    return TwoGridOperator(A, R, P, S, T_bar, Π_c, X)
end

"""
    TwoGridOperator(A, R, n_c::Int)

builds the optimal prolongation P_opt for a coarse space of size `n_c` \n
using the first `n_c` A-orthonormal eigenvectors of T_bar (solving A*R_bar*A*v = λ*A*v).
"""
function TwoGridOperator(A::AbstractMatrix, R::AbstractMatrix, n_c::Int)
    # we want to solve R_bar * A * v = λv
    # which is equivalent to A * R_bar * A * v = λA * v
    # but that's better since M is supposed to be symmetric
    R_bar = R' + R - R' * A * R
    M = Symmetric(A * R_bar * A) # should we symmetrize or take the real part of eigenvectors ?

    # Generalized eigenvalue problem: M v = λ A v
    F = eigen(M, A)
    idx_opt = sortperm(F.values)
    Q_opt = F.vectors[:, idx_opt]

    P_opt = Q_opt[:, 1:n_c]
    return TwoGridOperator(A, R, P_opt)
end

# -----------------------------------------------------------------------------
# 2. Spectral Analysis Functions & DataFrames
# -----------------------------------------------------------------------------

"""
    tg_spectrum(tg::TwoGridOperator, op::Symbol; absolute=true, sorted=true)

Extracts the eigenvalues of a specific operator stored in `tg` (`:S`, `:T_bar`, or `:X`).
"""
function tg_spectrum(tg::TwoGridOperator, op::Symbol; absolute::Bool=true, sorted::Bool=true, tol = 1e-10)
    M = getfield(tg, op)
    evs, evect = eigen(M)
    condv = cond(evect)
    if (maximum(imag(evs)) < tol) println("complex eigenvalues, condition number of V: $condv") end
    if absolute
        evs = abs.(evs)
    end
    return sorted ? sort(evs) : evs
end

"""
    tg_spectrum_dataframe(tgs::AbstractDict{String, <:TwoGridOperator}, op::Symbol; kwargs...)

Computes the spectrum of `op` across a collection of operators and returns a DataFrame.
"""
function tg_spectrum_dataframe(tgs::AbstractDict{String, <:TwoGridOperator}, op::Symbol; kwargs...)
    n = size(first(values(tgs)).A, 1) # dimension
    df = DataFrame(k = 1:n)
    for (name, tg) in tgs
        df[!, name] = tg_spectrum(tg, op; kwargs...)
    end
    return df
end

# -----------------------------------------------------------------------------
# 3. Two-Grid Solvers & Convergence Tracking
# -----------------------------------------------------------------------------

"""
    two_grid_step!(tg, v, r, f, A_c_fact)

Executes one two-grid iteration on vector `v` with residual `r`.
"""
function two_grid_step!(tg::TwoGridOperator, v::AbstractVector, r::AbstractVector, f::AbstractVector, A_c_fact)
    # Residual: r = f - A*v
    # r .= f .- tg.A * v # assumed to have been computed

    # Coarse-grid correction: v = v + P * (P'*A*P)^(-1) * P' * r
    v .+= tg.P * (A_c_fact \ (tg.P' * r))

    # Post-smoothing: v = v + R * (f - A*v)
    r .= f .- tg.A * v
    v .+= tg.R * r
    return nothing
end

"""
    solve(tg, f, v0; maxiter=25, tol=1e-12, normA=false)

Runs the two-grid cycle and returns the normalized residual history.
"""
function solve(tg::TwoGridOperator, f::AbstractVector, v0::AbstractVector;
               maxiter::Int=25, tol::Float64=1e-12, normA::Bool=false)
    v = copy(v0)
    n = length(f)
    r = f .- tg.A * v

    A_c = tg.P' * tg.A * tg.P
    A_c_fact = factorize(A_c)

    res_hist = Vector{Float64}(undef, maxiter)

    compute_norm = (res) -> normA ? sqrt(dot(res, tg.A, res)) : norm(res)
    
    r_norm_0 = compute_norm(r)
    res_hist[1] = 1.0

    if r_norm_0 < tol
        return res_hist[1:1]
    end

    for k in 2:maxiter
        two_grid_step!(tg, v, r, f, A_c_fact)
        r .= f .- tg.A * v
        rel_res = compute_norm(r) / r_norm_0
        res_hist[k] = rel_res

        if rel_res < tol
            return res_hist[1:k]
        end
    end

    return res_hist
end

"""
    solve_dataframe(tgs::AbstractDict{String, <:TwoGridOperator}, f, v0; maxiter=25, normA=false)

Runs two-grid solves across a collection of operators starting from the same initial guess.
"""
function solve_dataframe(tgs::AbstractDict{String, <:TwoGridOperator}, f::AbstractVector, v0::AbstractVector;
                         maxiter::Int=25, normA::Bool=false)
    df = DataFrame(k = 1:maxiter)
    for (name, tg) in tgs
        df[!, name] = solve(tg, f, v0; maxiter=maxiter, normA=normA)
    end
    return df
end