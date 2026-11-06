using LinearAlgebra
using DataFrames

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