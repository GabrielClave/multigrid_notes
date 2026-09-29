using LinearAlgebra

n = 50
n_c = 20

# 1D Laplacian
A = SymTridiagonal(2.0 * ones(n), -1.0 * ones(n - 1))

# eigenvalues and eigenvectors of A
k_indices = 1:n
λ_A = [4.0 * sin(k * π / (2 * (n + 1)))^2 for k in k_indices]

Q = [sqrt(2 / (n + 1)) * sin(i * k * π / (n + 1)) for i in 1:n, k in 1:n]
println("||Q'Q - I||_∞ = ", norm(Q'*Q - I, Inf)) # orthogonal basis: yes
println("||Q' A Q - I||_∞ = ", norm(Q'*A*Q - I, Inf)) # A-orthogonal basis: no
# do we want an A-orthogonal basis for the A geometry ?

# λ_A, Q = eigen(A)

# Gauss-Seidel Smoother
D_L = LowerTriangular(A)
R = inv(D_L)

R_exact = [i >= j ? (0.5)^(i - j + 1) : 0.0 for i in 1:n, j in 1:n]
println("||R - R_exact||_∞ = ", norm(R - R_exact, Inf))

R_bar = R' + R - R' * A * R

# Split Q into smooth (P) and high-frequency (Q_high) spaces
P = Q[:, 1:n_c] # using the coarse space range(P) = eigenvectors of A with the n_c smallest eigenvalues
Q_high = Q[:, (n_c + 1):n]       # Matrix of size n x (n - n_c)

Π_c = P*inv(P'*A*P)*P'*A

X = (I - Π_c)*R_bar * A  # here it is X: V -> V

# is X A-self-adjoint ?
is_A_self_adjoint = isapprox(A * X, X' * A, atol=1e-12)
println("X is A-self-adjoint: ", is_A_self_adjoint) # not self adjoint

λ_X = real(eigen(X).values) # n_c 0 eigenvalue: because nul(X) = n_c
sort!(λ_X, rev=true)

# we want a X: H -> H operator
# Q_high is the basis of H, Q_high'*Q_high = I
# any v ∈ H can be expressed as v = Q_high*c, c a vector of dim n - n_c 
# if Q_high' * A * Q_high = I, conversely c = Q_high' * A * v

X_H = (I - Π_c)*R_bar * A * Q_high
X_H*ones(n-n_c)

X_modal = Q_high' * (I - Π_c)*R_bar * A * Q_high # (n - n_c) x (n - n_c)

Λ_high = Diagonal(λ_A[(n_c + 1):n])
is_A_self_adjoint = isapprox(Λ_high * X_modal, X_modal' * Λ_high, atol=1e-12) # true
# println("X is A-self-adjoint: ", is_A_self_adjoint)

λ_X_modal = real(eigen(X_modal).values) # n_c 0 eigenvalue: because nul(X) = n_c
sort!(λ_X_modal, rev=true)

Ri = inv(R_bar)

# eigenvalues, eigenvectors = f(Ri)

# are eigenvectors orthonormal ?
# else orthogonalize

# using the optimal coarse space for Gauss-Seidel

# P = f1 to fn_c
# Π_c = P*inv(P'*A*P)*P'

# X = R_bar * A * (I - Π_c)

# is.symmetric(X)

# eigenvalues, eigenvectors = f(X)

# save eigenvalues

R_bar_test = R'*Diagonal(2*ones(n))*R
inv(R_bar_test)