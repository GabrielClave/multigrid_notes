using LinearAlgebra

n = 50
n_c = 20

# 1D Laplacian
A = SymTridiagonal(2.0 * ones(n), -1.0 * ones(n - 1))

# eigenvalues and eigenvectors of A
k_indices = 1:n
λ_A = [4.0 * sin(k * π / (2 * (n + 1)))^2 for k in k_indices]

Q = [sqrt(2 / (n + 1)) * sin(i * k * π / (n + 1)) for i in 1:n, k in 1:n]

# Gauss-Seidel Smoother
D_L = LowerTriangular(A)
R = inv(D_L)

R_exact = [i >= j ? (0.5)^(i - j + 1) : 0.0 for i in 1:n, j in 1:n]
println("||R - R_exact||_∞ = ", norm(R - R_exact, Inf))

R_bar = R' + R - R' * A * R

# using the coarse space range(P) = eigenvectors of A with the n_c smallest eigenvalues
P = Q[:, 1:n_c]

Π_c = P*inv(P'*A*P)*P'

X = R_bar * A * (I - Π_c)

# is X A-self-adjoint ?
is_A_self_adjoint = isapprox(A * X, X' * A, atol=1e-12)
println("X is A-self-adjoint: ", is_A_self_adjoint) # not self adjoint

λ_X = real(eigen(X).values) # negative ev ! and even a 0
sort!(λ_X, rev=true)

# exact eigenvectors of X
λ_X_analytical = [2.0 / λ_A[k] for k in (n_c + 1):n]

println("Max difference with analytical X spectrum: ", 
        norm(sort(λ_X_modal[1:(n - n_c)]) - sort(λ_X_analytical), Inf))

Ri = inv(R_bar)

# Verify S-GS exact identity: R_bar^-1 = 1/2 * A^2
Ri_exact = 0.5 * (A^2)
println("||R_bar^-1 - 0.5 A^2||_∞ = ", norm(Ri - Ri_exact, Inf))

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