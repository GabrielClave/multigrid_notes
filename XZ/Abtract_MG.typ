// Document setup & page styling
#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  numbering: "1",
)

// Font configuration (using XCharter / Charter equivalent)
#set text(
  font: "xCharter",
  size: 11pt,
  lang: "en",
)

// Paragraph spacing
#set par(justify: true, leading: 0.65em)

// Custom Math Operators
#let range = math.op("range")
#let nul = math.op("null")
#let im = math.op("Im")
#let ker = math.op("ker")
#let Ri = $overline(R)^(-1)$

// Document Title & Metadata
#align(center)[
  #text(size: 18pt, weight: "bold")[Abstract Multigrid notes] \
  #v(8pt)
  #text(size: 12pt)[Gabriel Clave] \
  #v(8pt)
  #datetime.today().display("[month repr:long] [day], [year]")
]

#v(1.5em)

// Table of Contents
#outline(indent: auto)

#v(2em)

= Variational Formulation

We consider the problem
$ A u = f $
with $A: V -> V'$ SSPD, $u in V$ and $f in V'$.

$u$ and $f$ do not live in the same vector spaces:

- $u$ is in $V tilde.equiv RR^n$, and represents a primal state (temperature, pressure), and the natural norm is the energy norm $\|x\|^2 = x^T A x$.
- $f$ is in $V' tilde.equiv (RR^n)'$, and represents a source (force, flux), and the natural norm is the dual norm $\|b\|^2 = b^T A^(-1) b$.

// In particular x^T b does not really make sense

It is this problem that arises naturally from the weak formulation.

We solve
$ "Find " u in V " such that " chevron.l A u, v chevron.r = chevron.l f, v chevron.r quad forall v in V $

// Noting that Au is in V', so <Au, v> is the real number (Au)(v)

== Smoother

=== general linear iteration

The first component of a multigrid solver is a smoother.\
It is a linear operator $B : V' arrow.r.long V$ that can be used to solve (equation) through the linear iteration:
$
  v^((k+1)) = v^((k)) + B(f - A v^((k)))
$
with associated error propagation operator
$
  E = I - B A
$
If $B$ is not symmetric, it can be symmetrized through the iteration:

with
$
  overline(E) = I - overline(B) A
$

=== a two-grid smoother

Let $ R : V' arrow.r.long V$ be a valid smoother. \
If $e in ker(A)$, then $E e = (I - R A)e = e$: nothing happens.\
The smoother alone is unable to reduce the composant of the error that are in ker(A)

similarly, if $e_lambda in E_lambda (A)$, then $E e_lambda = (I - lambda R)e_lambda$\
for a given $R$ with bounded norm,  $E e_lambda arrow.r.long_(lambda arrow.r 0) e_lambda $

It will be difficult for the smoother to remove composant of the error in directions associated with the smallest eigenvalues of $A$.

The role of R is to act as an approximate inverse of A, but not on the entire spectrum.\
The application of the smoother will damp error composants in the directions associated with the highest eigenvalues of $A$.\
We will call this vector space $H$, the space of high frequency error.\
The error not in $H$ will need to be addressed separately.

==== smoother geometry

$Ri$ is invertible by design and symmetric: it defines a scalar product on $V$ $(dot,dot)_(Ri).$

Let $overline(T) = overline(R)A$

We have $ker(A) subset ker(overline(T))$

$overline(T)$ acts as a bridge between the $A$ geometry and the $Ri$ geometry
$ 
  (u,v)_A = (overline(R) A u, v)_(Ri) = (overline(T) u,v)_(Ri)
$

we also have:

$
  x in ker(overline(T)) & arrow.r.double ||overline(T)x||_(Ri) = 0 \
                         & arrow.r.double ||x||_A = 0 \
                         & arrow.r.double x in ker(A)
$

which means $ker(A) = ker(overline(T))$ and so $range(overline(T)) = W$
// W is not defined yet

== Coarse space correction

In order to deal with the error composant unaffected by the smoother, the idea is to use a coarse space $V_c$, linked with V by the injective operator $P: V_c -> V$

We will then look to approximate the solution of (equation) in the range of $P$, which means finding a $u_P in range(P)$, such that $u approx u_P = P u_c, quad u_c in V_c$.
// need smoother link: Galerkin projection allows to remove composant in range(P)

We are doing a Galerkin projection: instead of solving the variational problem on the entire space V, we restrict the search space to $range(P) in V$:

Find $u_P in range(P)$ such that
$
                        &chevron.l A u_P, v_P chevron.r = chevron.l f, v_P chevron.r quad forall v_P in range(P) \
  arrow.l.r.double.long &chevron.l A P u_c, P v_c chevron.r = chevron.l f, P v_c chevron.r quad  forall v_c in V_c \
  arrow.l.r.double.long &chevron.l A P u_c, P v_c chevron.r = chevron.l f, P v_c chevron.r quad forall v_c in V_c \
  arrow.l.r.double.long &chevron.l P' A P u_c, v_c chevron.r = chevron.l P' f, v_c chevron.r quad forall v_c in V_c
$

Noting $A_c = P' A P$ and $f_c = P'f$, we are solving the equivalent system on $V_c$
$ A_c u_c = f_c $

If we define $e_("approx") = u - u_P$, we can note that $chevron.l A e_("approx"), v_P chevron.r = 0$ for all $v_P in range(P)$, i.e., the approximation error is orthogonal to $range(P)$.

// a word on solving the coarse problem, and designing an AMG method

= Convergence Theory

Let us consider $A$ SSPD. \
We will consider the case $V_c subset V$.\

A has a non trivial kernel $N$.\
The smoother $E = I - R A$ is blind to anything in $N$, so the coarse space correction has to deal with the error in $N$: we need
$
 N subset range(P)
$

Using W = $N^(perp_A)$ we decompose $V$ into
$
  V = N xor_A W 
$

and we define $ W_P = range(P) inter W$\
meaning we have
$
  range(P) = W_P xor_A N = P(V_c)
$

$P$ is injective from $V_c$ to $range(P)$, so $V_c$ can be decomposed into
$
  V_c = W_c xor_A N_c
$
with $W_c = P^(-1)(W_P)$ and $N_c = P^(-1)(N)$

#v(1cm)

We define the following $A$-orthogonal projections:

$
  Pi_c : V &arrow.r.long range(P)\
  Pi_W : V &arrow.r.long W
$

and we will note that when $w in W_P = range(P) inter W$, we have $Pi_c (w) = Pi_W (w)$

#v(1cm)

We define the following $Ri$-orthogonal projections:

$
 Q_c: V &arrow.r.long range(P)\
 Q_W: V &arrow.r.long W
$
#v(1cm)
=== Two-Grids Operator

A two-grid correction operator is an operator $B : V' arrow.r.long V$ that acts on $f in V'$ with the two steps:
- apply the subspace correction
$
  w &= P A_c^(-1)P' f // f or r ?
$
- apply the smoother
$
  B f &= w + R(f - A w)
$

and the associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

the convergence rate of the method is given by:
$ ||E||_A^2 = 1 - 1 / (K (V_c)) $

with
$
  K (V_c) = max_(w in W) (|| (I - Q_c)v||_(Ri)^2 ) / ( ||v||_A^2) =
  max_(w in W) min_(v_P in range(P)) (|| w - v_P||_(Ri)^2 ) / ( ||v||_A^2)
$
#v(1cm)
=== Proof

We notice that $||(I - R A) v||_A^2 = ((I - overline(R) A) v, v)_A$ for all $v in V$.

Then we have:
$
  ||E||_A^2 &= max_(w in W) (||(I - R A)(I - Pi_c) w||_A^2) / (||w||_A^2) \
  &= max_(w in W) (((I - overline(R) A)(I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) \
  &= 1 - min_(w in W) ((overline(R) A (I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) \
  &= 1 - min_(w in W) ((Q_W overline(R) A (I - Pi_c) w, (I - Pi_c) w)_A) / (||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2) \
  &= 1 - min_(v in W_P^(perp_A)) ((Q_W overline(R) A v, v)_A) / (||v||_A^2) \
  &= 1 - min_(v in W_P^(perp_A)) (((I - Pi_c) Q_W overline(R) A v, v)_A) / (||v||_A^2) \
  &= 1 - lambda_(min)(X),
$

where $ X = (I - Pi_c) Q_W overline(R) A quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A) $

#v(1cm)

additional details:

$
  ((overline(R) A (I - Pi_c) w, (I - Pi_c) w)_A) / (||w||_A^2) = ((Q_W overline(R) A (I - Pi_c) w, (I - Pi_c) w)_A) / (||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2)
$

$overline(R) A (I - Pi_c) w in range(overline(T)) = W$ so it is unchanged by $Q_W$ \
$||w||_A^2 = ||(I - Pi_c) w||_A^2 + ||Pi_c w||_A^2$ using the orthogonal decomposition of W

#v(1cm)

$
  (Q_W overline(R) A v, v)_A = ((I-Pi_c)Q_W overline(R) A v, v)_A
$

$Q_W overline(R) A v in W$, so it is composed of some part along $W_P$ that will vanish when applying the scalar product with $v in W_P^perp_A$: \
$ forall w in W, v in W_P^(perp_A), quad(w,v)_A = ((I-Pi_c)w,v)_A $ 

#v(1cm)

==== X

$
  X = (I - Pi_c) overline(T) quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A)
$

X is self-adjoint with respect to $(dot, dot)_A$, as $I-Pi_c$ is an $A$-orthogonal projection and $overline(R)$ is SPD.\
The $Q_W$ projection is redundant here as in our setting $im(overline(T)) = W$

X represents the effect of the (symmetrized) smoother on the non-smooth errors.

Let $H = range(P)^(perp_A)$ be the space of these high frequency errors.\
We have $H = W_P^(perp_A) inter W$, and the effect of the two-grids cycle on $e in H$ is:
$
  overline(E)|_H &= (I - overline(T))( I - Pi_c ) \
                 &= I - overline(T) quad (I - Pi_c = I "on" H)\
$

the composant along H are given by:
$
  (I- Pi_c)overline(E)|_H &= I - (I- Pi_c)overline(T) \
                          &= I - X
$

if $X approx I "on" H$, the smoother is effective and the error is almost entirely annihilated: \
$lambda_"min"(X) approx 1$ and $||E||_A^2 = 1 - lambda_"min"(X) approx 0$

if X has a small eigenvalue, the associated direction in H will escape both the effect of the smoother and the coarse space correction: converge will stagnate.

==== Z

The inverse of $X$ on $W_P^(perp_A)$ can be explicitly written as:
$
  Z = (overline(T))^(-1) (I - Q_c) quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A) 
$

using $(u,v)_A = (overline(T) u,v)_(Ri)$ we have indeed, $forall u in W, forall v_P in W_P$:

$
  (Z u, v)_A &= ( (I - Q_c) u, v )_(Ri) = 0
$
meaning $im(Z) subset W_P^(perp_A)$, and we have// don't we need the equality here ?

$
  X Z = (I - Pi_c)(I - Q_c) = I "on" W_P^(perp_A)
$

$"so" lambda_min (X) = 1 / (lambda_max (Z)) $

interpretation of $Z$:\
after you apply the smoother the remaining smooth error should be as close as possible to $range(P)$\
the distance is $ d_Ri (e,range(P)) = min_(e_P in range(P)) ||e - e_P||_Ri = ||(I - Q_c)v||^2_Ri$

$Z$ allows to express this distance in the $A$-geometry:
$
  (Z v, v)_A = (overline(T)^(-1)(I - Q_c)v, v)_A &= (overline(T) overline(T)^(-1)(I - Q_c)v, v)_Ri \
                                    &= ||(I - Q_c)v||_Ri^2 
$
and so in the worst case scenario we have:
$
  lambda_max (Z) &= max_(v in W_P^(perp_A)) r_Z (v) =  max_(v in W_P^(perp_A))( (Z v, v)_A ) / (v,v)_A\

                 &= max_(v in W_P^(perp_A)) (||(I - Q_c)v||_Ri^2) / (||v||_A^2 )
$

which is also the maximum for $w in W,"decomposed into" w = v + v_P$, $v in W_P^(perp_A) "and" v_P in W_P$
$
  r_Z (w) = (||(I - Q_c)(v + v_P)||_Ri^2) / (||v||_A^2 + ||v_P||_A^2) lt.eq (||(I - Q_c)v||_Ri^2) / (||v||_A^2)
$
$ "so " lambda_max (Z)= max_(w in W) r_Z (w) = max_(v in W_P^(perp_A)) r_Z (v) = K(V_c) $

$ ||E||_A^2 = 1 - 1 / (K (V_c)) $

which can be understood as follow: the effectiveness of a two-grids methods depends on 
how well the coarse space can approximate the smooth error (quantified by $Z$) that escaped the smoother (quantified by $X$)

=== Optimal Coarse Space

We want to chose the optimal coarse space $V_c$ so that $range(V_c)$ is the best complement to the space of high frequency errors H.

This is achieved when $K(V_c)$ is minimal.
$
  1 / (K(V_c)) = min_(w in W) max_(w_P in W_P) (||w||_A^2)/(||w - w_P||_Ri^2)
$
If we take the minimum over the space of high frequency error $H^Ri = (W_P^(perp_Ri) inter W) subset W$,\ for which $min_(w in H^Ri) ||w - w_P||_Ri^2 = ||w||_Ri^2 $ // dessin ?
$
  1 / (K(V_c)) &lt.eq min_(w in H^Ri) max_(w_P in W_P) (||w||_A^2)/(||w - w_P||_Ri^2) \
  &= min_(w in H^Ri)(||w||_A^2)/(||w||_Ri^2) \
  &= min_(w in H^Ri)(||overline(R)A w||_Ri^2)/(||w||_Ri^2) \
  &= min_(w in H^Ri) r_(overline(R)A) (w)
$

and so
$
  max_(dim V_c = n_c) 1 / (K(V_c)) &= max_(dim V_c = n_c) quad min_(w in H^Ri) quad r_(overline(R)A) (w)\
  &= max_(dim H^Ri = n - (n_c + 1) + 1) quad min_(w in H^Ri) quad r_(overline(R)A) (w)\
  &= mu_(n_c + 1)
$
by the max-min theorem (Courant-Fisher), with ${mu_j , q_j}$ the eigenpairs of $overline(R)A$ in increasing order.

The optimal two-grids convergence is then $||E||_A = 1 - mu_(n_c + 1)$\
which is achieved by choosing $V_c$ so that $range(P) = "vect"{q_1, ..., q_n_c}$, \ ie $H^Ri = "vect"{q_(n_c + 1), ..., q_n}$ since:
$
  min_(w in "vect"{q_(n_c + 1), ..., q_n})(||overline(R)A w||_Ri^2)/(||w||_Ri^2) = mu_(n_c + 1)
$

which can be obtained by setting $V_c = RR^(n_c)$ and $P = (q_1, q_2, ..., q_n)$
// un exemple pour rigoler ?