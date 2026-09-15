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
#let im = math.op("im")

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

Let $ R : V' arrow.r.long V$ be a valid smoother
// valid ? + role

Let $overline(T) = overline(R)A$

We have $nul(A) subset nul(overline(T))$

noting that $overline(T)$ acts as a bridge between the $A$ geometry and the $overline(R)^(-1)$ geometry
$ 
  (u,v)_A = (overline(R) A u, v)_(overline(R)^(-1)) = (overline(T) u,v)_(overline(R)^(-1))
$

we also have:

$
  x in nul(overline(T)) & arrow.r.double ||overline(T)x||_(overline(R)^(-1)) = 0 \
                         & arrow.r.double ||x||_A = 0 \
                         & arrow.r.double x in nul(A)
$

which means $nul(A) = nul(overline(T))$ and so $range(overline(T)) = W$

== Coarse space correction


The idea is to use a coarse space $V_c$, linked with V by the injective operator $P: V_c -> V$ \
We will then look to approximate the solution in the range of $P$, which means finding a $u_P in range(P)$, such that $u approx u_P = P u_c, quad u_c in V_c$.

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

$P$ is injective from $V_c$ to $range(P)$, so it is decomposed into
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

We define the following $overline(R)^(-1)$-orthogonal projections:

$
 Q_c: V &arrow.r.long range(P)\
 Q_W: V &arrow.r.long W
$
#v(1cm)
=== Operator

A two-grid correction operator is an operator $B : V' arrow.r.long V$ composed of the two steps:
- apply the subspace correction
- apply the smoother

$
  w &= P A_c^(-1)P' f\ // f or r ?
  B f &= w + R(f - A w)
$

and the associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

the convergence rate of the method is given by:
$ ||E||_A^2 = 1 - 1 / (K (V_c)) $

with
$
  K (V_c) = max_(w in W) (|| (I - Q_c)v||_(overline(R)^(-1))^2 ) / ( ||v||_A^2) =
  max_(w in W) min_(v_P in range(P)) (|| w - v_P||_(overline(R)^(-1))^2 ) / ( ||v||_A^2)
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

where $ X = (I - Pi_c) Q_W overline(R) A : W_P^(perp_A) arrow.r.long W_P^(perp_A) $

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

X is self-adjoint with respect to $(dot, dot)_A$.

One key observation is that the inverse of $X$ on $W_P^(perp_A)$ can be explicitly written as:
$
  Z = (Q_W overline(R) A)^(-1) (I - Q_c)
$
