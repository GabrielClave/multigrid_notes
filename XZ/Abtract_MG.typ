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
with $A: V -> V'$, $u in V$ and $f in V'$.

$u$ and $f$ do not live in the same vector spaces:

- $u$ is in $V tilde.equiv RR^n$, and represents a primal state (temperature, pressure), and the natural norm is the energy norm $\|x\|^2 = x^T A x$.
- $f$ is in $V' tilde.equiv (RR^n)'$, and represents a source (force, flux), and the natural norm is the dual norm $\|b\|^2 = b^T A^(-1) b$.

// In particular x^T b does not really make sense

It is this problem that arises naturally from the weak formulation.

We solve
$ "Find " u in V " such that " chevron.l A u, v chevron.r = chevron.l f, v chevron.r quad forall v in V $

// Noting that Au is in V', so <Au, v> is the real number (Au)(v)

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
