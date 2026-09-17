// Import the Touying package and the Metropolis theme
#import "@preview/touying:0.5.3": *
#import themes.metropolis: *

// The classic LaTeX Metropolis theme uses Fira Sans. 
// (Make sure the font is installed on your system if compiling locally)
// #set text(font: "Fira Sans")arrow.l.r.doubledouble

// Initialize the theme
#show: metropolis-theme.with(
  aspect-ratio: "16-9",
  // You can customize the theme colors here if desired, 
  // but it defaults to the classic dark teal, orange, and off-white.
  config-info(
    title: [Multigrid Convergence and XZ Theory],
    subtitle: [Smoother and Coarse Grid Correction Interplay],
    author: [Your Name],
    date: datetime.today().display(),
    institution: [Your Institution],
  ),
)

#let Ri = $overline(R)^(-1)$

// Generate the title slide
#title-slide()

= Classic Convergence Theory

== Slide 1: Smooth error
$ S e &approx e \
  arrow.l.r.double (I - M^(-1)N)e &approx 0 \
  arrow.l.r.double (M - N)e &approx 0 \
  arrow.l.r.double A e = r &approx 0 $

== Slide 2: Oscillatory error
$ ||e||_A^2 = ( A e, e ) $

$ ( A e, e ) = sum_k lambda_k c_k^2 $

Huge if components are along the eigenvectors associated with large eigenvalues.

#pagebreak()

Convergence of a multigrid method requires two hypotheses:

== Slide 3: Smoothing property
*The Smoothing Property*: The smoother $S$ reduces high-frequency errors.
$ ||S e||_A^2 <= ||e||_A^2 - alpha ||A e||_(D^(-1))^2 $

- This implies $S$ is a contraction in the $A$-norm, and the relaxation would converge.
- The negative term quantifies the drop: 
  - If the error is smooth, $A e approx 0$, and progress stalls.
  - If the error is oscillatory, $||r||_(D^(-1))^2 = ||A e||_(D^(-1))^2 = sum_k lambda_k^2 c_k^2$ is huge and progress is fast.

== Slide 4: The approximation property
*The Approximation Property*: The coarse grid can approximate smooth errors.
$ min_(u_c) ||e - P u_c||_D^2 <= beta ||e||_A^2 $

(The distance between the error and $op("range")(P)$)
- If the error is oscillatory, $||e||_A$ is huge, and the inequality is easily satisfied, regardless of the quality of the interpolation.
- If $e$ is a smooth error, the interpolation needs to accurately represent it: smooth errors must be near $op("range")(P)$.

== Slide 5: Convergence theorem
*Convergence theorem* \
If $alpha <= beta$, the iteration converges, with $||S T G||_A <= sqrt(1 - alpha/beta)$.

#pagebreak()

= XZ Theory

== Slide 6: The promise 
// (I'll handle this one)

== Slide 7: Some notations (smoothers)
== Slide 8: Some notations (projections)
== Slide 9: Some notations (subspace decomposition)

#pagebreak()

== Slide 10: Two-Grids Operator

A two-grid correction operator is an operator $B : V' -> V$ that acts on $f in V'$ with two steps:
- Apply the subspace correction:
  $ w = P A_c^(-1)P' f $
- Apply the smoother:
  $ B f = w + R(f - A w) $

The associated error transfer operator is:
$ E = (I - R A)(I - Pi_c) $

#pagebreak()

== Slides 11 and 12: The operator X

We have:
$ ||E||_A^2 = 1 - lambda_"min"(X) $
where $ X = (I - Pi_c) overline(T) quad : quad W_P^(perp_A) -> W_P^(perp_A) $

$X$ represents the effect of the (symmetrized) smoother on the non-smooth errors.

Let $H = op("range")(P)^(perp_A)$ be the space of these high-frequency errors. \
We have $H = W_P^(perp_A) inter W$, and the effect of the two-grids cycle on $e in H$ is:
$ overline(E)|_H &= (I - overline(T))( I - Pi_c ) \
                 &= I - overline(T) quad (I - Pi_c = I "on" H) $

The components along $H$ are given by:
$ (I - Pi_c)overline(E)|_H &= I - (I - Pi_c)overline(T) \
                            &= I - X $

- If $X approx I$ on $H$, the smoother is effective and the error is almost entirely annihilated: \
  $lambda_"min"(X) approx 1$ and $||E||_A^2 = 1 - lambda_"min"(X) approx 0$
- If $X$ has a small eigenvalue, the associated direction in $H$ will escape both the effect of the smoother and the coarse space correction: convergence will stagnate.

#pagebreak()

==  Z

The inverse of $X$ on $W_P^(perp_A)$ can be explicitly written as:
$ Z = (overline(T))^(-1) (I - Q_c) quad : quad W_P^(perp_A) -> W_P^(perp_A) $

So $lambda_"min"(X) = 1 / lambda_"max"(Z)$.

*Interpretation of $Z$:* \
After applying the smoother, the remaining smooth error should be as close as possible to $op("range")(P)$. \
The distance is $d_(Ri) (e, op("range")(P)) = min_(e_P in op("range")(P)) ||e - e_P||_(Ri)^2 = ||(I - Q_c)v||^2_(Ri)$

$Z$ allows us to express this distance in the $A$-geometry:
$ (Z v, v)_A =  ||(I - Q_c)v||_(Ri)^2 $

== Z

And so in the worst-case scenario we have:
$ lambda_"max"(Z) &= max_(v in W_P^(perp_A)) r_Z (v) \
                  &= max_(v in W_P^(perp_A)) (Z v, v)_A / (v,v)_A \
                  &= max_(v in W_P^(perp_A)) (||(I - Q_c)v||_(Ri)^2) / (||v||_A^2) $

Which is also the maximum for $w in W$. So:
$ lambda_"max"(Z) = max_(v in W_P^(perp_A)) r_Z (v) = K(V_c) $

$ ||E||_A^2 = 1 - 1 / K(V_c) $

This can be understood as follows: the effectiveness of a two-grids method depends on how well the coarse space can approximate the smooth error (quantified by $Z$) that escaped the smoother (quantified by $X$).

#pagebreak()

== Slide 15: Optimal coarse space

The optimal two-grids convergence is when $||E||_A = 1 - mu_(n_c + 1)$. \
This is achieved by choosing $V_c$ so that $op("range")(P) = op("span"){q_1, ..., q_(n_c)}$, i.e., $H^(Ri) = op("span"){q_(n_c + 1), ..., q_n}$, since:

$ min_(w in op("span"){q_(n_c + 1), ..., q_n}) (||overline(R)A w||_(Ri)^2) / (||w||_(Ri)^2) = mu_(n_c + 1) $

This can be obtained by setting $V_c = RR^(n_c)$ and $P = (q_1, q_2, ..., q_(n_c))$.

// general idea: smoother / CGC interplay
// let's skip this for now

// classic convergence theory

// smooth error & oscillatory error
// slide 1: smooth error
// \begin{align*}
//       & Se \approx  e \\ 
//  \arrow.l.r.double & \left(I - M^{-1}N\right)e \approx  0 \\ 
//  \arrow.l.r.double & \left(M - N\right)e \approx  0 \\ 
//  \arrow.l.r.double & Ae = r \approx  0 \\ 
// \end{align*}
// 
// // slide 2: oscillatory error
// $\|e\|_{A} = \langle Ae,e \rangle$
// \[
// \langle Ae,e \rangle = \sum_{k}{\lambda_{k}c_{k}^{2}}
// \]
// huge if component are along the eigenvectors associated with large eigenvalue

// Convergence of a multigrid method requires two hypothesis:
// slide 3: smoothing property
// \textbf{The Smoothing Property}: The smoother $S$ reduces high-frequency errors
// \[
//     \|S e\|_{A}^2 \le \|e\|_{A}^2 - \alpha \|A e\|_{D^{-1}}^2
// \]

// \begin{itemize}
//     \item This implies $S$ is a contraction in the $A$-norm, and the relaxation would converge
//     \item the negative term quantifies the drop: if the error is smooth $A e \approx 0$, and progress stalls \\
//     if the error is oscillatory, $\|r\|_{D^{-1}} = \|Ae\|_{D^{-1}} = \sum_{k}{\lambda_{k}^2c_{k}^{2}}$ is huge and progress is fast
// \end{itemize}
// 
// slide 4: the approximation property
// \textbf{The Approximation Property}: The coarse grid can approximate smooth errors

// \[
// \min_{u_c} \|e - P u_c\|_D^2 \le \beta \|e\|_{A}^2
// \]
// % the distance between the error and range(P)

// \begin{itemize}
//     \item if the error is oscillatory, $\|e\|_{A}$ is huge, and the inequality is easily satisfied, regardless of the quality of the interpolation
//     \item if $e$ is a smooth error, the interpolation needs to accurately represent it: smooth errors must be near $\range(P)$
// \end{itemize}
// 
// slide 5: convergence theorem
// textbf{Convergence theorem}

// if $\alpha \le \beta$ the iteration converges, with $\| STG \|_A \le \sqrt{1 - \frac{\alpha}{\beta}}$

// XZ theory

// slide 6: the promise (I'll handle this one)
// slide7: some notations (smoothers)
// slide8: some notations (projections)
// slide9: some notations (subspace decomposition)
// slide 10: convergence theorem
// === Two-Grids Operator

// A two-grid correction operator is an operator $B : V' arrow.r.long V$ that acts on $f in V'$ with the two steps:
// - apply the subspace correction
// $
//   w &= P A_c^(-1)P' f // f or r ?
// $
// - apply the smoother
// $
//   B f &= w + R(f - A w)
// $

// and the associated error transfer operator is:
// $ E = (I - R A)(I - Pi_c) $
 
// slide 11 and 12: the operator X

// // we have:
// $
//   ||E||_A^2 &= 1 - lambda_(min)(X)
// $
// where $
//   X = (I - Pi_c) overline(T) quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A)
// $

// X represents the effect of the (symmetrized) smoother on the non-smooth errors.

// Let $H = range(P)^(perp_A)$ be the space of these high frequency errors.\
// We have $H = W_P^(perp_A) inter W$, and the effect of the two-grids cycle on $e in H$ is:
// $
//   overline(E)|_H &= (I - overline(T))( I - Pi_c ) \
//                  &= I - overline(T) quad (I - Pi_c = I "on" H)\
// $

// the composant along H are given by:
// $
//   (I- Pi_c)overline(E)|_H &= I - (I- Pi_c)overline(T) \
//                           &= I - X
// $

// if $X approx I "on" H$, the smoother is effective and the error is almost entirely annihilated: \
// $lambda_"min"(X) approx 1$ and $||E||_A^2 = 1 - lambda_"min"(X) approx 0$

// if X has a small eigenvalue, the associated direction in H will escape both the effect of the smoother and the coarse space correction: converge will stagnate.

// // slides 13 and 14:

// ==== Z

// The inverse of $X$ on $W_P^(perp_A)$ can be explicitly written as:
// $
//   Z = (overline(T))^(-1) (I - Q_c) quad : quad W_P^(perp_A) arrow.r.long W_P^(perp_A) 
// $

// $"so" lambda_min (X) = 1 / (lambda_max (Z)) $

// interpretation of $Z$:\
// after you apply the smoother the remaining smooth error should be as close as possible to $range(P)$\
// the distance is $ d_Ri (e,range(P)) = min_(e_P in range(P)) ||e - e_P||_Ri = ||(I - Q_c)v||^2_Ri$

// $Z$ allows to express this distance in the $A$-geometry:
// $
//   (Z v, v)_A =  ||(I - Q_c)v||_Ri^2 
// $
// and so in the worst case scenario we have:
// $
//   lambda_max (Z) &= max_(v in W_P^(perp_A)) r_Z (v) =  max_(v in W_P^(perp_A))( (Z v, v)_A ) / (v,v)_A\

//                  &= max_(v in W_P^(perp_A)) (||(I - Q_c)v||_Ri^2) / (||v||_A^2 )
// $

// which is also the maximum for $w in W$
// $ "so " lambda_max (Z)= max_(v in W_P^(perp_A)) r_Z (v) = K(V_c) $

// $ ||E||_A^2 = 1 - 1 / (K (V_c)) $

// which can be understood as follow: the effectiveness of a two-grids methods depends on how well the coarse space can approximate the smooth error (quantified by $Z$) that escaped the smoother (quantified by $X$)

// // // slides 15: optimal coarse space

// The optimal two-grids convergence is when $||E||_A = 1 - mu_(n_c + 1)$\
// which is achieved by choosing $V_c$ so that $range(P) = "vect"{q_1, ..., q_n_c}$, \ ie $H^Ri = "vect"{q_(n_c + 1), ..., q_n}$ since:
// $
//   min_(w in "vect"{q_(n_c + 1), ..., q_n})(||overline(R)A w||_Ri^2)/(||w||_Ri^2) = mu_(n_c + 1)
// $

// which can be obtained by setting $V_c = RR^(n_c)$ and $P = (q_1, q_2, ..., q_n)$