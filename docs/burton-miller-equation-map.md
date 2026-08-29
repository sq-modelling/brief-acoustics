# Burton-Miller Equation-to-Code Map

**Release scope covered:** acoustic exterior Helmholtz problems with Dirichlet,
Neumann, or Robin boundary conditions
**Primary mathematical source:** Q. Sun and E. Klaseboer, "Non-Singular
Burton-Miller Boundary Element Method for Acoustics," *Fluids* 8 (2023), 56,
<https://doi.org/10.3390/fluids8020056>

This document describes the equations, signs, unknowns, and code locations used
by the Fortran solver's single-particle exterior Burton-Miller system.
Transmission and multilayer Burton-Miller systems are outside its supported
scope.

## 1. Implementation Summary

For a single exterior obstacle whose stored mesh normals point out of the solid
and into the acoustic fluid, the modern kernel and surface-operator modules
contain the seven raw operators required by Equations (15) and (16) of the 2023
paper.

The package solves this formulation using the `2N x 2N` augmented
system in Section 6. `fortran/test/TestBurtonMillerAssembly.f90` checks the
four blocks and boundary-condition elimination independently of mesh
integration, while the sphere validation scripts test the complete numerical
path.

The code-space coupling weight is
`w = i * exterior_free_term_sign * beta`, where `beta` is a length.
For the supported outward-from-solid mesh, `exterior_free_term_sign=-1`
and `w=-i*beta`. Section 7 describes the boundary-condition elimination used
for all three supported boundary conditions.

## 2. Conventions and Unknowns

Let

- `x0` be a collocation node;
- `x` be a source or quadrature point;
- `r_vector = x - x0` and `r = |r_vector|`;
- `n0` and `n` be the normals at `x0` and `x`;
- `phi` be the acoustic velocity potential;
- `q = dphi/dn`;
- `psi1` be the auxiliary harmonic field satisfying `psi1 = phi` on the
  boundary; and
- `s = dpsi1/dn` be the additional Burton-Miller unknown.

The code and paper use

$$
G_k(x_0,x)=\frac{\exp(ikr)}{r}, \qquad G_0(x_0,x)=\frac{1}{r}.
$$

There is no `1/(4*pi)` factor in either Green function. Consequently, the
smooth-boundary free terms appear as `4*pi`, not `1`.

The paper's normal points out of the exterior acoustic domain and into the
solid. The normal stored by the standard AU mesh points out of the solid and
into the acoustic domain. For that standard case,

$$
\boldsymbol n_{\rm mesh}=-\boldsymbol n_{\rm paper},\qquad
q_{\rm mesh}=-q_{\rm paper},\qquad s_{\rm mesh}=-s_{\rm paper},
$$

and `exterior_free_term_sign=-1`. In the rest of this document, unadorned operator symbols
refer to the stored-mesh-normal convention used by the modern code.

## 3. Published Equations in Operator Form

The following form is a reorganisation of Equations (4), (14), (15), and (16)
of the 2023 paper. The term proportional to
`phi(x0) * n0 dot (x-x0) * dG0/dn` on the right of Equation (14) has been moved
to its left. This is exactly how the code stores the operator.

### 3.1 Ordinary non-singular equation: Equation (4)

Define the regularised Helmholtz operators `H` and `G` such that

$$
H\phi=Gq. \tag{A}
$$

For one smooth exterior boundary in the paper convention, the left side
contains

$$
4\pi\phi_0+
\int_S\left[\phi\frac{\partial G_k}{\partial n}
-\phi_0\frac{\partial G_0}{\partial n}\right]dS,
$$

and the right side contains

$$
\int_S qG_k\,dS+
q_0\int_S\left[(\boldsymbol n_0\!\cdot\!\boldsymbol r)
\frac{\partial G_0}{\partial n}
-(\boldsymbol n_0\!\cdot\!\boldsymbol n)G_0\right]dS.
$$

### 3.2 Regularised normal-derivative equation: Equation (14)

Define `D`, `K`, and `L` so that

$$
D\phi=Kq+Ls. \tag{B}
$$

In the paper convention their actions are

$$
\begin{aligned}
D\phi={}&\int_S \phi
\left(\frac{\partial^2G_k}{\partial n\partial n_0}
-\frac{\partial^2G_0}{\partial n\partial n_0}\right)dS\\
&+\frac{k^2}{2}\phi_0\int_S
\left[(\boldsymbol n_0\!\cdot\!\boldsymbol r)
\frac{\partial G_0}{\partial n}
-(\boldsymbol n_0\!\cdot\!\boldsymbol n)G_0\right]dS,
\end{aligned}
$$

$$
Kq=\int_S q\frac{\partial G_k}{\partial n_0}dS
+q_0\int_S\frac{\partial G_0}{\partial n}dS-4\pi q_0,
$$

and

$$
Ls=-\int_S s\frac{\partial G_0}{\partial n_0}dS
-s_0\int_S\frac{\partial G_0}{\partial n}dS+4\pi s_0.
$$

### 3.3 Auxiliary Laplace equation: Equation (16)

The new unknown `s` is linked to `psi1=phi` through

$$
H_0\phi=G_0s. \tag{C}
$$

The operators are the non-singular Laplace counterparts of `H` and `G`:

$$
H_0\phi=4\pi\phi_0+
\int_S(\phi-\phi_0)\frac{\partial G_0}{\partial n}dS,
$$

$$
G_0s=\int_S sG_0\,dS+s_0\int_S
\left[(\boldsymbol n_0\!\cdot\!\boldsymbol r)
\frac{\partial G_0}{\partial n}
-(\boldsymbol n_0\!\cdot\!\boldsymbol n)G_0\right]dS.
$$

### 3.4 Burton-Miller combination: Equation (15)

Using residuals

$$
R_{\rm ordinary}=H\phi-Gq,
\qquad
R_{\rm normal}=D\phi-Kq-Ls,
$$

the paper combines its normal convention as

$$
R_{\rm ordinary,paper}+i\beta R_{\rm normal,paper}=0. \tag{D}
$$

The parameter `beta` has dimensions of length. The paper identifies the object
size and `1/k` as natural scales and uses `beta=1/k` in its rigid-sphere sweep.

## 4. Normal Transformation and Coupling Sign

For one boundary, write

$$
\boldsymbol n_{\rm mesh}=\sigma\boldsymbol n_{\rm paper},
\qquad \sigma\in\{-1,+1\}.
$$

In AU data, `sigma` is `layer(particle_id)%exterior_free_term_sign`. Reversing both source
and collocation normals gives the following operator transformation:

| Quantity | Mesh-normal value |
|---|---:|
| `q`, `s` | `sigma` times the paper value |
| `H`, `H0` | `sigma` times the paper operator |
| `G`, `G0`, `D` | unchanged |
| `K`, `L` | `sigma` times the paper operator |

The ordinary residual therefore changes by `sigma`, while the normal-derivative
residual does not. Equation (D), expressed entirely in code variables, is

$$
R_{\rm ordinary}+wR_{\rm normal}=0,
\qquad w=i\sigma\beta. \tag{E}
$$

For the supported outward-from-solid mesh, `sigma=-1`, so

$$
w=-i\beta.
$$

The function `burton_miller_coupling_weight` includes
`exterior_free_term_sign` to express this normal transformation explicitly.
The public Burton-Miller runner currently accepts only outward-from-solid meshes.

This simple row-wise transformation assumes that all boundary components of one
acoustic domain have a consistent relation to the domain normal. A nested
bounded region can contain components with opposite relations. That case needs
source-component orientation factors and is outside the supported scope.

## 5. Modern Operator Map

The modern implementation deliberately stores the unweighted raw operators.
The Burton-Miller coupling weight must left-multiply rows later in
`AU_Solver.f90`.

| Mathematical operator | Modern array | Meaning after global assembly |
|---|---|---|
| `G` | `external_g` | Regularised Helmholtz single-layer operator |
| `H` | `external_h` | Regularised Helmholtz double-layer operator plus the ordinary exterior free term |
| `K` | `bm_dg_dn0` | Collocation-normal derivative acting on `q`, including the oriented image of the paper's diagonal `-4*pi*q0` term |
| `D` | `bm_d2g_dndn0` | Regularised mixed second-normal-derivative operator acting on `phi` |
| `L` | `bm_minus_dg0_dn0` | Auxiliary static derivative operator acting on `s` |
| `G0` | `bm_g0` | Regularised static single-layer operator |
| `H0` | `bm_h0` | Regularised static double-layer operator plus its free term |
| `w_i` | `bm_coupling_weight(i)` | Coupling applied to collocation row `i` |

The array declarations are in
[`AU_SurfaceOperators.f90`](../fortran/src/AU_SurfaceOperators.f90).
Element kernels are accumulated in
[`AU_HelmholtzKernels.f90`](../fortran/src/AU_HelmholtzKernels.f90),
then scattered and corrected on the diagonal in
[`AU_SurfaceOperators.f90`](../fortran/src/AU_SurfaceOperators.f90).

### 5.1 Same-particle regularisation map

Let `same` mean that the source element belongs to the same particle as the
collocation node. The scalar corrections below are added to matrix entry
`(row_node,row_node)` after each same-particle element is integrated.

| Element result | Added to | Integral represented |
|---|---|---|
| `g_regularization` | `G` | `(n0 dot r) dG0/dn - (n0 dot n) G0` |
| `h_regularization` | `H` | `-dG0/dn` |
| `bm%dg0_dn_regularization` | `K` | `+dG0/dn` |
| `bm%h_regularization` | `D` | `(k^2/2)[(n0 dot r)dG0/dn-(n0 dot n)G0]` |
| `bm%minus_dg0_dn_regularization` | `L` | `-dG0/dn` |
| `bm%g0_regularization` | `G0` | `(n0 dot r)dG0/dn-(n0 dot n)G0` |
| `bm%h0_regularization` | `H0` | `-dG0/dn` |

The nodal parts paired with these corrections are:

- `g(:) = integral N_a Gk`;
- `h(:) = integral N_a dGk/dn`;
- `bm%dg_dn0(:) = integral N_a dGk/dn0`;
- `bm%d2g_dndn0(:) = integral N_a (d2Gk-d2G0)/(dn dn0)` on the same particle,
  and `integral N_a d2Gk/(dn dn0)` otherwise;
- `bm%minus_dg0_dn0(:) = -integral N_a dG0/dn0` on the same particle;
- `bm%g0(:) = integral N_a G0` on the same particle; and
- `bm%h0(:) = integral N_a dG0/dn` on the same particle.

The analytical diagonal free terms are assembled as follows. Define
`F_i = 4*pi*exterior_free_term_sign(i)`.

| Operator | Diagonal free term in modern code |
|---|---:|
| `H` | `+F_i` only when the exterior domain reaches infinity |
| `K` | `-F_i` |
| `L` | `+F_i` |
| `H0` | `+F_i` |
| `G`, `D`, `G0` | none |

For an isolated outward-normal solid, `F_i=-4*pi`. These placements agree with
the paper after the normal transformation in Section 4.

`internal_h` receives no free term in the existing ordinary bounded-domain
assembly, and no internal Burton-Miller operators are assembled. The solver rejects `use_burton_miller=.true.` for internal or transmission
equation rows. The table above specifies the external Burton-Miller path; it is
not an implicit derivation of an interior augmented formulation.

### 5.2 Kernel formula check

[`evaluate_green_values`](../fortran/src/AU_HelmholtzKernels.f90)
uses `r_vector=x-x0` and evaluates

$$
\frac{\partial G_k}{\partial n}
=-(\boldsymbol n\!\cdot\!\boldsymbol r)
(r^{-1}-ik)e^{ikr}r^{-2},
$$

$$
\frac{\partial G_k}{\partial n_0}
=(\boldsymbol n_0\!\cdot\!\boldsymbol r)
(r^{-1}-ik)e^{ikr}r^{-2}.
$$

Their opposite signs are required because differentiation with respect to `x0`
reverses the derivative of `r=|x-x0|`.

## 6. Augmented `2N x 2N` System

For `N` physical boundary nodes, define the row-diagonal matrix

$$
W=\operatorname{diag}(w_1,\ldots,w_N).
$$

First form

$$
H_c=H+WD,\qquad G_c=G+WK,\qquad L_c=WL. \tag{F}
$$

The two boundary equations are

$$
H_c\phi-G_cq-L_cs=0, \tag{G}
$$

$$
H_0\phi-G_0s=0. \tag{H}
$$

The multiplication by `W` is on the left because each collocation row has its
own coupling length and exterior free-term sign.

Let `u` contain the one physical unknown selected at every node. Express the
scattered boundary fields as

$$
\phi=\bar\phi+Pu,\qquad q=\bar q+Qu, \tag{I}
$$

where `P` and `Q` are diagonal coefficient matrices supplied by the boundary
condition. Substitution into Equations (G) and (H) gives

$$
\begin{bmatrix}
H_cP-G_cQ & -L_c\\
H_0P & -G_0
\end{bmatrix}
\begin{bmatrix}u\\s\end{bmatrix}
=
\begin{bmatrix}
-H_c\bar\phi+G_c\bar q\\
-H_0\bar\phi
\end{bmatrix}. \tag{J}
$$

These are the four `N x N` matrix blocks assembled by the solver. The
unknown order is

```text
x(1:N)       = physical boundary unknown u
x(N+1:2*N)   = auxiliary normal derivative s
```

The solved `s` field is mathematical support data, not acoustic pressure or
normal velocity. It may remain solver-local, but tests should retain it long
enough to check the Equation (H) residual.

## 7. Boundary-Condition Elimination

Boundary data are imposed on the total field,

$$
a(\phi_{sc}+\phi_{inc})+b(q_{sc}+q_{inc})=f.
$$

The table below gives `phi_bar`, `q_bar`, `P`, and `Q` in Equation (I).

| Condition and physical unknown `u` | `phi_bar` | `P` | `q_bar` | `Q` |
|---|---|---:|---|---:|
| Dirichlet, `u=q_sc` | `f-phi_inc` | `0` | `0` | `I` |
| Neumann, `u=phi_sc` | `0` | `I` | `f-q_inc` | `0` |
| Robin with `b != 0`, `u=phi_sc` | `0` | `I` | `(f-a*phi_inc-b*q_inc)/b` | `-(a/b)I` |
| Robin solved through `a != 0`, `u=q_sc` | `(f-a*phi_inc-b*q_inc)/a` | `-(b/a)I` | `0` | `I` |

Important special cases follow directly from Equation (J):

### 7.1 Dirichlet

$$
\begin{bmatrix}
-G_c & -L_c\\
0 & -G_0
\end{bmatrix}
\begin{bmatrix}q\\s\end{bmatrix}
=
\begin{bmatrix}-H_c\bar\phi\\-H_0\bar\phi\end{bmatrix}.
$$

### 7.2 Neumann

$$
\begin{bmatrix}
H_c & -L_c\\
H_0 & -G_0
\end{bmatrix}
\begin{bmatrix}\phi\\s\end{bmatrix}
=
\begin{bmatrix}G_c\bar q\\0\end{bmatrix}.
$$

### 7.3 Robin

No separate hand-written Burton-Miller Robin branch is needed. The existing
`phi_unknown_relation` and `dphi_dn_unknown_relation` results can be converted
to the diagonal `P`, `Q`, `phi_bar`, and `q_bar` values and passed through
Equation (J). This retains the `H0*P` contribution in the second block row
for either choice of the physical boundary unknown.

## 8. Supported Domain

The public Burton-Miller solver accepts:

- exterior Dirichlet, Neumann, or Robin equations only;
- one closed, connected particle;
- `exterior_free_term_sign=-1` with a verified outward-from-solid mesh;
- real, positive acoustic wavenumber; and
- no transmission condition.

The auxiliary operators `L`, `G0`, and `H0` are assembled only for source
elements on the collocation node's particle. The present validation covers
one particle; it does not establish a multiple-particle or multilayer
Burton-Miller formulation. Those cases are rejected by the solver.

## 9. Implementation and Validation

- [`AU_LayerTopology.f90`](../fortran/src/AU_LayerTopology.f90) defines the
  coupling length and the exterior-domain normal convention.
- [`AU_HelmholtzKernels.f90`](../fortran/src/AU_HelmholtzKernels.f90) evaluates
  the element kernels and regularisation terms.
- [`AU_SurfaceOperators.f90`](../fortran/src/AU_SurfaceOperators.f90)
  assembles the unweighted operators listed in Section 5.
- [`AU_Solver.f90`](../fortran/src/AU_Solver.f90) checks the supported domain,
  assembles Equation (J), and solves the augmented system.
- [`TestBurtonMillerAssembly.f90`](../fortran/test/TestBurtonMillerAssembly.f90)
  checks all four augmented blocks and the boundary-condition elimination.

The analytical sphere cases test Dirichlet, Neumann, and Robin data. The
manufactured ellipsoid checks a smooth non-spherical surface, and the
fictitious-frequency regression compares ordinary and Burton-Miller solutions
against the analytical sphere solution at the two published frequencies.

Run these checks through the public validation command:

```sh
bash scripts/validate.sh
```

Individual commands and the evidence limitations are listed in the
[README](../README.md#validation-results) and the case guides under `validation/`.
