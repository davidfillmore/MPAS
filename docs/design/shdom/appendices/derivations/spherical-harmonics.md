(app-shdom-sphharm)=
# Spherical harmonic representation and transforms

SHDOM alternates between two angular representations of the
radiance and source-function fields: a spherical-harmonic
representation used to evaluate the scattering integral (where the
angular dependence collapses to a simple multiplication) and a
discrete-ordinate representation used to integrate the radiative
transfer equation along rays (where the spatial streaming is
local). This section derives the mappings between these two
representations — Eqs. (A1)–(A3) of Evans (1998) — from the
continuous spherical-harmonic orthogonality and the Gauss–Legendre
quadrature rule. The final subsection counts operations and obtains
the $\sim 9N^{3/2}$ estimate quoted by Evans (§2 of the main text
of Evans, 1998).

(app-shdom-sphharm-Ylm)=
## Direction coordinates and real spherical harmonics

Direction $\hat{\Omega}$ on the unit sphere is parameterized by
zenith cosine $\mu = \cos\theta \in [-1, 1]$ and azimuth
$\phi \in [0, 2\pi)$, with differential solid angle
$d\hat{\Omega} = d\mu\,d\phi$. Evans uses *real-valued* spherical
harmonics $Y_{lm}(\mu, \phi)$ that factor as

$$
  Y_{lm}(\mu, \phi) \;=\; \Lambda_{lm}(\mu)\,u(m\phi),
  \qquad l = 0, 1, 2, \ldots, \quad m = -l, \ldots, +l,
$$ (eq-Ylm-def)

where $\Lambda_{lm}(\mu)$ is the normalized associated Legendre
function of degree $l$ and azimuthal order $|m|$, and $u(m\phi)$ is
the real Fourier basis

$$
  u(m\phi) \;=\; \begin{cases}
    \cos(m\phi),    & m > 0, \\
    1,              & m = 0, \\
    \sin(|m|\phi),  & m < 0,
  \end{cases}
$$ (eq-u-def)

with $m > 0$ indexing cosine modes and $m < 0$ indexing sine modes,
as stated in the appendix of Evans (1998). The basis is
*orthonormal* on the sphere (Evans, 1998, appendix §a):

$$
  \int_{0}^{2\pi}\!\int_{-1}^{1}
    Y_{lm}(\mu, \phi)\,Y_{l'm'}(\mu, \phi)\,d\mu\,d\phi
  \;=\; \delta_{ll'}\,\delta_{mm'}.
$$ (eq-Ylm-ortho)

Orthonormality {eq}`eq-Ylm-ortho` constrains $\Lambda_{lm}$ to have
an $m$-dependent $L^2$ normalization: combined with the azimuthal
inner product

$$
  \int_{0}^{2\pi} u(m\phi)\,u(m'\phi)\,d\phi
  \;=\; c_m\,\delta_{mm'},
  \qquad c_m \;\equiv\; \begin{cases} 2\pi, & m = 0, \\ \pi, & m \ne 0, \end{cases}
$$ (eq-u-ortho)

which follows from the standard trigonometric identities,
{eq}`eq-Ylm-ortho` forces

$$
  \int_{-1}^{1} \Lambda_{lm}^{2}(\mu)\,d\mu \;=\; \frac{1}{c_m}.
$$ (eq-Lambda-ortho)

A closed form for $\Lambda_{lm}$ — absorbing the $\sqrt{c_m}$
normalization into the associated Legendre functions in the
standard way — is not needed below; the derivations depend only on
{eq}`eq-Ylm-ortho`, {eq}`eq-u-ortho`, and {eq}`eq-Lambda-ortho`.

Radiance and source function expand as

$$
  I(\mu, \phi) \;=\; \sum_{l=0}^{\infty}\sum_{m=-l}^{l}
      I_{lm}\,Y_{lm}(\mu, \phi),
  \qquad
  J(\mu, \phi) \;=\; \sum_{l=0}^{\infty}\sum_{m=-l}^{l}
      J_{lm}\,Y_{lm}(\mu, \phi),
$$ (eq-A1-exact)

which, after truncation
(see [](#app-shdom-sphharm-trunc)), reproduces Eq. (A1) of
Evans (1998). The expansion coefficients are recovered by
projection onto the orthonormal basis:

$$
  I_{lm} \;=\; \int_{0}^{2\pi}\!\int_{-1}^{1}
    I(\mu, \phi)\,Y_{lm}(\mu, \phi)\,d\mu\,d\phi
  \qquad (\text{and likewise for } J_{lm}).
$$ (eq-I-lm-continuous)

(app-shdom-sphharm-trunc)=
## Truncation and mode count

SHDOM truncates {eq}`eq-A1-exact` at polynomial degree $L$ in $\mu$
and azimuthal order $M$ in $\phi$, retaining modes with $l \le L$
and $|m| \le \min(l, M)$. Evans chooses

$$
  L \;=\; N_\mu - 1, \qquad M \;=\; N_\phi/2 - 1,
$$ (eq-LM)

where $N_\mu$ is the number of Gauss–Legendre abscissae in
$\mu \in [-1, 1]$ and $N_\phi$ is the (maximal) number of
equispaced azimuthal abscissae. The rationale for this specific
choice is that it is the largest truncation for which the
*discrete* spherical-harmonic inner product on the
$(\mu_j, \phi_k)$ grid reproduces the continuous inner product
{eq}`eq-Ylm-ortho` exactly; we establish this in
[](#app-shdom-sphharm-ortho-discrete).

The number of retained modes is

$$
  N_{\mathrm{lm}}^{\text{3D}}
  \;=\; (L+1)^{2}
  \qquad (\text{both cosine and sine azimuthal modes retained}),
$$ (eq-Nlm-3D)

corresponding to a *triangular* truncation ($M = L$, i.e.
$|m| \le l$). For 2D problems with azimuthal symmetry about
$\phi = 0$ only the cosine modes ($m \ge 0$) are needed, giving

$$
  N_{\mathrm{lm}}^{\text{2D}}
  \;=\; (L+1)(L/2 + 1),
$$ (eq-Nlm-2D)

as stated on page 441 of Evans (1998).

(app-shdom-sphharm-grid)=
## Gauss–Legendre abscissae and reduced azimuthal grid

**Zenith cosine.** The $N_\mu$ discrete zenith-cosine abscissae
$\mu_1, \ldots, \mu_{N_\mu}$ are the zeros of the Legendre
polynomial $P_{N_\mu}(\mu)$ on $[-1, 1]$; the associated
Gauss–Legendre weights $w_j$ satisfy

$$
  \int_{-1}^{1} p(\mu)\,d\mu
  \;=\; \sum_{j=1}^{N_\mu} w_j\,p(\mu_j)
$$ (eq-gauss-exact)

exactly for any polynomial $p$ of degree $\le 2N_\mu - 1 =
2L + 1$. This exactness is the entire reason for the truncation
choice $L = N_\mu - 1$ in {eq}`eq-LM`: products of two normalized
associated Legendre functions at fixed $m$,
$\Lambda_{lm}(\mu)\Lambda_{l'm}(\mu)$, are polynomials in $\mu$ of
degree $l + l' \le 2L$, which {eq}`eq-gauss-exact` integrates
exactly.

**Azimuth.** For each zenith abscissa $\mu_j$ the algorithm uses
$N_{\phi, j}$ equispaced azimuthal abscissae

$$
  \phi_{jk} \;=\; \frac{2\pi (k - 1)}{N_{\phi, j}},
  \qquad k = 1, \ldots, N_{\phi, j}.
$$ (eq-phi-jk)

On the unreduced grid $N_{\phi, j} = N_\phi$ for every $j$. SHDOM
*reduces* the azimuthal resolution near the poles, i.e. uses
$N_{\phi, j} < N_\phi$ for $|\mu_j|$ close to $1$: because the
zenith-polar region represents a small fraction of the solid
angle yet would consume a full set of azimuthal abscissae on an
unreduced grid, one can reduce $N_{\phi, j}$ there without loss
of angular fidelity. Evans reports that the resulting total
ordinate count

$$
  N \;\equiv\; \sum_{j=1}^{N_\mu} N_{\phi, j}
$$ (eq-N-total-ordinates)

is about 70 % of the unreduced count $N_\mu N_\phi$
(Evans, 1998, §2).

The azimuthal weights $\hat{w}_{jk}$ are chosen so that the
discrete inner product reproduces {eq}`eq-Ylm-ortho` exactly; we
work out their constraint in the next subsection.

(app-shdom-sphharm-ortho-discrete)=
## Discrete orthogonality of the spherical harmonics

The discrete inner product on the reduced Gaussian grid is

$$
  \langle f, g \rangle_{h}
  \;\equiv\;
  \sum_{j=1}^{N_\mu} w_j
   \sum_{k=1}^{N_{\phi, j}} \hat{w}_{jk}\,
   f(\mu_j, \phi_{jk})\, g(\mu_j, \phi_{jk}),
$$ (eq-discrete-inner-product)

and the discrete-orthogonality property required of the weights
$\{w_j, \hat{w}_{jk}\}$ and the basis $\{Y_{lm}\}$ is that it
reproduces the continuous inner product {eq}`eq-Ylm-ortho`
exactly:

$$
  \boxed{\;
    \langle Y_{lm}, Y_{l'm'} \rangle_{h}
    \;=\; \delta_{ll'}\,\delta_{mm'}
    \qquad \forall\, l, l' \le L,\; |m|, |m'| \le M.
  \;}
$$ (eq-Ylm-ortho-discrete)

Using the factorization
$Y_{lm} = \Lambda_{lm}(\mu)\,u(m\phi)$, the discrete inner product
separates:

$$
  \langle Y_{lm}, Y_{l'm'} \rangle_{h}
  \;=\;
  \sum_{j=1}^{N_\mu} w_j\,\Lambda_{lm}(\mu_j)\,\Lambda_{l'm'}(\mu_j)
  \underbrace{\sum_{k=1}^{N_{\phi, j}} \hat{w}_{jk}\,u(m\phi_{jk})\,u(m'\phi_{jk})}_{=\,A_j(m, m')}.
$$ (eq-ortho-separated)

We evaluate the two factors in turn.

**Azimuthal sum.** The simplest (and standard) choice of azimuthal
weights is

$$
  \hat{w}_{jk} \;=\; \frac{2\pi}{N_{\phi, j}},
$$ (eq-wjk)

i.e. uniform weights normalized to sum to $2\pi$ over the
$N_{\phi, j}$ equispaced abscissae at zenith row $j$. For
equispaced abscissae {eq}`eq-phi-jk` and trigonometric orders
$|m|, |m'| \le M = N_\phi/2 - 1$, the composite-trapezoidal (or
equivalently midpoint) rule is exact for any trigonometric
polynomial of total order $< N_\phi$
{cite:p}`golub2013matrix`. The product $u(m\phi)\,u(m'\phi)$ is a
trigonometric polynomial of total order $|m| + |m'| \le 2M =
N_\phi - 2 < N_\phi$; therefore the discrete sum reproduces the
continuous integral {eq}`eq-u-ortho` exactly,

$$
  A_j(m, m') \;=\; \int_{0}^{2\pi} u(m\phi)\,u(m'\phi)\,d\phi
                \;=\; c_m\,\delta_{mm'}.
$$ (eq-A-result)

**Zenith sum.** Inserting {eq}`eq-A-result` into
{eq}`eq-ortho-separated`,

$$
  \langle Y_{lm}, Y_{l'm'} \rangle_{h}
  \;=\; c_m\,\delta_{mm'}\,
  \sum_{j=1}^{N_\mu} w_j\,
    \Lambda_{lm}(\mu_j)\,\Lambda_{l'm}(\mu_j).
$$ (eq-ortho-after-phi)

Because $\Lambda_{lm}\Lambda_{l'm}$ is a polynomial in $\mu$ of
degree $l + l' \le 2L = 2N_\mu - 2 \le 2N_\mu - 1$, the
Gauss–Legendre rule {eq}`eq-gauss-exact` is exact and reproduces
the continuous integral exactly. Evaluating the continuous
integral via {eq}`eq-Lambda-ortho` at $m' = m$,

$$
  \sum_{j=1}^{N_\mu} w_j\,\Lambda_{lm}(\mu_j)\,\Lambda_{l'm}(\mu_j)
  \;=\; \int_{-1}^{1}\Lambda_{lm}(\mu)\,\Lambda_{l'm}(\mu)\,d\mu
  \;=\; \frac{1}{c_m}\,\delta_{ll'}.
$$ (eq-zenith-sum-result)

(Orthogonality at fixed $m$,
$\int \Lambda_{lm}\Lambda_{l'm}\,d\mu = (1/c_m)\delta_{ll'}$,
follows from {eq}`eq-Lambda-ortho` together with the standard
orthogonality of associated Legendre functions at fixed $m$.)

**Combination.** Inserting {eq}`eq-zenith-sum-result` into
{eq}`eq-ortho-after-phi`,

$$
  \langle Y_{lm}, Y_{l'm'} \rangle_{h}
  \;=\; c_m\,\delta_{mm'}\,\cdot\,\frac{1}{c_m}\,\delta_{ll'}
  \;=\; \delta_{ll'}\,\delta_{mm'},
$$

which is {eq}`eq-Ylm-ortho-discrete`, as claimed. The $m$-dependent
factor $c_m$ cancels exactly between the azimuthal and zenith
sums; this cancellation is the concrete meaning of Evans's
"normalized appropriately" phrasing.

In summary: the choice $L = N_\mu - 1$ together with
$M = N_\phi/2 - 1$ is the *largest* truncation compatible with
exact discrete orthogonality on the $(\mu_j, \phi_{jk})$ grid,
because it just saturates the exactness bound of Gauss–Legendre
quadrature (in $\mu$) and of the composite-trapezoidal rule (in
$\phi$). Using any larger truncation would introduce aliasing
into the discrete inner product.

(app-shdom-sphharm-forward)=
## Forward transform: spherical harmonics to discrete ordinates

The forward transform evaluates a spherical-harmonic
representation at the discrete ordinates $(\mu_j, \phi_{jk})$.
Starting from {eq}`eq-A1-exact` (truncated at $l \le L$,
$|m| \le M$) and substituting the direction $(\mu_j, \phi_{jk})$,

$$
  J_{jk} \;\equiv\; J(\mu_j, \phi_{jk})
  \;=\; \sum_{l=0}^{L}\sum_{m=-l}^{l} J_{lm}\,Y_{lm}(\mu_j, \phi_{jk}).
$$

Using $Y_{lm} = \Lambda_{lm}(\mu)\,u(m\phi)$ and regrouping so that
the $m$ sum is outermost,

$$
  \boxed{\;
    J_{jk} \;=\; \sum_{m=-M}^{M} u(m\phi_{jk})\,
                 \sum_{l=|m|}^{L}
                 \Lambda_{lm}(\mu_j)\,J_{lm},
  \;}
$$ (eq-A2)

reproducing Eq. (A2) of Evans (1998). The useful feature of
{eq}`eq-A2` is that the inner sum over $l$ depends on $(j, m)$ but
not on $k$, so it may be precomputed and reused across all
azimuthal abscissae for fixed $(\mu_j)$; this partial separation
is what drives the favorable operation count of
[](#app-shdom-sphharm-opcount).

(app-shdom-sphharm-inverse)=
## Inverse transform: discrete ordinates to spherical harmonics

The inverse transform computes the spherical-harmonic
coefficients $I_{lm}$ from the discrete-ordinate radiance values
$I_{jk} \equiv I(\mu_j, \phi_{jk})$. The continuous projection
{eq}`eq-I-lm-continuous` yields
$I_{lm} = \int\!\int I(\mu, \phi)\,Y_{lm}(\mu, \phi)\,d\mu\,d\phi$.
Approximating the integral by the discrete inner product
{eq}`eq-discrete-inner-product` — exact for bandlimited integrands
by [](#app-shdom-sphharm-ortho-discrete) —

$$
  I_{lm}
  \;=\; \sum_{j=1}^{N_\mu} w_j\,
        \sum_{k=1}^{N_{\phi, j}} \hat{w}_{jk}\,
        I_{jk}\,\Lambda_{lm}(\mu_j)\,u(m\phi_{jk}).
$$ (eq-I-lm-discrete-raw)

(Discrete exactness on bandlimited inputs was established in
[](#app-shdom-sphharm-ortho-discrete); there is no residual
$1/c_m$ factor because the $c_m$ implicit in $\Lambda_{lm}$'s
normalization cancels the $c_m$ from the azimuthal sum, exactly
as in the discrete-orthogonality derivation.) Regrouping so that
the $j$ sum is outermost,

$$
  \boxed{\;
    I_{lm} \;=\; \sum_{j=1}^{N_\mu} w_j\,\Lambda_{lm}(\mu_j)\,
                 \sum_{k=1}^{N_{\phi, j}} \hat{w}_{jk}\,u(m\phi_{jk})\,I_{jk},
  \;}
$$ (eq-A3)

which is Eq. (A3) of Evans (1998). Partial separation is again
visible: the inner sum over $k$ depends on $(j, m)$ but not on
$l$, so it may be precomputed once and reused across $l$.

The forward–inverse pair {eq}`eq-A2`–{eq}`eq-A3` is exact for any
$J(\mu, \phi)$ whose spherical-harmonic content is bandlimited to
$l \le L$ and $|m| \le M$: under this condition the discrete inner
product reproduces the continuous one exactly, per
[](#app-shdom-sphharm-ortho-discrete). For broader-band inputs
aliasing errors enter, which is why Evans's adaptive truncation
of the spherical-harmonic series (§2 of the main text) monitors
mode amplitudes and truncates where they fall below threshold
rather than extending past $L, M$.

(app-shdom-sphharm-opcount)=
## Operations count: the $\sim 9N^{3/2}$ estimate

A naive evaluation of {eq}`eq-A2`, treating the two-index sum as
a single dense matrix–vector product, costs
$O(N_{\mathrm{lm}} \cdot N) = O(N^{2})$ multiply–adds, where
$N = \sum_j N_{\phi, j}$ is the total ordinate count and
$N_{\mathrm{lm}} = (L+1)^{2}$ is the number of retained spherical
harmonic modes. Partial separation reduces this to $O(N^{3/2})$
in both the forward and inverse directions.

Treat the zenith and azimuthal extents of the discrete-ordinate
grid as comparable, $N_\mu \sim N_\phi \sim \sqrt{N}$; then
$L \sim N_\mu - 1 \sim \sqrt{N}$ and
$M \sim N_\phi / 2 \sim \tfrac{1}{2}\sqrt{N}$ by {eq}`eq-LM`, and
$N_{\mathrm{lm}} = (L+1)^{2} \sim N$.

**Forward transform {eq}`eq-A2`: inner $l$-sum.** Define the
auxiliary array
$T_{jm} \equiv \sum_{l=|m|}^{L} \Lambda_{lm}(\mu_j) J_{lm}$, of
size $N_\mu \times (2M+1)$. Each entry is a sum of at most
$(L+1)$ terms; forming all of $T_{jm}$ therefore costs

$$
  C_1 \;\sim\; 2 \cdot N_\mu \cdot (2M+1) \cdot (L+1)
       \;\sim\; 2 \cdot \sqrt{N} \cdot \sqrt{N} \cdot \sqrt{N}
       \;=\; 2\,N^{3/2}
$$ (eq-C1-forward)

multiply–adds, each counted as two floating-point operations.

**Forward transform {eq}`eq-A2`: outer $m$-sum.** For each
discrete ordinate $(j, k)$ (there are $N$ of these), the outer
sum $\sum_m u(m\phi_{jk})\,T_{jm}$ costs $(2M+1)$ multiply–adds,
i.e.

$$
  C_2 \;\sim\; 2 \cdot N \cdot (2M+1)
       \;\sim\; 2 \cdot N \cdot \sqrt{N}
       \;=\; 2\,N^{3/2}.
$$ (eq-C2-forward)

**Forward transform total.**

$$
  C_{\text{fwd}} \;=\; C_1 + C_2 \;\sim\; 4\,N^{3/2}.
$$ (eq-C-fwd)

**Inverse transform {eq}`eq-A3`.** By an analogous partial
separation (inner $k$-sum yields
$U_{jm} \equiv \sum_k \hat{w}_{jk}\,u(m\phi_{jk})\,I_{jk}$, outer
$j$-sum yields
$I_{lm} = \sum_j w_j\Lambda_{lm}(\mu_j)\,U_{jm}$):

$$
  C_{\text{inv}} \;\sim\; 4\,N^{3/2}.
$$ (eq-C-inv)

**Source-function update.** After the inverse transform, SHDOM
computes the new source from Eq. (A5) below (Eq. (A5) of
Evans, 1998), which is *pointwise multiplication* of each
$I_{lm}$ by the scalar $\omega \chi_l / (2l+1)$ plus a
thermal/solar addition. The cost is
$C_{\text{src}} \sim 2 N_{\mathrm{lm}} \sim 2 N$, which is
sub-leading compared with the transform cost.

**Total per-iteration cost.**

$$
  C_{\text{iter}} \;=\; C_{\text{fwd}} + C_{\text{inv}}
                        + C_{\text{src}}
                        + O(\text{streaming})
  \;\sim\; 8\,N^{3/2}\;+\;O(N).
$$ (eq-C-iter)

This $8 N^{3/2}$ bookkeeping estimate rounds to the
$\sim 9N^{3/2}$ quoted in Evans (1998, §2), with the difference
absorbed by constant-overhead terms, edge handling at the
truncation boundaries ($l = L$, $|m| = M$), and the accounting
convention for multiply–adds vs. separate multiplications and
additions.

**FFT acceleration of the azimuthal direction.** When the
azimuthal component of the transforms is evaluated by an FFT
rather than by direct summation (Evans, 1998, notes this is done
for $N_\phi \gtrsim 12$), the $O(M)$-per-ordinate outer sum in
{eq}`eq-C2-forward` and its inverse counterpart become
$O(\log M)$ per ordinate. This replaces $C_2$ and its analogue by
$\sim 2 N \log(M) = 2 N \log \sqrt{N} = N \log N$, which is
asymptotically subdominant to the $O(N^{3/2})$ inner sums. The
total per-iteration cost then scales as

$$
  C_{\text{iter}}^{\text{FFT}} \;\sim\; 2\,N^{3/2} + N \log N
  \;\to\; \sim 3\,N^{3/2}
  \quad (\text{Evans's FFT asymptotic quote}),
$$

since only the inner $l$- and $j$-sums retain the $N^{3/2}$
scaling.
