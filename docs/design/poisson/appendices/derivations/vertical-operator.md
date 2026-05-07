(app-vertical)=
# Vertical stencil and boundary treatment

## Interior vertical flux

The upper face of the prism $\Omega_{i, k}$ is a horizontal
polygon at altitude $z^{\mathrm{f}}_{i, k + \hlf}$, area $A_i$,
outward normal $+\hat{\bm{z}}$. The integrated outward flux is

$$
  \mathcal{F}^{+}_{i, k}
  \;=\; \int_{\mathrm{top\,face}}
    \varepsilon (\nabla \varphi) \cdot \hat{\bm{z}} \; dA
  \;=\; \varepsilon \int_{\mathrm{top\,face}}
        \partial_z \varphi \; dA.
$$ (eq-B-top-flux-int)

With one-point quadrature at the face centroid and
finite-difference approximation for $\partial_z \varphi$ at the
midpoint between layer-midpoint altitudes $z_{i, k}$ and
$z_{i, k+1}$,

$$
  \mathcal{F}^{+}_{i, k}
  \;\approx\; \varepsilon \, A_i \,
  \frac{\varphi_{i, k+1} - \varphi_{i, k}}{\dz{k+\hlf}},
  \qquad \dz{k+\hlf} = z_{i, k+1} - z_{i, k}.
$$ (eq-B-top-flux)

Analogously, the bottom face has inward normal $+\hat{\bm{z}}$ at
altitude $z^{\mathrm{f}}_{i, k - \hlf}$, so the inward-flux form
(which is what enters the control-volume equation with a $-$
sign, cf. {eq}`eq-A-split`) is:

$$
  \mathcal{F}^{-}_{i, k}
  \;\approx\; \varepsilon \, A_i \,
  \frac{\varphi_{i, k} - \varphi_{i, k-1}}{\dz{k-\hlf}}.
$$ (eq-B-bot-flux)

Substituting into equation {eq}`eq-A-discrete-eq` and rearranging:

$$
\begin{aligned}
  \mathcal{F}^{+}_{i, k} - \mathcal{F}^{-}_{i, k}
  &\;=\; \varepsilon \, A_i \,
       \left[
         \frac{\varphi_{i, k+1} - \varphi_{i, k}}{\dz{k+\hlf}}
         - \frac{\varphi_{i, k} - \varphi_{i, k-1}}{\dz{k-\hlf}}
       \right] \\
  &\;=\; \varepsilon \, A_i \,
         \frac{\varphi_{i, k+1} - \varphi_{i, k}}{\dz{k+\hlf}}
         + \varepsilon \, A_i \,
         \frac{\varphi_{i, k-1} - \varphi_{i, k}}{\dz{k-\hlf}}.
\end{aligned}
$$ (eq-B-vert-combined)

The last form, which groups the contributions symmetrically in
$\varphi_{i, k+1}$ and $\varphi_{i, k-1}$, is the form used in
equation {eq}`eq-L-full`.

(app-vertical-matrix)=
## Stencil contributions to $\mathbb{L}$

From {eq}`eq-B-vert-combined`, the vertical contribution to row
$(i, k)$ of $\mathbb{L}$ is:

- off-diagonal $(i, k+1)$: entry
  $+\varepsilon A_i / \dz{k+\hlf}$;
- off-diagonal $(i, k-1)$: entry
  $+\varepsilon A_i / \dz{k-\hlf}$;
- diagonal $(i, k)$: additional contribution of
  $-\varepsilon A_i / \dz{k+\hlf} - \varepsilon A_i / \dz{k-\hlf}$.

## Ground-Dirichlet boundary ($k = 1$)

At the ground face $z = z_{\mathrm{g}}$, the Dirichlet condition
{eq}`eq-bc-ground` prescribes $\varphi = 0$. The bottom face of
cell $(i, 1)$ has face altitude
$z^{\mathrm{f}}_{i, \tfrac{1}{2}} = z_{\mathrm{g}}$; let
$\dz{\tfrac{1}{2}}^{\mathrm{half}} \equiv z_{i, 1} - z_{\mathrm{g}}$
denote the half-layer distance from the ground face to the
layer-1 midpoint. The inward flux at the ground face becomes

$$
  \mathcal{F}^{-}_{i, 1} \;=\;
  \int_{\mathrm{ground\,face}}
     \varepsilon (\nabla \varphi) \cdot \hat{\bm{z}} \, dA
  \;\approx\; \varepsilon A_i \,
  \frac{\varphi_{i, 1} - \varphi_{\mathrm{g}}}{\dz{\tfrac{1}{2}}^{\mathrm{half}}}
  \;=\; \varepsilon A_i \,
  \frac{\varphi_{i, 1} - 0}{\dz{\tfrac{1}{2}}^{\mathrm{half}}}
  \;=\; \varepsilon A_i \,
  \frac{\varphi_{i, 1}}{\dz{\tfrac{1}{2}}^{\mathrm{half}}}.
$$ (eq-B-ground-flux)

Substituting this into {eq}`eq-A-discrete-eq` for $k = 1$:

$$
  \mathcal{F}^{h}_{i, 1} + \mathcal{F}^{+}_{i, 1}
  - \varepsilon A_i \, \frac{\varphi_{i, 1}}{\dz{\tfrac{1}{2}}^{\mathrm{half}}}
  \;=\; -\rho_{i, 1} V_{i, 1}.
$$

Moving the Dirichlet contribution to the left side:

$$
  \mathcal{F}^{h}_{i, 1} + \mathcal{F}^{+}_{i, 1}
  \;+\; \mathcal{F}^{(\mathrm{D})}_{i, 1} \varphi_{i, 1}
  \;=\; -\rho_{i, 1} V_{i, 1},
  \qquad
  \mathcal{F}^{(\mathrm{D})}_{i, 1}
  \;\equiv\; -\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}}.
$$ (eq-B-ground-L-RHS)

Therefore the $k = 1$ row of $\mathbb{L}$ gets an additional
diagonal contribution
$-\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}}$ compared to
the interior-row formula. No off-diagonal contribution to a
$k = 0$ ghost cell is needed, since the Dirichlet value
$\varphi_{\mathrm{g}} = 0$ contributes nothing to any row.

**Generalization for $\varphi_{\mathrm{g}} \ne 0$.** For a
prescribed nonzero ground potential $\varphi_{\mathrm{g}}$, the
ground-face flux becomes
$\varepsilon A_i (\varphi_{i, 1} - \varphi_{\mathrm{g}}) / \dz{\tfrac{1}{2}}^{\mathrm{half}}$,
and equation {eq}`eq-B-ground-L-RHS` acquires a right-hand-side
correction:

$$
  \mathcal{F}^{h}_{i, 1} + \mathcal{F}^{+}_{i, 1}
  \;+\; \mathcal{F}^{(\mathrm{D})}_{i, 1} \varphi_{i, 1}
  \;=\; -\rho_{i, 1} V_{i, 1}
  \;-\; \varepsilon A_i \, \varphi_{\mathrm{g}}
        / \dz{\tfrac{1}{2}}^{\mathrm{half}}.
$$ (eq-B-ground-L-RHS-nonzero)

The right-hand side gains the term
$-\varepsilon A_i \varphi_{\mathrm{g}} / \dz{\tfrac{1}{2}}^{\mathrm{half}}$.
This is implemented generically in the code so that a nonzero
ground potential can be set at runtime.

## Top Neumann boundary ($k = K$)

At the model top $z = z_{\mathrm{top}}$, the Neumann condition
$\partial_z \varphi = 0$ implies zero outward flux. Therefore

$$
  \mathcal{F}^{+}_{i, K} \;=\; 0.
$$

The $k = K$ row of $\mathbb{L}$ drops the entry connecting
$(i, K)$ to $(i, K+1)$ (there is no such cell) and also drops the
corresponding diagonal contribution
$-\varepsilon A_i / \dz{K+\hlf}$. Only the lower vertical face
and all horizontal faces contribute.

## Summary

Combining all contributions, the assembled form of row $(i, k)$
is exactly the form given in the main-text stencil summary
appendix and repeated here:

$$
  (\mathbb{L} \bphi)_{i, k}
  \;=\; \sum_{e \in \partial i}
          \varepsilon_e \, \frac{\elen \dz{k}}{\dcen}
          \, (\varphi_{j, k} - \varphi_{i, k})
        + \varepsilon \, A_i \,
          \frac{\varphi_{i, k+1} - \varphi_{i, k}}{\dz{k+\hlf}}
        + \varepsilon \, A_i \,
          \frac{\varphi_{i, k-1} - \varphi_{i, k}}{\dz{k-\hlf}},
$$

with the boundary adjustments at $k = 1$ and $k = K$ described
above.
