import numpy as np
from fecore.verify import mms

def _sphere_points(n=200, R=1.0):
    rng = np.random.default_rng(0)
    p = rng.normal(size=(n,3)); p /= np.linalg.norm(p, axis=1)[:,None]
    return p*R

def test_surface_source_is_eigenpair_of_laplace_beltrami():
    # -Δ_sphere Y_l^m = l(l+1)/R^2 Y_l^m ; surface_source = eps*l(l+1)/R^2 * y42
    R, eps = 6.371e6, 8.854e-12
    p = _sphere_points(R=R)
    y = mms.y42_cart(p, R)
    s = mms.surface_source(p, R, eps)
    np.testing.assert_allclose(s, eps*20.0/R**2*y, rtol=1e-12)
    assert np.ptp(y) > 0  # non-trivial field

def test_shell_exact_satisfies_bcs():
    R, H = 6.371e6, 2.0e4
    ground = _sphere_points(R=R); top = _sphere_points(R=R+H)
    np.testing.assert_allclose(mms.shell_exact(ground, R, H), 0.0, atol=1e-9)  # Dirichlet
    # Neumann at top: radial derivative ~0 -> finite-diff check
    eps_r = 1.0
    outer = top*(1+eps_r/(R+H)); inner = top*(1-eps_r/(R+H))
    dphidr = (mms.shell_exact(outer, R, H) - mms.shell_exact(inner, R, H))/(2*eps_r)
    assert np.max(np.abs(dphidr)) < 1e-6 * np.max(np.abs(mms.shell_exact(_sphere_points(R=R+H/2), R, H)) + 1e-30)
