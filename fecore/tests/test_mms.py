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

def test_y42_matches_independent_spherical_harmonic():
    # scipy is not available in this env, so we verify y42_cart is the true
    # l=4, m=2 Laplace-Beltrami eigenfunction via a finite-difference spherical
    # Laplacian:  Δ_sphere f = -20/R^2 * f  must hold pointwise.
    #
    # Δ_sphere f = (1/sinθ) d/dθ(sinθ df/dθ) + (1/sin²θ) d²f/dφ²
    # where θ = colatitude (0 at N pole) and φ = longitude.
    # We evaluate on a grid away from the poles (0.15π < θ < 0.85π) and
    # away from nodes of f (|f| large enough to make the ratio stable).
    R = 6.371e6
    # Build a lat/lon grid away from poles.
    # 200x400 gives O(dth^2) FD error ~2.6e-4 < rtol=1e-3 (fast, ~0.1s).
    n_theta, n_phi = 200, 400
    theta = np.linspace(0.15*np.pi, 0.85*np.pi, n_theta)  # colatitude
    phi   = np.linspace(0.0, 2.0*np.pi, n_phi, endpoint=False)
    dth = theta[1] - theta[0]
    dph = phi[1]   - phi[0]

    # Helper: evaluate y42_cart from (theta, phi) on sphere of radius R
    def y42_tp(th, ph):
        sin_th = np.sin(th); cos_th = np.cos(th)
        x = R * sin_th * np.cos(ph)
        y = R * sin_th * np.sin(ph)
        z = R * cos_th
        xyz = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
        return mms.y42_cart(xyz, R).reshape(th.shape)

    TH, PH = np.meshgrid(theta, phi, indexing='ij')  # (n_theta, n_phi)
    f  = y42_tp(TH, PH)

    # First and second θ-derivatives via centered differences (interior only)
    # Use [1:-1] in theta direction
    f_th_fwd = y42_tp(TH+dth, PH)
    f_th_bwd = y42_tp(TH-dth, PH)
    d2f_dth2 = (f_th_fwd - 2.0*f + f_th_bwd) / dth**2
    df_dth   = (f_th_fwd - f_th_bwd) / (2.0*dth)

    # Second φ-derivative (periodic, centered differences)
    f_ph_fwd = y42_tp(TH, PH+dph)
    f_ph_bwd = y42_tp(TH, PH-dph)
    d2f_dph2 = (f_ph_fwd - 2.0*f + f_ph_bwd) / dph**2

    sin_th = np.sin(TH)
    cos_th = np.cos(TH)

    # Spherical Laplacian on sphere of radius R:
    # R^2 Δ_sphere f = d2f/dth2 + cos_th/sin_th * df/dth + 1/sin_th^2 * d2f/dph2
    lap_f = (1.0/R**2) * (d2f_dth2 + (cos_th/sin_th)*df_dth + (1.0/sin_th**2)*d2f_dph2)

    # For a Y_4^2 eigenfunction: Δ_sphere f = -20/R^2 * f
    expected = -20.0/R**2 * f

    # Only test at points where |f| is large enough that ratio is meaningful
    threshold = 0.1 * np.max(np.abs(f))
    mask = np.abs(f) > threshold
    assert mask.sum() > 50, "Too few non-zero test points — grid problem"

    np.testing.assert_allclose(
        lap_f[mask], expected[mask], rtol=1e-3,
        err_msg="y42_cart does NOT satisfy Δ_sphere y42 = -20/R^2 * y42: "
                "formula is NOT a true l=4,m=2 eigenfunction"
    )


def test_shell_exact_satisfies_bcs():
    R, H = 6.371e6, 2.0e4
    ground = _sphere_points(R=R); top = _sphere_points(R=R+H)
    np.testing.assert_allclose(mms.shell_exact(ground, R, H), 0.0, atol=1e-9)  # Dirichlet
    # Neumann at top: radial derivative ~0 -> finite-diff check
    eps_r = 1.0
    outer = top*(1+eps_r/(R+H)); inner = top*(1-eps_r/(R+H))
    dphidr = (mms.shell_exact(outer, R, H) - mms.shell_exact(inner, R, H))/(2*eps_r)
    assert np.max(np.abs(dphidr)) < 1e-6 * np.max(np.abs(mms.shell_exact(_sphere_points(R=R+H/2), R, H)) + 1e-30)
