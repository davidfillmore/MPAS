"""Manufactured solutions for the FEM Poisson verification (Y_4^2)."""
import numpy as np

def y42_cart(xyz, R=1.0):
    """Real, unnormalized Y_4^2 sampled at Cartesian points on radius R.
    Y_4^2 ∝ cos^2(lat)(7 sin^2(lat)-1) cos(2 lon) (the prototype's convention)."""
    x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
    r = np.sqrt(x*x+y*y+z*z); lat = np.arcsin(np.clip(z/r, -1, 1)); lon = np.arctan2(y, x)
    return np.cos(lat)**2*(7.0*np.sin(lat)**2-1.0)*np.cos(2.0*lon)

def surface_source(xyz, R=1.0, eps=1.0):
    """-Δ_sphere Y_4^2 = l(l+1)/R^2 Y_4^2, l=4 -> source = eps*20/R^2 * Y_4^2."""
    return eps*20.0/R**2*y42_cart(xyz, R)

def _radial(xyz, R, H):
    r = np.linalg.norm(xyz, axis=1)
    zeta = (r - R)/H                                  # 0 at ground, 1 at top
    # Rr(0)=0, Rr'(1)=0: choose Rr = sin(pi/2 * zeta)
    return zeta, np.sin(0.5*np.pi*zeta)

def shell_exact(xyz, R, H):
    _, Rr = _radial(xyz, R, H)
    return y42_cart(xyz, R)*Rr

def shell_source(xyz, R, H, eps=1.0):
    """-∇^2(Y_4^2 Rr(r)) = [20/r^2 Rr - (Rr'' + 2/r Rr')] Y_4^2 (spherical Laplacian)."""
    r = np.linalg.norm(xyz, axis=1); zeta = (r-R)/H
    Rr = np.sin(0.5*np.pi*zeta)
    Rrp = (0.5*np.pi/H)*np.cos(0.5*np.pi*zeta)
    Rrpp = -(0.5*np.pi/H)**2*np.sin(0.5*np.pi*zeta)
    y = y42_cart(xyz, R)
    return eps*(20.0/r**2*Rr - (Rrpp + 2.0/r*Rrp))*y
