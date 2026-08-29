"""
analytic.py
===========

Analytical solution for acoustic scattering from a rigid sphere.

This module provides the reference solution used to validate the mini-AU
boundary element solver.

Problem
-------
A time-harmonic incident plane wave scatters from a rigid sphere of radius a.

The incident wave is

    phi_inc(x) = exp(i k d · x),

where

    k = wavenumber,
    d = propagation direction.

For the default mini-AU example, d = (1, 0, 0), so the incident wave is

    phi_inc = exp(i k x).

The total field is

    phi_total = phi_inc + phi_scattered.

Rigid boundary condition
------------------------
For a rigid sphere, the total normal derivative vanishes on r = a:

    d(phi_total) / dr = 0.

This is the sound-hard Neumann condition.

Analytical expansion
--------------------
Use spherical coordinates centered on the sphere.  Let theta be the polar angle
measured from the incident-wave direction.  Then

    mu = cos(theta).

The incident wave expansion is

    exp(i k r cos(theta))
    = sum_{n=0}^{infinity} (2n + 1) i^n j_n(k r) P_n(mu),

where

    j_n = spherical Bessel function,
    P_n = Legendre polynomial.

The outgoing scattered field is

    phi_scattered
    = sum_{n=0}^{infinity} A_n h_n^{(1)}(k r) P_n(mu),

where h_n^{(1)} is the outgoing spherical Hankel function.

The rigid boundary condition gives

    A_n = -(2n + 1) i^n j_n'(k a) / h_n^{(1)'}(k a).

This module uses scipy.special for spherical Bessel and Legendre functions.

Sign convention
---------------
This analytical solution uses the outgoing Green-function convention consistent
with exp(i k r) / r and incident wave exp(i k d · x).

If your numerical code uses the opposite time or spatial convention, the complex
conjugate or sign of k may need to be adjusted.

Run this file directly to perform simple self-checks:

    python analytic.py
"""

from __future__ import annotations

import numpy as np
from scipy.special import eval_legendre, spherical_jn, spherical_yn
import matplotlib.pyplot as plt
from matplotlib import cm

Array = np.ndarray


def plot_rigid_sphere_surface_field(
    quantity: str = "phi_total",
    component: str = "real",
    k: float = 1.0,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
    n_theta: int = 120,
    n_phi: int = 240,
) -> None:
    """
    Visualize analytical fields on the sphere surface.

    This routine is for human inspection only.  The validation scripts below use
    the numerical arrays returned by rigid_sphere_scattered_field and
    rigid_sphere_total_field instead of this interactive plot.

    Parameters
    ----------
    quantity:
        "phi_inc", "phi_sca", or "phi_total".

    component:
        "real", "imag", "abs", or "phase".
    """
    theta = np.linspace(0.0, np.pi, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)

    # Build a latitude-longitude grid on the sphere:
    # theta is the polar angle, phi is the azimuthal angle.
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")

    x = radius * np.sin(theta_grid) * np.cos(phi_grid)
    y = radius * np.sin(theta_grid) * np.sin(phi_grid)
    z = radius * np.cos(theta_grid)

    points = np.column_stack(
        [x.ravel(), y.ravel(), z.ravel()]
    )

    if quantity == "phi_inc":
        values = incident_plane_wave(
            points,
            k=k,
            direction=direction,
        )
    elif quantity == "phi_sca":
        values = rigid_sphere_scattered_field(
            points,
            k=k,
            radius=radius,
            direction=direction,
            n_terms=n_terms,
        )
    elif quantity == "phi_total":
        values = rigid_sphere_total_field(
            points,
            k=k,
            radius=radius,
            direction=direction,
            n_terms=n_terms,
        )
    else:
        raise ValueError("quantity must be phi_inc, phi_sca, or phi_total")

    # Convert the one-dimensional list of point values back to the same shape as
    # the surface grid so Matplotlib can color the sphere.
    values = values.reshape(theta_grid.shape)

    if component == "real":
        field = np.real(values)
    elif component == "imag":
        field = np.imag(values)
    elif component == "abs":
        field = np.abs(values)
    elif component == "phase":
        field = np.angle(values)
    else:
        raise ValueError("component must be real, imag, abs, or phase")

    fmin = np.min(field)
    fmax = np.max(field)

    if abs(fmax - fmin) < 1.0e-14:
        colors = np.zeros_like(field)
    else:
        # Matplotlib surface colors expect values between 0 and 1.
        colors = (field - fmin) / (fmax - fmin)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(
        x,
        y,
        z,
        facecolors=cm.coolwarm(colors),
        linewidth=0,
        antialiased=True,
        shade=False,
    )

    ax.set_title(f"{quantity} ({component})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_box_aspect([1, 1, 1])

    mappable = cm.ScalarMappable(cmap=cm.coolwarm)
    mappable.set_array(field)
    fig.colorbar(mappable, ax=ax, shrink=0.75)

    plt.tight_layout()
    plt.show()

def plot_surface_meridian(
    quantity: str = "phi_sca",
    component: str = "real",
    k: float = 1.0,
    radius: float = 1.0,
    n_terms: int = 50,
) -> None:
    """
    Plot analytical solution along a meridian.

    phi = 0
    theta = 0 ... pi

    This is often much more informative than a 3D sphere plot.
    """
    theta = np.linspace(0.0, np.pi, 1000)

    points = np.column_stack(
        [
            radius * np.cos(theta),
            radius * np.sin(theta),
            np.zeros_like(theta),
        ]
    )

    if quantity == "phi_inc":
        values = incident_plane_wave(points, k=k)

    elif quantity == "phi_sca":
        values = rigid_sphere_scattered_field(
            points,
            k=k,
            radius=radius,
            n_terms=n_terms,
        )

    elif quantity == "phi_total":
        values = rigid_sphere_total_field(
            points,
            k=k,
            radius=radius,
            n_terms=n_terms,
        )

    else:
        raise ValueError("invalid quantity")

    print()
    print("Meridian diagnostics")
    print("--------------------")
    print(f"max |value| = {np.max(np.abs(values)):.6e}")
    print(f"min Re      = {np.min(np.real(values)):.6e}")
    print(f"max Re      = {np.max(np.real(values)):.6e}")
    print(f"min Im      = {np.min(np.imag(values)):.6e}")
    print(f"max Im      = {np.max(np.imag(values)):.6e}")
    print()

    if component == "real":
        field = np.real(values)

    elif component == "imag":
        field = np.imag(values)

    elif component == "abs":
        field = np.abs(values)

    elif component == "phase":
        field = np.angle(values)

    else:
        raise ValueError("invalid component")

    plt.figure(figsize=(8, 5))
    plt.plot(theta, field)

    plt.xlabel("theta [rad]")
    plt.ylabel(f"{quantity} ({component})")
    plt.title(f"{quantity} ({component}) on sphere surface")

    plt.grid(True)
    plt.tight_layout()
    plt.show()


def spherical_hankel1(n: int | Array, z: Array) -> Array:
    """
    Evaluate the spherical Hankel function of the first kind.

    Formula
    -------
        h_n^{(1)}(z) = j_n(z) + i y_n(z),

    where j_n is the spherical Bessel function and y_n is the spherical Neumann
    function.
    """
    return spherical_jn(n, z) + 1j * spherical_yn(n, z)


def spherical_hankel1_derivative(n: int | Array, z: Array) -> Array:
    """
    Evaluate derivative of the spherical Hankel function of the first kind.

    Formula
    -------
        d h_n^{(1)} / dz = d j_n / dz + i d y_n / dz.
    """
    return spherical_jn(n, z, derivative=True) + 1j * spherical_yn(n, z, derivative=True)


def _legendre_derivative(n: int, mu: Array) -> Array:
    """Evaluate dP_n(mu)/dmu, including the two endpoint limits."""

    values = np.asarray(mu, dtype=float)
    if n == 0:
        return np.zeros_like(values)

    derivative = np.empty_like(values)
    interior = np.abs(values) < 1.0 - 1.0e-12
    pn = eval_legendre(n, values[interior])
    pn_minus_one = eval_legendre(n - 1, values[interior])
    derivative[interior] = n * (pn_minus_one - values[interior] * pn) / (
        1.0 - values[interior] ** 2
    )

    positive_endpoint = values >= 1.0 - 1.0e-12
    negative_endpoint = values <= -1.0 + 1.0e-12
    endpoint_magnitude = 0.5 * n * (n + 1)
    derivative[positive_endpoint] = endpoint_magnitude
    derivative[negative_endpoint] = ((-1) ** (n + 1)) * endpoint_magnitude
    return derivative


def rigid_sphere_coefficients(
    k: float,
    radius: float = 1.0,
    n_terms: int = 40,
) -> Array:
    """
    Compute scattered-field expansion coefficients for a rigid sphere.

    Parameters
    ----------
    k:
        Wavenumber.

    radius:
        Sphere radius.

    n_terms:
        Number of partial-wave terms.  Terms n = 0, ..., n_terms - 1 are used.

    Returns
    -------
    Array
        Complex coefficients A_n with shape (n_terms,).

    Notes
    -----
    The scattered field is

        phi_sc = sum A_n h_n^{(1)}(k r) P_n(cos(theta)).
    """
    if k <= 0.0:
        raise ValueError("k must be positive")
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if n_terms <= 0:
        raise ValueError("n_terms must be positive")

    ka = k * radius
    n = np.arange(n_terms)

    # The rigid boundary condition is d(phi_total)/dr = 0 at r = a.
    # Since d/d r of j_n(k r) gives k j_n'(k r), the factor k appears in both
    # incident and scattered terms and cancels from this coefficient formula.
    jn_prime = spherical_jn(n, ka, derivative=True)
    hn_prime = spherical_hankel1_derivative(n, ka)

    return -(2 * n + 1) * (1j ** n) * jn_prime / hn_prime


def robin_sphere_coefficients(
    k: float,
    radius: float = 1.0,
    robin_a: complex = 1.0 + 0.0j,
    robin_b: complex = 0.5 + 0.0j,
    n_terms: int = 40,
) -> Array:
    """
    Compute scattered-field coefficients for a Robin sphere.

    The boundary condition is imposed on the total field at r = radius:

        robin_a * phi_total + robin_b * d(phi_total)/dr = 0.

    The derivative is taken along the outward radial normal of the sphere.
    """
    if k <= 0.0:
        raise ValueError("k must be positive")
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if n_terms <= 0:
        raise ValueError("n_terms must be positive")

    ka = k * radius
    n = np.arange(n_terms)
    # incident_scale is the coefficient multiplying j_n(k r) P_n(mu) in the
    # plane-wave expansion exp(i k r mu).
    incident_scale = (2 * n + 1) * (1j ** n)

    jn = spherical_jn(n, ka)
    jn_prime = spherical_jn(n, ka, derivative=True)
    hn = spherical_hankel1(n, ka)
    hn_prime = spherical_hankel1_derivative(n, ka)

    # The Robin condition is a*phi + b*dphi/dr = 0.
    # The radial derivative introduces k because z = k r, so
    # d j_n(k r)/dr = k*j_n'(k r).
    numerator = robin_a * jn + robin_b * k * jn_prime
    denominator = robin_a * hn + robin_b * k * hn_prime

    if np.any(np.abs(denominator) <= 1.0e-14):
        raise ValueError("Robin coefficient denominator is near zero")

    return -incident_scale * numerator / denominator


def dirichlet_sphere_coefficients(
    k: float,
    radius: float = 1.0,
    n_terms: int = 40,
) -> Array:
    """Compute coefficients for a pressure-release (Dirichlet) sphere.

    The boundary condition is ``phi_total = 0`` at ``r = radius``.  This is the
    exact Robin special case ``robin_a=1`` and ``robin_b=0``.  Reusing the Robin
    formula keeps one source of truth for all linearly related surface data.
    """

    return robin_sphere_coefficients(
        k=k,
        radius=radius,
        robin_a=1.0 + 0.0j,
        robin_b=0.0 + 0.0j,
        n_terms=n_terms,
    )


def plot_scattering_coefficients(
    k: float,
    radius: float = 1.0,
    n_terms: int = 50,
) -> None:
    """
    Plot scattered-field expansion coefficients.

    Useful for identifying which spherical-harmonic modes dominate.
    """
    coeffs = rigid_sphere_coefficients(
        k=k,
        radius=radius,
        n_terms=n_terms,
    )

    n = np.arange(n_terms)

    plt.figure(figsize=(8, 5))
    plt.semilogy(n, np.abs(coeffs), "o-")

    plt.xlabel("Mode n")
    plt.ylabel("|A_n|")
    plt.title(f"Rigid-sphere scattering coefficients (ka={k*radius:.2f})")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

def incident_plane_wave(
    points: Array,
    k: float,
    direction: Array | None = None,
    potential_amplitude: complex = 1.0 + 0.0j,
) -> Array:
    """
    Evaluate the incident plane wave.

    Formula
    -------
        phi_inc(x) = potential_amplitude * exp(i k d · x).

    Parameters
    ----------
    points:
        Point array with shape (3,) or (n_points, 3).

    k:
        Wavenumber.

    direction:
        Incident propagation direction.  If None, (1, 0, 0) is used.

    potential_amplitude:
        Complex velocity-potential amplitude.

    Returns
    -------
    Array
        Complex incident field values.
    """
    pts = _as_points(points, "points")
    d = np.array([1.0, 0.0, 0.0]) if direction is None else np.asarray(direction, dtype=float)
    if d.shape != (3,):
        raise ValueError("direction must have shape (3,)")
    norm = np.linalg.norm(d)
    if norm <= 1.0e-14:
        raise ValueError("direction must be non-zero")
    d = d / norm
    # Matrix-vector product gives d dot x for every point x.
    phase = pts @ d
    return potential_amplitude * np.exp(1j * k * phase)


def rigid_sphere_scattered_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
) -> Array:
    """
    Evaluate the analytical scattered field outside a rigid sphere.

    Parameters
    ----------
    points:
        Field points with shape (3,) or (n_points, 3).  Points should satisfy
        |x| >= radius.

    k:
        Wavenumber.

    radius:
        Sphere radius.

    direction:
        Incident wave direction.  If None, (1, 0, 0) is used.

    n_terms:
        Number of partial-wave terms.

    Returns
    -------
    Array
        Complex scattered field values.

    Raises
    ------
    ValueError
        If any field point is inside the sphere.
    """
    pts = _as_points(points, "points")
    single_point = np.asarray(points).shape == (3,)

    d = np.array([1.0, 0.0, 0.0]) if direction is None else np.asarray(direction, dtype=float)
    if d.shape != (3,):
        raise ValueError("direction must have shape (3,)")
    d = d / np.linalg.norm(d)

    r = np.linalg.norm(pts, axis=1)
    if np.any(r < radius * (1.0 - 1.0e-12)):
        raise ValueError("analytical exterior solution requested inside the sphere")

    # mu = cos(theta), where theta is measured from the incident-wave direction.
    mu = (pts @ d) / r
    mu = np.clip(mu, -1.0, 1.0)

    coeffs = rigid_sphere_coefficients(k=k, radius=radius, n_terms=n_terms)
    kr = k * r

    # Sum the partial-wave expansion term by term.
    phi = np.zeros(pts.shape[0], dtype=complex)
    for n in range(n_terms):
        phi += coeffs[n] * spherical_hankel1(n, kr) * eval_legendre(n, mu)

    if single_point:
        return phi[0]
    return phi


def robin_sphere_scattered_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    robin_a: complex = 1.0 + 0.0j,
    robin_b: complex = 0.5 + 0.0j,
    n_terms: int = 40,
) -> Array:
    """
    Evaluate the analytical scattered field outside a Robin sphere.
    """
    pts = _as_points(points, "points")
    single_point = np.asarray(points).shape == (3,)

    d = np.array([1.0, 0.0, 0.0]) if direction is None else np.asarray(direction, dtype=float)
    if d.shape != (3,):
        raise ValueError("direction must have shape (3,)")
    d = d / np.linalg.norm(d)

    r = np.linalg.norm(pts, axis=1)
    if np.any(r < radius * (1.0 - 1.0e-12)):
        raise ValueError("analytical exterior solution requested inside the sphere")

    # mu = cos(theta), where theta is measured from the incident-wave direction.
    mu = (pts @ d) / r
    mu = np.clip(mu, -1.0, 1.0)

    coeffs = robin_sphere_coefficients(
        k=k,
        radius=radius,
        robin_a=robin_a,
        robin_b=robin_b,
        n_terms=n_terms,
    )
    kr = k * r

    # The field shape is the same as the rigid case; only the coefficients
    # change because the boundary condition has changed.
    phi = np.zeros(pts.shape[0], dtype=complex)
    for n in range(n_terms):
        phi += coeffs[n] * spherical_hankel1(n, kr) * eval_legendre(n, mu)

    if single_point:
        return phi[0]
    return phi


def dirichlet_sphere_scattered_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
) -> Array:
    """Evaluate the scattered field outside a pressure-release sphere."""

    return robin_sphere_scattered_field(
        points=points,
        k=k,
        radius=radius,
        direction=direction,
        robin_a=1.0 + 0.0j,
        robin_b=0.0 + 0.0j,
        n_terms=n_terms,
    )


def dirichlet_sphere_scattered_normal_derivative(
    points: Array,
    normals: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
) -> Array:
    """Evaluate the scattered-field derivative along supplied unit normals.

    The surface solver differentiates along its stored nodal normal.  Those
    normals are almost radial on a spherical mesh but need not be exactly so.
    This routine evaluates the full spherical-coordinate gradient and projects
    it onto the same supplied normals, avoiding a hidden normal-direction error
    in the Dirichlet validation.
    """

    pts = _as_points(points, "points")
    single_point = np.asarray(points).shape == (3,)
    normal_array = _as_points(normals, "normals")
    if normal_array.shape != pts.shape:
        raise ValueError("normals must have the same shape as points")

    normal_lengths = np.linalg.norm(normal_array, axis=1)
    if np.any(normal_lengths <= 1.0e-14):
        raise ValueError("normals must be non-zero")
    unit_normals = normal_array / normal_lengths[:, None]

    d = np.array([1.0, 0.0, 0.0]) if direction is None else np.asarray(direction, dtype=float)
    if d.shape != (3,):
        raise ValueError("direction must have shape (3,)")
    d = d / np.linalg.norm(d)

    r = np.linalg.norm(pts, axis=1)
    if np.any(r < radius * (1.0 - 1.0e-12)):
        raise ValueError("analytical exterior solution requested inside the sphere")
    radial_direction = pts / r[:, None]
    mu = np.clip(radial_direction @ d, -1.0, 1.0)

    # grad(mu) = (d - mu*e_r)/r.  The two projections below convert the
    # spherical radial and angular derivatives to the stored mesh normal.
    radial_projection = np.sum(radial_direction * unit_normals, axis=1)
    mu_gradient_projection = np.sum(
        (d[None, :] - mu[:, None] * radial_direction) * unit_normals,
        axis=1,
    ) / r

    coefficients = dirichlet_sphere_coefficients(k=k, radius=radius, n_terms=n_terms)
    kr = k * r
    derivative = np.zeros(pts.shape[0], dtype=complex)
    for n in range(n_terms):
        pn = eval_legendre(n, mu)
        dpn_dmu = _legendre_derivative(n, mu)
        hn = spherical_hankel1(n, kr)
        dhn_dkr = spherical_hankel1_derivative(n, kr)
        derivative += coefficients[n] * (
            k * dhn_dkr * pn * radial_projection
            + hn * dpn_dmu * mu_gradient_projection
        )

    if single_point:
        return derivative[0]
    return derivative


def rigid_sphere_total_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
) -> Array:
    """
    Evaluate the analytical total field outside a rigid sphere.

    Formula
    -------
        phi_total = phi_incident + phi_scattered.
    """
    return incident_plane_wave(points, k=k, direction=direction) + rigid_sphere_scattered_field(
        points=points,
        k=k,
        radius=radius,
        direction=direction,
        n_terms=n_terms,
    )


def robin_sphere_total_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    robin_a: complex = 1.0 + 0.0j,
    robin_b: complex = 0.5 + 0.0j,
    n_terms: int = 40,
) -> Array:
    """
    Evaluate the analytical total field outside a Robin sphere.
    """
    return incident_plane_wave(points, k=k, direction=direction) + robin_sphere_scattered_field(
        points=points,
        k=k,
        radius=radius,
        direction=direction,
        robin_a=robin_a,
        robin_b=robin_b,
        n_terms=n_terms,
    )


def dirichlet_sphere_total_field(
    points: Array,
    k: float,
    radius: float = 1.0,
    direction: Array | None = None,
    n_terms: int = 40,
) -> Array:
    """Evaluate the total field outside a pressure-release sphere."""

    return (
        incident_plane_wave(points, k=k, direction=direction)
        + dirichlet_sphere_scattered_field(
            points=points,
            k=k,
            radius=radius,
            direction=direction,
            n_terms=n_terms,
        )
    )


def estimate_required_terms(k: float, radius: float = 1.0, safety: int = 20) -> int:
    """
    Estimate a reasonable number of partial-wave terms.

    This simple rule is enough for small and moderate k a values:

        n_terms ≈ k a + safety.

    Parameters
    ----------
    k:
        Wavenumber.

    radius:
        Sphere radius.

    safety:
        Extra number of terms.

    Returns
    -------
    int
        Suggested number of terms.
    """
    # A larger ka needs more spherical-harmonic modes.  The safety margin is
    # deliberately conservative for validation, where the analytical reference
    # should be much more accurate than the numerical BEM result.
    return max(10, int(np.ceil(k * radius + safety)))


def convergence_check(
    k: float = 1.0,
    radius: float = 1.0,
) -> None:
    """
    Check convergence of the partial-wave expansion.
    """
    point = np.array([2.0 * radius, 0.0, 0.0])

    terms = [10, 20, 40, 60, 80]

    values = []

    for n_terms in terms:
        # Recompute the same field with more and more partial-wave terms.  The
        # change between successive values should decrease when the expansion is
        # converged.
        value = rigid_sphere_scattered_field(
            point,
            k=k,
            radius=radius,
            n_terms=n_terms,
        )

        values.append(value)

        print(
            f"n_terms={n_terms:3d} "
            f"phi={value.real:+.8e} "
            f"{value.imag:+.8e}j"
        )

    print()

    for i in range(1, len(values)):
        diff = abs(values[i] - values[i - 1])

        print(
            f"{terms[i-1]:3d} -> {terms[i]:3d} "
            f"change = {diff:.6e}"
        )

def audit_scattered_field(
    k: float = 1.0,
    radius: float = 1.0,
    n_terms: int = 20,
) -> None:
    """
    Audit the rigid-sphere scattered-field expansion.
    """
    ka = k * radius

    coeffs = rigid_sphere_coefficients(
        k=k,
        radius=radius,
        n_terms=n_terms,
    )

    print()
    print("Scattered-field audit")
    print("---------------------")
    print(f"ka = {ka:.6f}")

    total = 0.0 + 0.0j

    for n in range(n_terms):

        hn = spherical_hankel1(n, np.array([ka]))[0]

        contribution = coeffs[n] * hn

        total += contribution

        print(
            f"n={n:2d} "
            f"|A_n|={abs(coeffs[n]):.6e} "
            f"|h_n|={abs(hn):.6e} "
            f"|A_n h_n|={abs(contribution):.6e}"
        )

    print()
    print(
        f"|sum(A_n h_n)| = "
        f"{abs(total):.6e}"
    )

def _as_points(points: Array, name: str) -> Array:
    """
    Convert input to a two-dimensional point array.

    Valid input shapes are:

    - (3,)
    - (n_points, 3)

    Returns
    -------
    Array
        Array with shape (n_points, 3).
    """
    array = np.asarray(points, dtype=float)
    if array.shape == (3,):
        # A single point is reshaped into one row so the rest of the code can
        # use the same vectorized operations for one point and many points.
        return array.reshape(1, 3)
    if array.ndim == 2 and array.shape[1] == 3:
        return array
    raise ValueError(f"{name} must have shape (3,) or (n_points, 3)")


def _self_check_rigid_boundary_condition() -> None:
    """
    Check the rigid boundary condition by one-sided finite difference outside the sphere.

    The analytical total field should satisfy d(phi_total)/dr = 0 at r = a.
    """
    k = 1.0
    radius = 1.0
    n_terms = 50
    eps = 1.0e-5

    directions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.3, 0.4, 0.5],
        ]
    )
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]

    surface = radius * directions
    outer = (radius + eps) * directions

    phi_surface = rigid_sphere_total_field(surface, k=k, radius=radius, n_terms=n_terms)
    phi_outer = rigid_sphere_total_field(outer, k=k, radius=radius, n_terms=n_terms)

    radial_derivative = (phi_outer - phi_surface) / eps
    error = np.max(np.abs(radial_derivative))

    print("Rigid sphere boundary-condition self-check")
    print(f"  max |d(phi_total)/dr| at r=a = {error:.6e}")
    assert error < 1.0e-4


if __name__ == "__main__":

    _self_check_rigid_boundary_condition()

    ktest = 1.0
    plot_scattering_coefficients(
        k=ktest,
        radius=1.0,
        n_terms=50,
    )

    audit_scattered_field(
        k=ktest,
        radius=1.0,
        n_terms=20,
    )

    plot_surface_meridian(
        quantity="phi_sca",
        component="real",
        k=ktest,
        radius=1.0,
        n_terms=50,
    )

    plot_surface_meridian(
        quantity="phi_sca",
        component="imag",
        k=ktest,
        radius=1.0,
        n_terms=50,
    )

    convergence_check(
        k=ktest,
        radius=1.0,
    )

    print("All analytical self-checks passed.")
