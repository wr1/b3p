import numpy as np


def failure_calc(
    STRESS,
    UPSTRAN,
    Xt,
    Xc,
    Yt,
    Yc,
    Zt,
    Zc,
    S12,
    S13,
    S23,
    THETAF,
    MGF,
    ANU12,
    ANU12f,
    E11,
    E11f,
):
    """
    Python translation of Fortran subroutine failure_calc, implementing Puck failure criterion.
    Fully vectorized using NumPy for multiple stress states and fracture plane angles.

    Parameters:
    - STRESS: Array of stress tensors, shape (n, 6) [σ1, σ2, σ3, τ23, τ13, τ12] (MPa)
    - UPSTRAN: Array of strains, shape (n, 6) (not used, included for compatibility)
    - Xt, Xc: Tensile/compressive strength in fiber direction (MPa)
    - Yt, Yc: Tensile/compressive strength in transverse direction (y-axis) (MPa)
    - Zt, Zc: Tensile/compressive strength in transverse direction (z-axis) (MPa)
    - S12, S13, S23: Shear strengths (MPa)
    - THETAF: Fracture plane angle (degrees)
    - MGF: Fiber-matrix interaction factor
    - ANU12, ANU12f: Poisson's ratios (matrix and fiber)
    - E11, E11f: Young's moduli in fiber direction (matrix and fiber) (MPa)

    Returns:
    - e: Failure indices, shape (n, 2) [FF_index, IFF_index]
    - THETAMAX: Critical fracture plane angles for IFF, shape (n,)
    """
    # Ensure STRESS is a 2D array with shape (n, 6)
    STRESS = np.atleast_2d(STRESS)
    n = STRESS.shape[0]

    # Initialize output arrays
    e = np.zeros((n, 2))  # [FF_index, IFF_index]
    THETAMAX = np.zeros(n)  # Critical fracture plane angles

    # Constants
    Pi = np.pi
    deg2rad = Pi / 180.0

    # Fiber Failure (FF)
    # A = σ1 - (ν12 - ν12f * MGF * (E11/E11f)) * (σ2 + σ3)
    A = STRESS[:, 0] - (ANU12 - ANU12f * MGF * (E11 / E11f)) * (
        STRESS[:, 1] + STRESS[:, 2]
    )

    # Vectorized FF index calculation
    e[:, 0] = np.where(A >= 0, A / Xt, np.abs(A) / Xc)

    # Inter-Fiber Failure (IFF)
    # Strength parameters
    Rn = Yt
    Rn1 = S12
    Rnt = Yc / (2.0 * np.tan(deg2rad * THETAF))

    # Slope parameters
    Pnt = -1.0 / (2.0 * np.tan(2.0 * deg2rad * THETAF))
    Pn1 = Pnt * (Rn1 / Rnt)

    # Vectorized angle array (-90° to 90°, 181 points)
    J = np.linspace(0, 180, 19)  # np.arange(181)  # 0 to 180
    # print(J)
    THETA = -Pi / 2.0 + J * deg2rad  # Shape (m,), m = 181
    cos_theta = np.cos(THETA)  # Shape (m,)
    sin_theta = np.sin(THETA)  # Shape (m,)
    cos_theta2 = cos_theta**2  # Shape (m,)
    sin_theta2 = sin_theta**2  # Shape (m,)
    sin_cos_theta = sin_theta * cos_theta  # Shape (m,)

    # Broadcast stresses to shape (n, m) for vectorized computation
    # Stresses on fracture plane
    SFP = (
        STRESS[:, 1, np.newaxis] * cos_theta2  # Shape (n, m)
        + STRESS[:, 2, np.newaxis] * sin_theta2
        + 2.0 * STRESS[:, 5, np.newaxis] * sin_cos_theta
    )

    TNT = (
        -STRESS[:, 1, np.newaxis] * sin_cos_theta
        + STRESS[:, 2, np.newaxis] * sin_cos_theta
        + STRESS[:, 5, np.newaxis] * (cos_theta2 - sin_theta2)
    )

    TN1 = STRESS[:, 3, np.newaxis] * cos_theta + STRESS[:, 4, np.newaxis] * sin_theta

    # IFF calculation
    # For SFP >= 0 (tensile normal stress)
    IFF_tensile = np.sqrt(
        (SFP / Rn) ** 2
        + (TN1 / (Rn1 - Pn1 * SFP)) ** 2
        + (TNT / (Rnt - Pnt * SFP)) ** 2
    )
    # Shape (n, m)

    # For SFP < 0 (compressive normal stress)
    IFF_compressive = np.sqrt(
        (TN1 / (Rn1 - Pn1 * SFP)) ** 2 + (TNT / (Rnt - Pnt * SFP)) ** 2
    )
    # Shape (n, m)

    # Combine tensile and compressive cases
    IFF = np.where(SFP >= 0, IFF_tensile, IFF_compressive)  # Shape (n, m)

    # Find maximum IFF and corresponding angle
    e[:, 1] = np.max(IFF, axis=1)  # Shape (n,)
    max_indices = np.argmax(IFF, axis=1)  # Shape (n,)
    THETAMAX = THETA[max_indices] / deg2rad  # Convert back to degrees, shape (n,)

    return e, THETAMAX


# Example usage
if __name__ == "__main__":
    # Material properties (example values)
    Xt = 2000.0  # MPa
    Xc = 1500.0
    Yt = 50.0
    Yc = 200.0
    Zt = 50.0
    Zc = 200.0
    S12 = 80.0
    S13 = 80.0
    S23 = 80.0
    THETAF = 53.0  # degrees
    MGF = 0.5
    ANU12 = 0.3
    ANU12f = 0.25
    E11 = 150000.0  # MPa
    E11f = 200000.0

    # Example stress states [σ1, σ2, σ3, τ23, τ13, τ12] (MPa)
    STRESS = np.array(
        [[1000.0, -30.0, 0.0, 0.0, 0.0, 40.0], [1500.0, 20.0, 10.0, 0.0, 0.0, 50.0]]
    )
    # .repeat(1000000, axis=0)  # Repeat for 10 stress states

    # Strains (not used, included for compatibility)
    UPSTRAN = np.zeros_like(STRESS)

    # Calculate failure indices
    e, THETAMAX = failure_calc(
        STRESS,
        UPSTRAN,
        Xt,
        Xc,
        Yt,
        Yc,
        Zt,
        Zc,
        S12,
        S13,
        S23,
        THETAF,
        MGF,
        ANU12,
        ANU12f,
        E11,
        E11f,
    )

    # Output results
    for i in range(len(STRESS)):
        print(f"Stress State {i+1}:")
        print(f"  Fiber Failure Index: {e[i, 0]:.3f}")
        print(f"  Inter-Fiber Failure Index: {e[i, 1]:.3f}")
        print(f"  Critical Fracture Plane Angle: {THETAMAX[i]:.1f} degrees")
        print(f"  Failure Occurred: {np.any(e[i] >= 1.0)}")
