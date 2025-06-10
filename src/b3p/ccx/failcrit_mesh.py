import numpy as np
import pyvista as pv


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
    Vectorized Puck failure criterion for multiple stress states and fracture plane angles.
    (Copied from previous response for completeness.)
    """
    STRESS = np.atleast_2d(STRESS)
    n = STRESS.shape[0]
    e = np.zeros((n, 2))  # [FF_index, IFF_index]
    THETAMAX = np.zeros(n)

    Pi = np.pi
    deg2rad = Pi / 180.0

    # Fiber Failure
    A = STRESS[:, 0] - (ANU12 - ANU12f * MGF * (E11 / E11f)) * (
        STRESS[:, 1] + STRESS[:, 2]
    )
    e[:, 0] = np.where(A >= 0, A / Xt, np.abs(A) / Xc)

    # Inter-Fiber Failure
    Rn = Yt
    Rn1 = S12
    Rnt = Yc / (2.0 * np.tan(deg2rad * THETAF))
    Pnt = -1.0 / (2.0 * np.tan(2.0 * deg2rad * THETAF))
    Pn1 = Pnt * (Rn1 / Rnt)

    J = np.arange(181)
    THETA = -Pi / 2.0 + J * deg2rad
    cos_theta = np.cos(THETA)
    sin_theta = np.sin(THETA)
    cos_theta2 = cos_theta**2
    sin_theta2 = sin_theta**2
    sin_cos_theta = sin_theta * cos_theta

    SFP = (
        STRESS[:, 1, np.newaxis] * cos_theta2
        + STRESS[:, 2, np.newaxis] * sin_theta2
        + 2.0 * STRESS[:, 5, np.newaxis] * sin_cos_theta
    )
    TNT = (
        -STRESS[:, 1, np.newaxis] * sin_cos_theta
        + STRESS[:, 2, np.newaxis] * sin_cos_theta
        + STRESS[:, 5, np.newaxis] * (cos_theta2 - sin_theta2)
    )
    TN1 = STRESS[:, 3, np.newaxis] * cos_theta + STRESS[:, 4, np.newaxis] * sin_theta

    IFF_tensile = np.sqrt(
        (SFP / Rn) ** 2
        + (TN1 / (Rn1 - Pn1 * SFP)) ** 2
        + (TNT / (Rnt - Pnt * SFP)) ** 2
    )
    IFF_compressive = np.sqrt(
        (TN1 / (Rn1 - Pn1 * SFP)) ** 2 + (TNT / (Rnt - Pnt * SFP)) ** 2
    )
    IFF = np.where(SFP >= 0, IFF_tensile, IFF_compressive)

    e[:, 1] = np.max(IFF, axis=1)
    max_indices = np.argmax(IFF, axis=1)
    THETAMAX = THETA[max_indices] / deg2rad

    return e, THETAMAX


def build_stiffness_matrix(E11, E22, E33, G12, G13, G23, nu12, nu13, nu23):
    """
    Build 6x6 stiffness matrix [C] for an orthotropic material in material coordinates.
    """
    # Compliance matrix components
    S11 = 1 / E11
    S22 = 1 / E22
    S33 = 1 / E33
    S12 = -nu12 / E11
    S13 = -nu13 / E11
    S23 = -nu23 / E22
    S44 = 1 / G23
    S55 = 1 / G13
    S66 = 1 / G12

    # Compliance matrix
    S = np.array(
        [
            [S11, S12, S13, 0, 0, 0],
            [S12, S22, S23, 0, 0, 0],
            [S13, S23, S33, 0, 0, 0],
            [0, 0, 0, S44, 0, 0],
            [0, 0, 0, 0, S55, 0],
            [0, 0, 0, 0, 0, S66],
        ]
    )

    # Stiffness matrix: C = S^-1
    C = np.linalg.inv(S)
    return C


def strain_to_stress(strain, C):
    """
    Convert strain vector [εxx, εyy, εzz, εxy, εyz, εxz] to stress [σxx, σyy, σzz, σxy, σyz, σxz].
    """
    # Strain in Voigt notation: [εxx, εyy, εzz, γxy, γyz, γxz]
    # γij = 2εij for shear components
    strain_voigt = strain.copy()
    strain_voigt[:, 3:] *= 2  # Convert εxy, εyz, εxz to γxy, γyz, γxz
    # Stress = C * strain
    stress = strain_voigt @ C.T  # Shape (n, 6)
    return stress


def transform_stress_to_material(stress, theta):
    """
    Transform stress [σxx, σyy, σzz, σxy, σyz, σxz] to material coordinates (1,2,3)
    where 1-direction is at angle theta (degrees) to x-axis in x-y plane.
    """
    c = np.cos(np.radians(theta))
    s = np.sin(np.radians(theta))
    c2 = c**2
    s2 = s**2
    cs = c * s

    # Transformation matrix for stress (Voigt notation)
    T = np.array(
        [
            [c2, s2, 0, 2 * cs, 0, 0],
            [s2, c2, 0, -2 * cs, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [-cs, cs, 0, c2 - s2, 0, 0],
            [0, 0, 0, 0, c, -s],
            [0, 0, 0, 0, s, c],
        ]
    )

    # Transform stress: σ_material = T * σ_global
    stress_material = stress @ T.T  # Shape (n, 6)
    return stress_material


def compute_laminate_failure(vtp_file, ply_stack, output_vtp="output_failure.vtp"):
    """
    Compute Puck failure indices for a 3D mesh with strains and a stack of plies.

    Parameters:
    - vtp_file: Path to input VTP file with strain data [εxx, εyy, εzz, εxy, εyz, εxz]
    - ply_stack: List of dicts, each with 'theta' (degrees), 'material' (dict of properties)
    - output_vtp: Path to save output VTP file with failure indices
    """
    # Read VTP file
    mesh = pv.read(vtp_file)
    strainid = "TOSTRAIN_1.000"
    if strainid not in mesh.point_data:
        raise ValueError(
            f"VTP file must contain point data array '{strainid}' with 6 components"
        )
    strains = mesh.point_data[strainid]  # Shape (n_points, 6)
    n_points = strains.shape[0]

    # Initialize output arrays
    failure_data = {}

    # Process each ply
    for ply_idx, ply in enumerate(ply_stack):
        theta = ply["theta"]
        mat = ply["material"]

        # Material properties
        E11, E22, E33 = mat["E11"], mat["E22"], mat["E33"]
        G12, G13, G23 = mat["G12"], mat["G13"], mat["G23"]
        nu12, nu13, nu23 = mat["nu12"], mat["nu13"], mat["nu23"]
        Xt, Xc = mat["Xt"], mat["Xc"]
        Yt, Yc = mat["Yt"], mat["Yc"]
        Zt, Zc = mat["Zt"], mat["Zc"]
        S12, S13, S23 = mat["S12"], mat["S13"], mat["S23"]
        THETAF = mat["THETAF"]
        MGF = mat["MGF"]
        ANU12, ANU12f = mat["ANU12"], mat["ANU12f"]
        E11_puck, E11f = mat["E11_puck"], mat["E11f"]

        # Build stiffness matrix
        C = build_stiffness_matrix(E11, E22, E33, G12, G13, G23, nu12, nu13, nu23)

        # Convert strains to stresses in global coordinates
        stress_global = strain_to_stress(strains, C)

        # Transform stresses to material coordinates
        stress_material = transform_stress_to_material(stress_global, theta)

        # Compute Puck failure indices
        UPSTRAN = np.zeros_like(stress_material)  # Placeholder (unused)
        e, THETAMAX = failure_calc(
            stress_material,
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
            E11_puck,
            E11f,
        )

        # Store results
        failure_data[f"FF_ply{ply_idx + 1}"] = e[:, 0]
        failure_data[f"IFF_ply{ply_idx + 1}"] = e[:, 1]
        failure_data[f"THETAMAX_ply{ply_idx + 1}"] = THETAMAX

        # Print summary
        print(f"Ply {ply_idx + 1} (θ={theta}°):")
        print(f"  Max FF Index: {np.max(e[:, 0]):.3f}")
        print(f"  Max IFF Index: {np.max(e[:, 1]):.3f}")
        print(f"  Failure Points: {np.sum(np.any(e >= 1.0, axis=1))}/{n_points}")

    # Save results to VTP file
    for key, data in failure_data.items():
        mesh.point_data[key] = data
    mesh.save(output_vtp)
    print(f"Results saved to {output_vtp}")


# Example usage
if __name__ == "__main__":
    # Example ply stack: 3 plies with different angles and materials
    ply_stack = [
        {
            "theta": 0.0,  # Fiber angle (degrees)
            "material": {
                # Stiffness properties
                "E11": 150000.0,
                "E22": 10000.0,
                "E33": 10000.0,  # MPa
                "G12": 5000.0,
                "G13": 5000.0,
                "G23": 4000.0,  # MPa
                "nu12": 0.3,
                "nu13": 0.3,
                "nu23": 0.4,
                # Strength properties
                "Xt": 2000.0,
                "Xc": 1500.0,  # MPa
                "Yt": 50.0,
                "Yc": 200.0,
                "Zt": 50.0,
                "Zc": 200.0,
                "S12": 80.0,
                "S13": 80.0,
                "S23": 80.0,
                # Puck parameters
                "THETAF": 53.0,  # degrees
                "MGF": 0.5,
                "ANU12": 0.3,
                "ANU12f": 0.25,
                "E11_puck": 150000.0,
                "E11f": 200000.0,
            },
        },
        {
            "theta": 45.0,
            "material": {  # Same material for simplicity
                "E11": 150000.0,
                "E22": 10000.0,
                "E33": 10000.0,
                "G12": 5000.0,
                "G13": 5000.0,
                "G23": 4000.0,
                "nu12": 0.3,
                "nu13": 0.3,
                "nu23": 0.4,
                "Xt": 2000.0,
                "Xc": 1500.0,
                "Yt": 50.0,
                "Yc": 200.0,
                "Zt": 50.0,
                "Zc": 200.0,
                "S12": 80.0,
                "S13": 80.0,
                "S23": 80.0,
                "THETAF": 53.0,
                "MGF": 0.5,
                "ANU12": 0.3,
                "ANU12f": 0.25,
                "E11_puck": 150000.0,
                "E11f": 200000.0,
            },
        },
        {
            "theta": 90.0,
            "material": {  # Different material
                "E11": 120000.0,
                "E22": 8000.0,
                "E33": 8000.0,
                "G12": 4000.0,
                "G13": 4000.0,
                "G23": 3000.0,
                "nu12": 0.28,
                "nu13": 0.28,
                "nu23": 0.35,
                "Xt": 1800.0,
                "Xc": 1300.0,
                "Yt": 40.0,
                "Yc": 180.0,
                "Zt": 40.0,
                "Zc": 180.0,
                "S12": 70.0,
                "S13": 70.0,
                "S23": 70.0,
                "THETAF": 50.0,
                "MGF": 0.4,
                "ANU12": 0.28,
                "ANU12f": 0.23,
                "E11_puck": 120000.0,
                "E11f": 180000.0,
            },
        },
    ]

    # Input VTP file (replace with your file path)
    vtp_file = "__temp.vtp"

    # Compute failure indices
    compute_laminate_failure(vtp_file, ply_stack, output_vtp="output_failure.vtp")
