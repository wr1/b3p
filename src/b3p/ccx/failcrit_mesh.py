import numpy as np
import pyvista as pv
import logging
import os

logger = logging.getLogger(__name__)


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
    """Vectorized Puck failure criterion for multiple stress states and fracture plane angles."""
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
    """Build 6x6 stiffness matrix [C] for an orthotropic material in material coordinates."""
    S11 = 1 / E11
    S22 = 1 / E22
    S33 = 1 / E33
    S12 = -nu12 / E11
    S13 = -nu13 / E11
    S23 = -nu23 / E22
    S44 = 1 / G23
    S55 = 1 / G13
    S66 = 1 / G12

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

    C = np.linalg.inv(S)
    return C


def strain_to_stress(strain, C):
    """Convert strain vector [εxx, εyy, εzz, εxy, εyz, εxz] to stress [σxx, σyy, σzz, σxy, σyz, σxz]."""
    strain_voigt = strain.copy()
    strain_voigt[:, 3:] *= 2
    stress = strain_voigt @ C.T
    return stress


def transform_stress_to_material(stress, theta):
    """Transform stress [σxx, σyy, σzz, σxy, σyz, σxz] to material coordinates at angle theta (degrees)."""
    c = np.cos(np.radians(theta))
    s = np.sin(np.radians(theta))
    c2 = c**2
    s2 = s**2
    cs = c * s

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

    stress_material = stress @ T.T
    return stress_material


def get_ply_stack():
    """Return default ply stack configuration with angles relative to local_x."""
    return [
        {
            "theta": 0.0,  # 0° relative to local_x
            "material": {
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
            "theta": 45.0,  # 45° relative to local_x
            "material": {
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
            "theta": 90.0,  # 90° relative to local_x
            "material": {
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


def compute_laminate_failure(mesh, ply_stack, alignment_threshold=0.9):
    """Compute Puck failure indices for a surface mesh with cell-based strains and local ply angles."""
    # Ensure cell normals and strains
    mesh = mesh.compute_normals(point_normals=False, cell_normals=True)
    strainid = "TOSTRAIN_1.000"
    if strainid not in mesh.point_data:
        raise ValueError(
            f"Mesh must contain point data array '{strainid}' with 6 components"
        )

    # Interpolate nodal strains to cell centers
    mesh = mesh.point_data_to_cell_data()
    strains = mesh.cell_data[strainid]
    n_cells = strains.shape[0]

    # Compute local x and y vectors
    normals = mesh.cell_data["Normals"]
    global_z = np.array([0, 0, 1])
    dot_z_n = np.sum(normals * global_z, axis=1, keepdims=True)
    local_x = global_z - dot_z_n * normals
    local_x_norm = np.linalg.norm(local_x, axis=1, keepdims=True)
    local_x = np.where(local_x_norm > 1e-10, local_x / local_x_norm, 0)
    local_y = np.cross(normals, local_x)
    local_y_norm = np.linalg.norm(local_y, axis=1, keepdims=True)
    local_y = np.where(local_y_norm > 1e-10, local_y / local_y_norm, 0)

    # Filter cells by alignment of local_x with global z
    dot_products = np.sum(local_x * global_z, axis=1)
    keep_mask = dot_products <= alignment_threshold

    if not np.any(keep_mask):
        logger.info("No cells meet the alignment threshold. Output file not created.")
        return None

    indices_to_keep = np.where(keep_mask)[0]
    new_mesh = mesh.extract_cells(indices_to_keep)
    new_strains = new_mesh.cell_data[strainid]
    new_local_x = local_x[keep_mask]
    new_local_y = local_y[keep_mask]
    n_cells = new_strains.shape[0]

    failure_data = {}

    for ply_idx, ply in enumerate(ply_stack):
        theta_global = np.zeros(n_cells)
        for i in range(n_cells):
            # Define ply angle relative to local_x
            local_x_i = new_local_x[i]
            local_y_i = new_local_y[i]
            # Fiber direction: rotate local_x by theta degrees in local x-y plane
            theta_rad = np.radians(ply["theta"])
            fiber_dir = np.cos(theta_rad) * local_x_i + np.sin(theta_rad) * local_y_i
            # Compute global angle of fiber direction relative to global x-axis
            theta_global[i] = np.degrees(np.arctan2(fiber_dir[1], fiber_dir[0]))

        mat = ply["material"]
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

        # Compute stresses for each cell
        stresses_material = np.zeros((n_cells, 6))
        for i in range(n_cells):
            C = build_stiffness_matrix(E11, E22, E33, G12, G13, G23, nu12, nu13, nu23)
            stress_global = strain_to_stress(new_strains[i : i + 1], C)
            stresses_material[i] = transform_stress_to_material(
                stress_global, theta_global[i]
            )

        UPSTRAN = np.zeros_like(stresses_material)
        e, THETAMAX = failure_calc(
            stresses_material,
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

        failure_data[f"FF_ply{ply_idx + 1}"] = e[:, 0]
        failure_data[f"IFF_ply{ply_idx + 1}"] = e[:, 1]
        failure_data[f"THETAMAX_ply{ply_idx + 1}"] = THETAMAX

        logger.info(f"Ply {ply_idx + 1} (θ={ply['theta']}° relative to local_x):")
        logger.info(f"  Max FF Index: {np.max(e[:, 0]):.3f}")
        logger.info(f"  Max IFF Index: {np.max(e[:, 1]):.3f}")
        logger.info(f"  Failure Cells: {np.sum(np.any(e >= 1.0, axis=1))}/{n_cells}")

    # Add local vectors to mesh
    new_mesh.cell_data["Local_X"] = new_local_x
    new_mesh.cell_data["Local_Y"] = new_local_y

    return new_mesh, failure_data


def compute_failure_for_meshes(mesh_files):
    """Process multiple mesh files and compute failure indices for each."""
    for mesh_file in mesh_files:
        try:
            mesh = pv.read(mesh_file).extract_surface()
            ply_stack = get_ply_stack()
            result = compute_laminate_failure(mesh, ply_stack, 0.9)

            if result is None:
                logger.info(f"No valid cells for mesh {mesh_file}; skipping.")
                continue

            new_mesh, failure_data = result
            for key, data in failure_data.items():
                new_mesh.cell_data[key] = data
            output_file = os.path.splitext(mesh_file)[0] + "_fail.vtu"
            new_mesh.save(output_file)
            logger.info(f"Written output to {output_file}")
        except Exception as e:
            logger.error(f"Error processing {mesh_file}: {e}")
