import warnings

import numpy as np

la = np.linalg


def align_vectors(v1, v2):
    """Euler-Rodrigues rotation of vector 1 to align with vector 2.
    
    Parameters
    ----------
    v1 : `np.ndarray`
        Vector that will be rotated
    v2 : `np.ndarray`
        Vector that we will rotate to (i.e. we will make v1 || to v2)
    
    Returns
    -------
    `np.ndarray`
        3x3 rotation matrix that should be applied to v1.
    
    """

    # Vector we will rotate about
    k = np.cross(v1, v2)
    k /= la.norm(k)

    # Angle we need to rotate
    th = np.arccos(np.dot(v1, v2) / (la.norm(v1) * la.norm(v2)))
    return rotate(k, th)


def rotate(k: np.ndarray, th: float) -> np.ndarray:
    """Rotation around an axis k using Euler-Rodrigues parameterization for rotation.
    
    Parameters
    ----------
    k : np.ndarray
        axis to rotate around
    th : float
        angle of rotation in radians
    
    Returns
    -------
    np.ndarray
        3x3 rotation matrix
    """

    # Make sure k is a unit vector
    k /= la.norm(k)

    # Euler/Rodrigues params
    # See https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
    a = np.cos(th / 2.0)
    b = k[0] * np.sin(th / 2.0)
    c = k[1] * np.sin(th / 2.0)
    d = k[2] * np.sin(th / 2.0)
    r = np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c],
        ]
    )

    return r


def calc_center_of_mass(masses: np.array, coords: np.array) -> np.array:
    """Calculate the center of mass for a given geometry.

    Parameters
    ----------
    masses : np.array
        The masses of the atoms.
    coords : np.array
        The coordinates of molecule.

    Returns
    -------
    np.array
        The coordinates for the center of mass.
    """
    mwc = coords * masses[:, np.newaxis]
    numerator = np.sum(mwc, axis=0)
    denominator = np.sum(masses)

    return numerator / denominator


def make_moi_tensor(masses: np.array, coords: np.array) -> np.array:
    """Calculate the moment of inertia tensor for a given set of masses
    and coordinates.

    Parameters
    ----------
    masses : np.array
        The masses of the atoms.
    coords : np.array
        The coordinates of molecule.

    Returns
    -------
    np.array
        The MOI tensor.
    """
    # https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
    moi_tensor = np.empty((3, 3))

    moi_tensor[0][0] = np.sum(masses * (coords[:, 1] ** 2 + coords[:, 2] ** 2))
    moi_tensor[1][1] = np.sum(masses * (coords[:, 0] ** 2 + coords[:, 2] ** 2))
    moi_tensor[2][2] = np.sum(masses * (coords[:, 0] ** 2 + coords[:, 1] ** 2))

    moi_tensor[0][1] = np.sum(masses * coords[:, 0] * coords[:, 1])
    moi_tensor[0][2] = np.sum(masses * coords[:, 0] * coords[:, 2])
    moi_tensor[1][2] = np.sum(masses * coords[:, 1] * coords[:, 2])

    moi_tensor[1][0] = moi_tensor[0][1]
    moi_tensor[2][0] = moi_tensor[0][2]
    moi_tensor[2][1] = moi_tensor[1][2]
    return moi_tensor


def calc_principal_moi(
    masses: np.array, coords: np.array, tol: float = 1e-8, max_iter: int = 20
) -> np.array:
    """Orient a rigid body and calculate the principal moments of inertia for a molecule.
    To handle linear and other symmetric bodies, we align the smallest principle moment
    with the z-axis.

    Parameters
    ----------
    masses : np.array
        The masses of the atoms.
    coords : np.array
        The coordinates of molecule.
    tol : float, optional
        The tolerance for convergence of the principal moments of inertia, by default 1e-8
    max_iter : int, optional
        The maximum number of iterations for the standardization procedure, by default 20

    Returns
    -------
    np.array
        The principal moments of inertia.
    """

    geom = coords.copy()
    com = calc_center_of_mass(masses, geom)
    geom -= com
    moi_tensor = make_moi_tensor(masses, geom)
    moi_old, v_old = np.linalg.eigh(moi_tensor)

    # If any of the MOI are already 0, return them
    if abs(moi_old.min()) < tol ** 2:
        moi_old[moi_old / moi_old.max() < tol ** 2] = 0
        return moi_old

    # Self Consistent Standardization of Orientation
    moi_diff = 1.0
    z_vec = np.array([0, 0, 1])
    for _ in range(max_iter):
        # Align smallest moment of inertia with z-axis
        v0 = v_old[:, 0]
        rot_mat = align_vectors(v0, z_vec)
        geom = geom.dot(rot_mat)
        moi_tensor_i = make_moi_tensor(masses, geom)
        moi_i, v_old = np.linalg.eigh(moi_tensor_i)

        moi_diff = np.linalg.norm(moi_i - moi_old)

        # Convergence check
        if moi_diff < tol:
            # If any of the values are substantially smaller
            # than the largest principal moment, set them to
            # zero
            print(moi_i)
            moi_i[moi_i / moi_i.max() < tol ** 2] = 0
            return moi_i
        else:
            moi_old = moi_i.copy()

    warnings.warn(
        f"SEE BELOW\nThe reorientation of the moleucle didn't converge (MOI diff={moi_diff})."
        + "\nThis may cause your rotational thermodynamic quantities to have large errors!",
        RuntimeWarning,
    )
    moi_old[moi_old / moi_old.max() < tol ** 2] = 0
    print(moi_old)
    return moi_old

