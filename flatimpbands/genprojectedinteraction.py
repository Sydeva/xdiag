import kwant as kw
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import opt_einsum as oe
from functools import lru_cache

sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])
ident = sigma_z @ sigma_z

rng = np.random.default_rng(42)


def GTAI(L, t2=0.1):
    lat = kw.lattice.square(a=1, norbs=2)
    Ringsym = kw.TranslationalSymmetry([L, 0], [0, L])

    syst = kw.Builder(Ringsym)

    syst[(lat(x, y) for x in range(L) for y in range(L))] = sigma_x

    syst[kw.builder.HoppingKind((1, 0), lat, lat)] = 1 / 2 * (sigma_x + 1j * sigma_y) + 1j * t2 * sigma_z
    syst[kw.builder.HoppingKind((0, 1), lat, lat)] = 1 / 2 * (sigma_x - 1j * sigma_y) + 1j * t2 * sigma_z
    syst[kw.builder.HoppingKind((1, 1), lat, lat)] = -1j * t2 * sigma_z
    syst = kw.wraparound.wraparound(syst)

    return syst, lat


def GTAI3(L, t2=0.1, nhal=0, comb=None):
    lat = kw.lattice.square(a=1, norbs=2)
    Ringsym = kw.TranslationalSymmetry([L, 0], [0, L])

    syst = kw.Builder(Ringsym)

    zero = np.zeros((2, 2), dtype=np.complex128)

    if comb is None:
        comb = np.zeros((L * L))
        comb[rng.choice(L ** 2, nhal, replace=False)] = 1
        comb = comb.reshape((L, L))

    syst[(lat(x, y) for x in range(L) for y in range(L))] = zero

    syst[kw.builder.HoppingKind((0, 1), lat, lat)] = zero
    syst[kw.builder.HoppingKind((1, 1), lat, lat)] = zero
    syst[kw.builder.HoppingKind((1, 0), lat, lat)] = zero

    for x, y in zip(np.argwhere(comb).T[0], np.argwhere(comb).T[1]):
        print(x, y)
        syst[lat(x, (y + 1)), lat(x, y)] = syst[lat(x, (y + 1)), lat(x, y)] + -1j * t2 * (ident - sigma_z) / 2
        syst[lat((x + 1), (y + 1)), lat(x, y)] = syst[lat((x + 1), (y + 1)), lat(x, y)] + 1j * t2 * (
            ident - sigma_z
        ) / 2
        syst[lat(x, (y + 1)), lat((x + 1), (y + 1))] = syst[lat(x, (y + 1)), lat((x + 1), (y + 1))] + 1j * t2 * (
            ident - sigma_z
        ) / 2

        syst[lat((x + 1), y), lat(x, y)] = syst[lat((x + 1), y), lat(x, y)] + 1j * t2 * (ident + sigma_z) / 2
        syst[lat((x + 1), (y + 1)), lat(x, y)] = syst[lat((x + 1), (y + 1)), lat(x, y)] + -1j * t2 * (
            ident + sigma_z
        ) / 2
        syst[lat((x + 1), y), lat((x + 1), (y + 1))] = syst[lat((x + 1), y), lat((x + 1), (y + 1))] + -1j * t2 * (
            ident + sigma_z
        ) / 2

    syst = kw.wraparound.wraparound(syst)

    return syst, lat, comb


def createP2(phix, phiy, syst1, syst2, N, nsub, emin, emax=0, sparse=True):
    H = syst1.hamiltonian_submatrix(params={"k_x": phix, "k_y": phiy}, sparse=sparse) + syst2.hamiltonian_submatrix(
        params={"k_x": phix, "k_y": phiy}, sparse=sparse
    )
    if sparse:
        E, V = sp.sparse.linalg.eigsh(H, k=sparse, sigma=0)
    else:
        E, V = np.linalg.eigh(H)

    V = V[:, (emin < E) & (E < emax)]
    return V


def genPlatt2(syst1, syst2, Nkpoint, nsites, nstates, emin, emax=0):
    Platt = np.zeros((Nkpoint, Nkpoint, nsites, nstates), dtype=np.complex128)

    for i, phix in enumerate((2 * np.pi) / Nkpoint * np.arange(0, Nkpoint, step=1)):
        for j, phiy in enumerate((2 * np.pi) / Nkpoint * np.arange(0, Nkpoint, step=1)):
            Platt[i, j, :, :] = createP2(phix, phiy, syst1, syst2, Nkpoint, nsites, emin, emax=emax, sparse=2 * nstates + 2)
    return Platt


def givebloch(Nkpoint, return_metadata=False):
    # returns array where first 2 indices are for momentum,
    # and third is the wavefunction in orbital basis and fourth the eigen index
    rows = cols = 32
    comb = np.zeros((rows) * (cols // 2))
    comb = comb.reshape((cols // 2, rows))
    comb[rows // 4 - 1, rows // 4 + 1] = 1
    comb[rows // 8 - 1, rows // 2 - 1] = 1
    comb[2 * rows // 5 - 1, 3 * rows // 4 - 3] = 1
    L = 32

    testsys1, lat = GTAI(L, t2=0.0)
    syst1 = testsys1.finalized()

    testsys2, _, comb = GTAI3(L, t2=0.93, nhal=6, comb=comb)
    syst2 = testsys2.finalized()

    bloch = genPlatt2(syst1, syst2, Nkpoint, nsites=2 * 32 ** 2, nstates=1, emin=-0.015)
    if not return_metadata:
        return bloch

    # Wraparound momenta are defined by the translational symmetry vectors [L,0],[0,L].
    trans_vecs = np.array([[L, 0.0], [0.0, L]], dtype=float)
    metadata = {
        "L": L,
        "trans_vecs": trans_vecs,
        "lat_prim_vecs": np.array(lat.prim_vecs, dtype=float),
        "dxs": dxs.copy(),
    }
    return bloch, metadata


subv = np.array([-0.16666667, 0.16666667])
L = 32
nxcos = np.array(list((i // 2) % L + (-1) ** i * subv[0] for i in range(2 * L ** 2))) / L
nycos = np.array(list((i // 2) // L + (-1) ** i * subv[1] for i in range(2 * L ** 2))) / L
dxs = np.array([nxcos, nycos])


def _to_pair(k):
    if isinstance(k, (tuple, list, np.ndarray)) and len(k) == 2:
        return int(k[0]), int(k[1])
    raise ValueError("Momentum must be a 2-tuple/list/array of grid indices.")


def _wrap_pair(k, nk):
    kx, ky = _to_pair(k)
    return kx % nk, ky % nk


def _infer_k4_and_delta_g(k1, k2, k3, nk):
    i1, j1 = _to_pair(k1)
    i2, j2 = _to_pair(k2)
    i3, j3 = _to_pair(k3)
    s0 = i1 + i2 - i3
    s1 = j1 + j2 - j3
    k4 = (s0 % nk, s1 % nk)
    # delta_G indices satisfy: k1+k2-k3-k4 = nk * delta_g
    delta_g = ((s0 - k4[0]) // nk, (s1 - k4[1]) // nk)
    return k4, delta_g


def _reciprocal_vectors(trans_vecs):
    # trans_vecs stores real-space vectors as rows.
    a1 = np.array(trans_vecs[0], dtype=float)
    a2 = np.array(trans_vecs[1], dtype=float)
    amat = np.column_stack((a1, a2))
    bmat = 2.0 * np.pi * np.linalg.inv(amat).T
    b1 = bmat[:, 0]
    b2 = bmat[:, 1]
    return b1, b2


def _magnetic_length_sq(trans_vecs):
    a1 = np.array(trans_vecs[0], dtype=float)
    a2 = np.array(trans_vecs[1], dtype=float)
    omega = abs(np.cross(a1, a2))
    return omega / (2.0 * np.pi)


def _kvec_from_index(k, nk, b1, b2):
    i, j = _wrap_pair(k, nk)
    return (i / nk) * b1 + (j / nk) * b2


def _laguerre_eval(m, x):
    # scipy.special.eval_laguerre accepts scalar x.
    return float(sp.special.eval_laguerre(int(m), float(x)))


def haldane_pseudopotential_q(qvec, pseudo_coeffs, lb2):
    q2 = float(np.dot(qvec, qvec))
    x = 0.5 * q2 * lb2
    envelope = np.exp(-x)
    val = 0.0
    for m, coeff in pseudo_coeffs.items():
        val += coeff * _laguerre_eval(int(m), x)
    return val * envelope


def make_projected_interaction_callable(
    bloch,
    dxs,
    trans_vecs,
    pseudo_coeffs,
    band=0,
    gmax=0,
    antisymmetrized=True,
):
    """
    Build a lazy callable H(k1,k2,k3[,k4]) for projected interaction matrix elements.

    Inputs:
      - bloch[kx, ky, orb, band]
      - dxs[2, orb]: orbital positions in fractional coordinates of the wraparound cell
      - trans_vecs: real-space translation vectors used by kwant wraparound
      - pseudo_coeffs: dict-like {m: V_m} for Haldane pseudopotential channels

    Momentum inputs are integer grid indices (ix, iy).
    If k4 is omitted, momentum conservation is used: k4 = k1 + k2 - k3 (mod Nk).
    By default the returned matrix element is antisymmetrized in both pairs:
      H_AS(1,2;3,4) = H(1,2;3,4) - H(2,1;3,4) - H(1,2;4,3) + H(2,1;4,3).
    """
    if bloch.ndim != 4:
        raise ValueError("bloch must have shape (Nk, Nk, Norb, Nband)")

    nkx, nky, norb, nband = bloch.shape
    if nkx != nky:
        raise ValueError("Only square k-grids are supported: bloch.shape[0] == bloch.shape[1]")
    nk = nkx
    if band < 0 or band >= nband:
        raise ValueError(f"band index {band} out of range [0, {nband})")

    dxs = np.asarray(dxs, dtype=float)
    if dxs.shape != (2, norb):
        raise ValueError(f"dxs must have shape (2, {norb}), got {dxs.shape}")
    pseudo_coeffs = dict(pseudo_coeffs)

    b1, b2 = _reciprocal_vectors(np.asarray(trans_vecs, dtype=float))
    lb2 = _magnetic_length_sq(np.asarray(trans_vecs, dtype=float))
    # convert fractional orbital coordinates to real-space vectors in the wraparound cell
    r_orb = np.einsum("ij,jk->ik", np.asarray(trans_vecs, dtype=float).T, dxs)

    g_list = [(gx, gy) for gx in range(-gmax, gmax + 1) for gy in range(-gmax, gmax + 1)]

    @lru_cache(maxsize=None)
    def _lambda_cached(k_a, k_b, g_idx):
        kax, kay = _wrap_pair(k_a, nk)
        kbx, kby = _wrap_pair(k_b, nk)
        gx, gy = g_idx
        gvec = gx * b1 + gy * b2
        phase = np.exp(-1j * (gvec[0] * r_orb[0] + gvec[1] * r_orb[1]))
        ua = bloch[kax, kay, :, band]
        ub = bloch[kbx, kby, :, band]
        return np.sum(phase * np.conj(ua) * ub)

    def _raw_matrix_element(k1, k2, k3, k4=None):
        k1 = _to_pair(k1)
        k2 = _to_pair(k2)
        k3 = _to_pair(k3)

        if k4 is None:
            k4, delta_g = _infer_k4_and_delta_g(k1, k2, k3, nk)
        else:
            k4 = _to_pair(k4)
            k4_wrapped, delta_g = _infer_k4_and_delta_g(k1, k2, k3, nk)
            if _wrap_pair(k4, nk) != k4_wrapped:
                return 0.0 + 0.0j

        k1vec = _kvec_from_index(k1, nk, b1, b2)
        k4vec = _kvec_from_index(k4, nk, b1, b2)

        total = 0.0 + 0.0j
        for gx, gy in g_list:
            g_idx = (gx, gy)
            gvec = gx * b1 + gy * b2
            qvec = k1vec - k4vec - gvec
            vq = haldane_pseudopotential_q(qvec, pseudo_coeffs, lb2)
            lam1 = _lambda_cached(k1, k4, g_idx)
            g2_idx = (delta_g[0] - gx, delta_g[1] - gy)
            lam2 = _lambda_cached(k2, k3, g2_idx)
            total += vq * lam1 * lam2
        return total

    def element(k1, k2, k3, k4=None):
        if k4 is None:
            k4, _ = _infer_k4_and_delta_g(k1, k2, k3, nk)

        direct = _raw_matrix_element(k1, k2, k3, k4)
        if not antisymmetrized:
            return direct

        swap12 = _raw_matrix_element(k2, k1, k3, k4)
        swap34 = _raw_matrix_element(k1, k2, k4, k3)
        swap12_swap34 = _raw_matrix_element(k2, k1, k4, k3)
        return direct - swap12 - swap34 + swap12_swap34

    # expose useful metadata for downstream scripts
    element.nk = nk
    element.band = band
    element.gmax = gmax
    element.pseudo_coeffs = pseudo_coeffs
    element.lb2 = lb2
    element.trans_vecs = np.asarray(trans_vecs, dtype=float)
    return element


def save_bloch_data(nkpoint, filepath):
    """Generate Bloch wavefunctions for *nkpoint* x *nkpoint* k-grid and save to a .npy file.

    The file stores a dict with keys: 'bloch', 'L', 'trans_vecs', 'lat_prim_vecs', 'dxs'.
    Reload with ``np.load(filepath, allow_pickle=True).item()``.

    Parameters
    ----------
    nkpoint : int
        Number of k-points along each direction.
    filepath : str or path-like
        Destination file path (should end in .npy).

    Returns
    -------
    str
        The filepath that was written.
    """
    bloch, metadata = givebloch(nkpoint, return_metadata=True)
    data = {"bloch": bloch, **metadata}
    np.save(filepath, data, allow_pickle=True)
    return filepath


def make_haldane_interaction_from_givebloch(
    nkpoint,
    pseudo_coeffs,
    band=0,
    gmax=0,
    antisymmetrized=True,
    bloch_file=None,
):
    """Build a projected-interaction callable, optionally loading Bloch data from a file.

    Parameters
    ----------
    nkpoint : int
        Number of k-points along each direction (used only when *bloch_file* is None).
    pseudo_coeffs : dict
        Haldane pseudopotential coefficients ``{m: V_m}``.
    band : int
        Band index within the Bloch array.
    gmax : int
        Maximum reciprocal-lattice shell index for umklapp sums.
    antisymmetrized : bool
        If True, returned matrix elements are antisymmetrized in both
        first and second momentum pairs.
    bloch_file : str or path-like, optional
        Path to a pre-generated ``.npy`` file produced by :func:`save_bloch_data`.
        When provided, Bloch data is loaded from disk instead of being recomputed.
    """
    if bloch_file is not None:
        data = np.load(bloch_file, allow_pickle=True).item()
        bloch = data["bloch"]
        metadata = {k: v for k, v in data.items() if k != "bloch"}
    else:
        bloch, metadata = givebloch(nkpoint, return_metadata=True)
    return make_projected_interaction_callable(
        bloch=bloch,
        dxs=metadata["dxs"],
        trans_vecs=metadata["trans_vecs"],
        pseudo_coeffs=pseudo_coeffs,
        band=band,
        gmax=gmax,
        antisymmetrized=antisymmetrized,
    )

