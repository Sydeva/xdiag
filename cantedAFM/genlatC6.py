#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 The xdiag Contributors
"""
genlatC6.py -- Generate a hexagonal (honeycomb) lattice torus with C6
rotation symmetry and xdiag-compatible TOML irrep sections.

Torus parametrisation (same (a, b) convention as genlat.py)
------------------------------------------------------------
Bravais torus matrix (columns = torus vectors in lattice coords):

    T = [[a,  b ],
         [-b, a-b]]

    Ncell = a² - a·b + b²   (Loeschian integer, number of Bravais cells)
    Nsite = 2·Ncell          (honeycomb: two atoms A, B per cell)

Primitive lattice vectors (Cartesian):
    a1 = (1, 0),   a2 = (1/2, √3/2)

Honeycomb basis in fractional lattice coordinates (hexagon-centre convention):
    τ_A = (1/3, 1/3)    τ_B = (2/3, 2/3)

C6 rotation (60° CCW) about hexagon centre (origin) in lattice coordinates:
    R6 = [[0, -1], [1, 1]]

    R6 τ_A = τ_B + (-1, 0)    (R6 maps A-site at cell n → B-site at n + (-1, 0))
    R6 τ_B = τ_A + (-1, 1)

    Six-site orbit of A(0,0) under R6:
        A(0,0) → B(-1,0) → A(-1,0) → B(-1,-1) → A(0,-1) → B(0,-1) → A(0,0)
    All six sites lie on a regular hexagon of radius 1/√3 centred at the origin.

Reciprocal / momentum theory
-----------------------------
Reciprocal primitive vectors (a_i·b_j = 2π δ_{ij}):
    b1 = 2π(1, -1/√3),   b2 = 2π(0, 2/√3)

Reduced coordinates q = (q1, q2)  give Cartesian  k = q1·b1 + q2·b2.

Dual torus matrix TinvN^T = [[a-b, b], [-b, a]] maps integer pairs (m1, m2)
to integer reduced reciprocal coordinates  fq = Ncell · q.
Allowed momenta: 0 ≤ fq1, fq2 < Ncell.

R6 dual action on integer reciprocal coordinates:
    R6* : (fq1, fq2) → ((fq1 - fq2) % N, fq1 % N)

C6 character table (ζ₆ = exp(2πi/6), ω = exp(2πi/3)):
-------------------------------------------------------
Irrep    E    C6    C3    C2    C3²   C6⁵    little-group label
  A      1     1     1     1     1     1      (trivial)
  B      1    -1     1    -1     1    -1
  E1a    1    ζ₆    ω    -1   -ζ₆  -(ζ₆)*
  E1b    1   ζ₆*   ω*   -1   -(ζ₆)*  ζ₆   (complex conjugate of E1a)
  E2a    1    ω    ω*    1     ω    ω*     (ω = e^{2πi/3})
  E2b    1    ω*   ω     1    ω*    ω     (complex conjugate of E2a)

Character formula:  χ(g) = exp(2πi q·r_l) × ζ₆^(ir_k · k)
where g has rotation index k ∈ {0,...,5} and translation r_l = cells[l].

High-symmetry points and their little groups:
----------------------------------------------
  Γ = (0, 0)         always;   little group C6  →  6 irreps A,B,E1a,E1b,E2a,E2b
  K = (2N/3, N/3)   if N%3=0;  little group C3  →  3 irreps A,Ea,Eb
  K'= (N/3, 2N/3)   if N%3=0;  little group C3  →  3 irreps A,Ea,Eb
  M (3 equiv pts)    if N%2=0;  little group C2  →  2 irreps A,B   (single key)
  generic            C6 orbit;  little group C1  →  1 irrep A      (one key per orbit)

C3 character table at K/K' (ω = e^{2πi/3}, generator = R6²):
  A  : χ(id)=1, χ(R6²)=1,  χ(R6⁴)=1
  Ea : χ(id)=1, χ(R6²)=ω,  χ(R6⁴)=ω²
  Eb : χ(id)=1, χ(R6²)=ω², χ(R6⁴)=ω

C2 character table at M (generator = R6³ = C2):
  A  : χ(id)=1, χ(C2)=1
  B  : χ(id)=1, χ(C2)=-1

Key/naming conventions
-----------------------
- K and K' are labelled separately (distinct momenta).
- The 3 equivalent M points share one key (lexicographic representative).
- Generic C6 orbits use one key per orbit (lexicographic representative).
- Symmetry ordering: block k=0 (pure translations), k=1 (R6∘trans), ..., k=5 (R6⁵∘trans).

Usage
-----
    python genlatC6.py [a [b]] [--toml FILE] [--verify]

    Defaults: a=3, b=0  (Ncell=9, Nsite=18  hexagonal torus)

Examples
--------
    python genlatC6.py 2 1 --toml honeycomb.6.toml   # Nsite=6, K+K' at Γ
    python genlatC6.py 3 0 --toml honeycomb.18.toml  # Nsite=18
    python genlatC6.py 4 2 --toml honeycomb.24.toml  # Nsite=24, K+K'+M
"""

import argparse
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_a1 = np.array([1.0, 0.0])
_a2 = np.array([0.5, np.sqrt(3.0) / 2.0])
_SQRT3 = np.sqrt(3.0)

# 60-degree CCW rotation in lattice coordinates
_R6_lat = np.array([[0, -1], [1, 1]], dtype=int)

# Honeycomb basis (fractional lattice coordinates) around hexagon center.
# Sublattice 0 (A): τ_A = (1/3, 1/3)
# Sublattice 1 (B): τ_B = (2/3, 2/3)
_TAU = {
    0: np.array([1.0 / 3.0, 1.0 / 3.0]),
    1: np.array([2.0 / 3.0, 2.0 / 3.0]),
}
_SUBLABEL = {0: "A", 1: "B"}

# Integer cell shifts produced by R6 acting on basis offsets:
#   R6 τ_A = τ_B + (-1, 0)   →  sub 0 maps to sub 1 in cell + (-1, 0)
#   R6 τ_B = τ_A + (-1, 1)   →  sub 1 maps to sub 0 in cell + (-1, 1)
_ROT_SHIFT = {
    0: np.array([-1, 0], dtype=int),
    1: np.array([-1, 1], dtype=int),
}

# Primitive sixth root of unity and cube root of unity
_ZETA6 = np.exp(2j * np.pi / 6.0)
_OMEGA = np.exp(2j * np.pi / 3.0)

# C6 irreps: (label, exponent k) such that character of R6^n is ζ₆^(k·n)
_C6_IRREPS = [
    ("A",   0),   # trivial
    ("B",   3),   # alternating ±1
    ("E1a", 1),   # ζ₆^n
    ("E1b", 5),   # ζ₆^{-n}  (conj of E1a)
    ("E2a", 2),   # ω^n       (ω = ζ₆²)
    ("E2b", 4),   # ω^{-n}    (conj of E2a)
]

# C3 irreps at K/K': (label, exponent k) such that χ(R6^{2n}) = ω^(k·n)
_C3_IRREPS = [
    ("A",  0),
    ("Ea", 1),
    ("Eb", 2),
]

# C2 irreps at M: (label, parity p) such that χ(R6^{3n}) = (-1)^(p·n)
_C2_IRREPS = [
    ("A", 0),
    ("B", 1),
]


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def lat_to_cart(v):
    """Lattice coordinates (x1, x2), integer or float, to Cartesian vector."""
    return float(v[0]) * _a1 + float(v[1]) * _a2


def site_to_cart(site):
    """Honeycomb site (n1, n2, sub) to Cartesian vector."""
    n1, n2, sub = site
    uv = np.array([n1, n2], dtype=float) + _TAU[int(sub)]
    return lat_to_cart(uv)


def _recip_cart(q1, q2):
    """Reduced reciprocal coordinates (q1, q2) to Cartesian (kx, ky)."""
    kx = 2.0 * np.pi * q1
    ky = 2.0 * np.pi * (-q1 / _SQRT3 + 2.0 * q2 / _SQRT3)
    return kx, ky


# ---------------------------------------------------------------------------
# Torus construction
# ---------------------------------------------------------------------------

def _enumerate_cells(a, b, TinvN, ncell):
    """Enumerate and sort Bravais cells in the torus fundamental domain."""
    rng = abs(a) + abs(b) + 1
    cells = []
    for n2 in range(-rng, rng + 1):
        for n1 in range(-rng, rng + 1):
            f = TinvN @ np.array([n1, n2], dtype=int)
            if 0 <= int(f[0]) < ncell and 0 <= int(f[1]) < ncell:
                cells.append((n1, n2))

    if len(cells) != ncell:
        raise RuntimeError(f"Expected {ncell} Bravais cells, found {len(cells)}.")

    cells.sort(key=lambda v: tuple(int(x) for x in TinvN @ np.array(v, dtype=int)))
    return cells


def build_torus(a, b):
    """
    Build a strict honeycomb torus based on the (a, b) triangular Bravais torus.

    Returns
    -------
    cells  : list[(n1, n2)]               length Ncell
    sites  : list[(n1, n2, sub)]          length Nsite = 2·Ncell
    T      : torus matrix in lattice coordinates
    TinvN  : Ncell · T^{-1} as integer matrix
    Ncell  : number of Bravais cells
    Nsite  : number of honeycomb sites
    """
    ncell = a * a - a * b + b * b
    if ncell <= 0:
        raise ValueError(f"(a={a}, b={b}) gives Ncell = {ncell} <= 0.")

    T = np.array([[a, b], [-b, a - b]], dtype=int)
    TinvN = np.array([[a - b, -b], [b, a]], dtype=int)

    cells = _enumerate_cells(a, b, TinvN, ncell)
    sites = []
    for n1, n2 in cells:
        sites.append((n1, n2, 0))
        sites.append((n1, n2, 1))

    return cells, sites, T, TinvN, ncell, 2 * ncell


# ---------------------------------------------------------------------------
# Permutation helpers
# ---------------------------------------------------------------------------

def _compose(p, q):
    """Compose permutations: result[i] = q[p[i]]."""
    return [q[p[i]] for i in range(len(p))]


def _frac_mod_cell(v, TinvN, ncell):
    """Cell fractional-coordinate tuple of v = (n1, n2), component-wise mod Ncell."""
    f = TinvN @ np.array(v, dtype=int)
    return (int(f[0]) % ncell, int(f[1]) % ncell)


def _build_site_lookup(sites, TinvN, ncell):
    """Map (cell fractional coordinate, sub) to site index."""
    return {
        (_frac_mod_cell((n1, n2), TinvN, ncell), sub): idx
        for idx, (n1, n2, sub) in enumerate(sites)
    }


def compute_translation_perm(sites, TinvN, ncell, delta):
    """Permutation for translation by Bravais vector delta = (dn1, dn2)."""
    lookup = _build_site_lookup(sites, TinvN, ncell)
    perm = []
    for n1, n2, sub in sites:
        key = (_frac_mod_cell((n1 + delta[0], n2 + delta[1]), TinvN, ncell), sub)
        perm.append(lookup[key])
    return perm


def compute_rotation_perm(sites, TinvN, ncell):
    """
    Permutation for the C6 (60° CCW) rotation about the hexagon centre.

    Real-space action:  (n, A) → (R6·n + shift_A, B)
                        (n, B) → (R6·n + shift_B, A)
    where R6 = [[0,-1],[1,1]] and shift_A=(-1,0), shift_B=(-1,1).
    """
    lookup = _build_site_lookup(sites, TinvN, ncell)
    perm = []

    for n1, n2, sub in sites:
        nvec = np.array([n1, n2], dtype=int)
        rot_n = _R6_lat @ nvec + _ROT_SHIFT[sub]
        sub_rot = 1 - sub
        key = (_frac_mod_cell((int(rot_n[0]), int(rot_n[1])), TinvN, ncell), sub_rot)
        perm.append(lookup[key])

    return perm


def compute_nn_bonds(sites, TinvN, ncell):
    """
    Honeycomb nearest-neighbour bonds (A–B only), counted once each: 3·Ncell bonds.

    Each A(n) connects to B(n), B(n − a1), B(n − a2).
    """
    lookup = _build_site_lookup(sites, TinvN, ncell)
    bonds = []
    b_deltas = [(0, 0), (-1, 0), (0, -1)]

    for i, (n1, n2, sub) in enumerate(sites):
        if sub != 0:
            continue
        for d1, d2 in b_deltas:
            key = (_frac_mod_cell((n1 + d1, n2 + d2), TinvN, ncell), 1)
            j = lookup[key]
            bonds.append((i, j))

    return bonds


def generate_all_translations(sites, cells, TinvN, ncell):
    """All Ncell translation permutations, indexed by displacement = cells[l]."""
    lookup = _build_site_lookup(sites, TinvN, ncell)
    all_trans = []

    for d1, d2 in cells:
        perm = []
        for n1, n2, sub in sites:
            key = (_frac_mod_cell((n1 + d1, n2 + d2), TinvN, ncell), sub)
            perm.append(lookup[key])
        all_trans.append(perm)

    return all_trans


# ---------------------------------------------------------------------------
# Symmetry group with metadata
# ---------------------------------------------------------------------------

def _rotation_powers(perm_rot, nsite, order):
    """Return [id, R, R², ..., R^{order-1}] as permutations."""
    powers = [list(range(nsite))]
    for _ in range(1, order):
        powers.append(_compose(perm_rot, powers[-1]))
    return powers


def build_sym_group(sites, cells, TinvN, ncell, nsite, perm_rot):
    """
    Build the full T_{Ncell} × C6 symmetry set with per-element metadata.

    Ordering: block k = 0..5 (rotation power), then all Ncell translations.
    Duplicate permutations (e.g. when C6 acts with order < 6 on a small cluster)
    are silently dropped; first occurrence wins.

    Returns
    -------
    all_syms : list of permutations (each a list of ints)
    sym_meta : list of (k, l)
        k ∈ {0,...,5}   — rotation power index
        l ∈ {0,...,Ncell-1} — translation index (displacement = cells[l])
    """
    all_trans = generate_all_translations(sites, cells, TinvN, ncell)
    rot_powers = _rotation_powers(perm_rot, nsite, 6)

    all_syms = []
    sym_meta = []
    seen = set()

    for k, rot in enumerate(rot_powers):
        for l, trans in enumerate(all_trans):
            key = tuple(_compose(trans, rot))
            if key not in seen:
                seen.add(key)
                all_syms.append(list(key))
                sym_meta.append((k, l))

    return all_syms, sym_meta


# ---------------------------------------------------------------------------
# Momentum helpers
# ---------------------------------------------------------------------------

def _r6_star_int(fq1, fq2, ncell):
    """
    R6* dual action in integer reduced reciprocal coordinates.

    Derivation: R6*q = (R6^{-T})q.  R6^{-1} = [[1,1],[-1,0]] in lattice coords.
    Transposing: (R6^{-T}) maps (q1, q2) → (q1 - q2, q1).
    Integer form: ((fq1 - fq2) % N, fq1 % N).
    """
    return ((fq1 - fq2) % ncell, fq1 % ncell)


def _r6_star_pow(fq, power, ncell):
    """Apply R6*^{power} to integer reduced reciprocal coordinates."""
    out = fq
    for _ in range(power):
        out = _r6_star_int(out[0], out[1], ncell)
    return out


def _orbit_r6(fq, ncell):
    """C6 orbit of fq (unique elements in traversal order)."""
    orbit = []
    cur = fq
    for _ in range(6):
        if cur in orbit:
            break
        orbit.append(cur)
        cur = _r6_star_int(cur[0], cur[1], ncell)
    return orbit


def _is_k_point(fq1, fq2, ncell):
    """True if (fq1, fq2) is the K = (2N/3, N/3) momentum."""
    return (ncell % 3 == 0) and (3 * fq1 == 2 * ncell) and (3 * fq2 == ncell)


def _is_kp_point(fq1, fq2, ncell):
    """True if (fq1, fq2) is the K' = (N/3, 2N/3) momentum."""
    return (ncell % 3 == 0) and (3 * fq1 == ncell) and (3 * fq2 == 2 * ncell)


def _is_m_point(fq1, fq2, ncell):
    """True if (fq1, fq2) is one of the three M points (C2-fixed, ncell even)."""
    if ncell % 2 != 0:
        return False
    h = ncell // 2
    return (fq1, fq2) in {(h, 0), (0, h), (h, h)}


def enumerate_momenta(a, b, ncell):
    """
    Enumerate all Ncell allowed momenta using the dual torus matrix.

    Returns sorted list of (fq1, fq2, q1, q2, kx, ky).
    """
    rng = abs(a) + abs(b) + 1
    result = []
    for m2 in range(-rng, rng + 1):
        for m1 in range(-rng, rng + 1):
            fq1 = (a - b) * m1 + b * m2
            fq2 = -b * m1 + a * m2
            if 0 <= fq1 < ncell and 0 <= fq2 < ncell:
                q1, q2 = fq1 / ncell, fq2 / ncell
                kx, ky = _recip_cart(q1, q2)
                result.append((fq1, fq2, q1, q2, kx, ky))

    assert len(result) == ncell, f"Expected {ncell} momenta, got {len(result)}"
    result.sort(key=lambda x: (x[0], x[1]))
    return result


# ---------------------------------------------------------------------------
# Irreducible representation sections
# ---------------------------------------------------------------------------

def _allowed_indices(sym_meta, allowed_rot_powers):
    """Indices into all_syms whose rotation power is in allowed_rot_powers."""
    return [j for j, (k, _) in enumerate(sym_meta) if k in allowed_rot_powers]


def _phase_for_sym(q1, q2, cell):
    """Translation phase exp(2πi (q1·n1 + q2·n2)) for cell = (n1, n2)."""
    return np.exp(2j * np.pi * (q1 * cell[0] + q2 * cell[1]))


def _clean_char(z, tol=1e-12):
    """Round real and imaginary parts to exact rationals when very close."""
    re = z.real
    im = z.imag
    # snap to common exact values: 0, ±1, ±0.5, ±√3/2
    for exact in (0.0, 1.0, -1.0, 0.5, -0.5, _SQRT3 / 2, -_SQRT3 / 2):
        if abs(re - exact) < tol:
            re = exact
        if abs(im - exact) < tol:
            im = exact
    return complex(re, im)


def _chars_for_section(cells, sym_meta, allowed, q1, q2, rot_char_fn):
    """
    Build the character list for a symmetry section.

    For each symmetry j in `allowed`, the character is:
        χ_j = exp(2πi q·r_{l_j}) × rot_char(k_j)
    where sym_meta[j] = (k_j, l_j) and cells[l_j] = r_{l_j}.
    """
    chars = []
    for j in allowed:
        k, l = sym_meta[j]
        phase = _phase_for_sym(q1, q2, cells[l])
        chars.append(_clean_char(phase * rot_char_fn(k)))
    return chars


def compute_irrep_sections(a, b, ncell, cells, all_syms, sym_meta):
    """
    Compute xdiag-style irrep sections for all momentum orbits.

    Strategy:
    ---------
    Γ = (0,0)       → C6 little group  → 6 sections (A, B, E1a, E1b, E2a, E2b)
    K  = (2N/3,N/3) → C3 little group  → 3 sections (A, Ea, Eb)
    K' = (N/3,2N/3) → C3 little group  → 3 sections (A, Ea, Eb)
    M  (3 equiv pts) → C2 little group  → 2 sections, single M representative
    generic orbits   → C1 little group  → 1 section per C6 orbit (lex representative)

    Returns list of dicts:
        name, characters, allowed_symmetries, momentum [kx, ky]
    """
    del all_syms  # not needed here; characters are built from sym_meta + cells

    momenta = enumerate_momenta(a, b, ncell)
    by_fq = {(fq1, fq2): (q1, q2, kx, ky) for fq1, fq2, q1, q2, kx, ky in momenta}

    sections_gamma = []
    sections_k = []
    sections_kp = []
    sections_m = []
    sections_generic = []
    used = set()

    allowed_c6 = _allowed_indices(sym_meta, {0, 1, 2, 3, 4, 5})
    allowed_c3 = _allowed_indices(sym_meta, {0, 2, 4})
    allowed_c2 = _allowed_indices(sym_meta, {0, 3})
    allowed_c1 = _allowed_indices(sym_meta, {0})

    # ---- first pass: Γ, K, K' ----
    for fq1, fq2, q1, q2, kx, ky in momenta:
        fq = (fq1, fq2)
        if fq in used:
            continue

        if fq == (0, 0):
            for lbl, ir_k in _C6_IRREPS:
                sections_gamma.append({
                    "name": f"Gamma.C6.{lbl}",
                    "characters": _chars_for_section(
                        cells, sym_meta, allowed_c6, q1, q2,
                        rot_char_fn=lambda rpow, k=ir_k: _ZETA6 ** (k * rpow),
                    ),
                    "allowed_symmetries": allowed_c6,
                    "momentum": [kx, ky],
                })
            used.add(fq)
            continue

        if _is_k_point(fq1, fq2, ncell):
            for lbl, ir_k in _C3_IRREPS:
                sections_k.append({
                    "name": f"K.C3.{lbl}",
                    "characters": _chars_for_section(
                        cells, sym_meta, allowed_c3, q1, q2,
                        # C3 generator is R6²; (rpow//2) % 3 gives C3-power index
                        rot_char_fn=lambda rpow, k=ir_k: _OMEGA ** (k * ((rpow // 2) % 3)),
                    ),
                    "allowed_symmetries": allowed_c3,
                    "momentum": [kx, ky],
                })
            used.add(fq)
            continue

        if _is_kp_point(fq1, fq2, ncell):
            for lbl, ir_k in _C3_IRREPS:
                sections_kp.append({
                    "name": f"Kp.C3.{lbl}",
                    "characters": _chars_for_section(
                        cells, sym_meta, allowed_c3, q1, q2,
                        rot_char_fn=lambda rpow, k=ir_k: _OMEGA ** (k * ((rpow // 2) % 3)),
                    ),
                    "allowed_symmetries": allowed_c3,
                    "momentum": [kx, ky],
                })
            used.add(fq)
            continue  # ← was missing before; add it for clarity

    # ---- M points: single representative for all 3 equivalent M momenta ----
    m_points = [(fq1, fq2) for fq1, fq2, *_ in momenta if _is_m_point(fq1, fq2, ncell)]
    if m_points:
        # All 3 M points form one C6 orbit of length 3.
        m_orbit_all = set()
        for m in m_points:
            m_orbit_all.update(_orbit_r6(m, ncell))
        # Use lexicographic minimum as stable representative.
        fq_m = sorted(m for m in m_orbit_all if m in by_fq)[0]
        q1, q2, kx, ky = by_fq[fq_m]

        for lbl, parity in _C2_IRREPS:
            sections_m.append({
                "name": f"M.C2.{lbl}",
                "characters": _chars_for_section(
                    cells, sym_meta, allowed_c2, q1, q2,
                    # C2 generator is R6³; (rpow//3) % 2 gives C2-power index
                    rot_char_fn=lambda rpow, p=parity: (-1.0 + 0j) ** (p * ((rpow // 3) % 2)),
                ),
                "allowed_symmetries": allowed_c2,
                "momentum": [kx, ky],
            })

        for fq in m_orbit_all:
            if fq in by_fq:
                used.add(fq)

    # ---- generic orbits: one section per C6 orbit, lexicographic representative ----
    x_count = 0
    for fq1, fq2, *_ in momenta:
        fq = (fq1, fq2)
        if fq in used:
            continue

        orbit = [o for o in _orbit_r6(fq, ncell) if o in by_fq]
        for o in orbit:
            used.add(o)

        fq_rep = sorted(orbit)[0]
        q1, q2, kx, ky = by_fq[fq_rep]
        name = "X.C1.A" if x_count == 0 else f"X{x_count}.C1.A"

        sections_generic.append({
            "name": name,
            "characters": _chars_for_section(
                cells, sym_meta, allowed_c1, q1, q2,
                rot_char_fn=lambda _k: 1.0 + 0j,
            ),
            "allowed_symmetries": allowed_c1,
            "momentum": [kx, ky],
        })
        x_count += 1

    return sections_gamma + sections_k + sections_kp + sections_m + sections_generic


# ---------------------------------------------------------------------------
# Verification helpers
# ---------------------------------------------------------------------------

def verify_perm(perm, name="perm"):
    """Assert perm is a valid permutation of {0, ..., n-1}."""
    n = len(perm)
    assert sorted(perm) == list(range(n)), f"'{name}' is not a valid permutation of 0..{n-1}."


def verify_c6(perm_rot, nsite):
    """Assert R6^6 = identity on sites."""
    r = list(range(nsite))
    for _ in range(6):
        r = _compose(perm_rot, r)
    if r != list(range(nsite)):
        raise AssertionError("perm_rot does not satisfy R6^6 = id.")


def verify_schur_orthogonality(irrep_sections, tol=1e-8):
    """
    Check Schur orthogonality of characters at each momentum.

    For irreps ρ, σ at the same k-point with the same allowed_symmetries:
        (1/|H|) Σ_g χ_ρ(g) χ_σ(g)* = δ_{ρσ}
    where |H| is the number of allowed symmetries.

    Raises AssertionError if any violation exceeds tol.
    """
    from collections import defaultdict
    by_mom = defaultdict(list)
    for sec in irrep_sections:
        key = (round(sec["momentum"][0], 8), round(sec["momentum"][1], 8))
        by_mom[key].append(sec)

    for mom, secs in by_mom.items():
        if len(secs) < 2:
            continue
        ref_allowed = secs[0]["allowed_symmetries"]
        if not all(s["allowed_symmetries"] == ref_allowed for s in secs):
            continue
        n = len(ref_allowed)
        chi = np.array([[s["characters"][j] for j in range(n)] for s in secs])
        gram = (chi @ chi.conj().T).real / n
        off_diag_err = np.max(np.abs(gram - np.eye(len(secs))))
        assert off_diag_err < tol, (
            f"Schur orthogonality violated at mom={mom}: max off-diag = {off_diag_err:.3e}"
        )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def print_lattice_info(sites, T, ncell, nsite, a, b):
    """Print a human-readable summary of the honeycomb torus."""
    t1c = lat_to_cart(T[:, 0])
    t2c = lat_to_cart(T[:, 1])

    print("# Hexagonal (honeycomb) lattice — C6 torus")
    print(f"# Parameters : a = {a}, b = {b}")
    print(f"# Cells      : Ncell = {ncell}")
    print(f"# Sites      : Nsite = {nsite}")
    print("# Torus vectors:")
    print(f"#   t1 = {T[0,0]}·a1 + ({T[1,0]})·a2  →  Cartesian {t1c}")
    print(f"#   t2 = {T[0,1]}·a1 + ({T[1,1]})·a2  →  Cartesian {t2c}")
    print("# Rotation centre: hexagon centre at origin")
    print("# R6 in lattice coords: [[0,-1],[1,1]];  τ_A=(1/3,1/3), τ_B=(2/3,2/3)")
    print()
    print(f"# {'idx':>4}   {'(n1,n2,sub)':>14}   {'Cartesian (x, y)':>28}")
    print(f"# {'-'*4}   {'-'*14}   {'-'*28}")

    for idx, s in enumerate(sites):
        c = site_to_cart(s)
        print(
            f"  {idx:4d}   ({s[0]:3d},{s[1]:3d},{_SUBLABEL[s[2]]:>1})   "
            f"({c[0]: .10f}, {c[1]: .10f})"
        )


def write_toml(
    filename,
    cells,
    sites,
    T,
    TinvN,
    ncell,
    nsite,
    a,
    b,
    perm_ta1,
    perm_ta2,
    perm_rot,
    coupling_name="J1",
):
    """Write an xdiag-compatible TOML file for the honeycomb C6 torus."""
    t1, t2 = T[:, 0], T[:, 1]

    all_syms, sym_meta = build_sym_group(sites, cells, TinvN, ncell, nsite, perm_rot)
    irrep_sections = compute_irrep_sections(a, b, ncell, cells, all_syms, sym_meta)
    nn_bonds = compute_nn_bonds(sites, TinvN, ncell)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("# Hexagonal (honeycomb) lattice with C6 rotation symmetry\n")
        f.write(f"# Parameters: a={a}, b={b}, Ncell={ncell}, Nsite={nsite}\n")
        f.write(f"# Torus vectors: t1=({t1[0]},{t1[1]}) lat, t2=({t2[0]},{t2[1]}) lat\n")
        f.write("# Generated by genlatC6.py\n\n")

        f.write("Coordinates = [\n")
        for s in sites:
            c = site_to_cart(s)
            f.write(f"  [{c[0]:.15f}, {c[1]:.15f}],\n")
        f.write("]\n\n")

        f.write(f"# Nearest-neighbour Heisenberg interactions ({len(nn_bonds)} bonds)\n")
        f.write("Interactions = [\n")
        for i, j in nn_bonds:
            f.write(f"  ['{coupling_name}', 'SdotS', {i}, {j}],\n")
        f.write("]\n\n")

        f.write("# Translation by primitive vector a1 = (1,0) in lattice coords\n")
        f.write(f"Translation_a1 = {perm_ta1}\n\n")
        f.write("# Translation by primitive vector a2 = (0,1) in lattice coords\n")
        f.write(f"Translation_a2 = {perm_ta2}\n\n")

        f.write("# 60-degree CCW rotation about hexagon centre (C6 generator)\n")
        f.write(f"Rotation_C6 = {perm_rot}\n\n")

        f.write(f"# Full T_{{Ncell}} × C6 symmetry group ({len(all_syms)} elements)\n")
        f.write("Symmetries = [\n")
        for sym in all_syms:
            f.write(f"  {sym},\n")
        f.write("]\n\n")

        f.write("# Irreducible representations\n")
        for sec in irrep_sections:
            f.write(f"\n[{sec['name']}]\n")
            f.write("characters = [\n")
            for chi in sec["characters"]:
                f.write(f"  [{chi.real:.16f}, {chi.imag:.16f}],\n")
            f.write("]\n")
            f.write(f"allowed_symmetries = {sec['allowed_symmetries']}\n")
            kx, ky = sec["momentum"]
            f.write(f"momentum = [{kx:.16f}, {ky:.16f}]\n")

    print(
        f"\n# Wrote: {filename} "
        f"({len(all_syms)} symmetry elements, {len(irrep_sections)} irrep sections)"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("a", nargs="?", type=int, default=3,
                        help="Torus parameter a (default: 3)")
    parser.add_argument("b", nargs="?", type=int, default=0,
                        help="Torus parameter b (default: 0)")
    parser.add_argument("--toml", metavar="FILE", default=None,
                        help="Write xdiag-compatible TOML")
    parser.add_argument("--coupling-name", metavar="NAME", default="J1",
                        help="Coupling label for NN bonds in TOML (default: J1)")
    parser.add_argument("--verify", action="store_true",
                        help="Run Schur-orthogonality self-test after computing irreps")
    args = parser.parse_args()

    a, b = args.a, args.b
    cells, sites, T, TinvN, ncell, nsite = build_torus(a, b)

    perm_ta1 = compute_translation_perm(sites, TinvN, ncell, (1, 0))
    perm_ta2 = compute_translation_perm(sites, TinvN, ncell, (0, 1))
    perm_rot = compute_rotation_perm(sites, TinvN, ncell)

    verify_perm(perm_ta1, "Translation a1")
    verify_perm(perm_ta2, "Translation a2")
    verify_perm(perm_rot, "Rotation C6")
    verify_c6(perm_rot, nsite)

    print_lattice_info(sites, T, ncell, nsite, a, b)
    print("\n# Translation by a1 = (1,0) [lattice]:")
    print(f"  perm_ta1 = {perm_ta1}")
    print("\n# Translation by a2 = (0,1) [lattice]:")
    print(f"  perm_ta2 = {perm_ta2}")
    print("\n# 60-degree CCW rotation (C6):")
    print(f"  perm_rot = {perm_rot}")
    print("\n# Verification: R6^6 = identity OK")

    if args.verify or args.toml:
        all_syms, sym_meta = build_sym_group(sites, cells, TinvN, ncell, nsite, perm_rot)
        irrep_sections = compute_irrep_sections(a, b, ncell, cells, all_syms, sym_meta)
        section_names = [s["name"] for s in irrep_sections]
        print(f"\n# Irrep sections ({len(irrep_sections)}): {section_names}")

        if args.verify:
            verify_schur_orthogonality(irrep_sections)
            print("# Schur orthogonality: OK")

    if args.toml:
        write_toml(
            args.toml,
            cells,
            sites,
            T,
            TinvN,
            ncell,
            nsite,
            a,
            b,
            perm_ta1,
            perm_ta2,
            perm_rot,
            coupling_name=args.coupling_name,
        )

    return cells, sites, T, ncell, nsite, perm_ta1, perm_ta2, perm_rot


if __name__ == "__main__":
    main()

