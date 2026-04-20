#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 The xdiag Contributors
"""
genlat.py -- Generate a triangular (hexagonal) lattice torus that is
             compatible with C3 (threefold) rotation symmetry, label each
             site with an integer, and compute translation / rotation
             symmetry permutations together with irreducible representation
             (irrep) sections at every momentum point on the torus.

Real-space theory
-----------------
A triangular lattice torus with C3 symmetry is parameterised by two integers
(a, b) via the torus matrix

    T = [[a,  b ],
         [-b, a-b]]          (columns = torus vectors in lattice coords)

Number of sites: N = a²−a·b+b² (a Loeschian integer).

Primitive lattice vectors (Cartesian):
    a1 = (1, 0),   a2 = (1/2, √3/2)

120° CCW rotation in lattice coordinates:
    R3 = [[-1,-1],[1,0]]     a1→-a1+a2,  a2→-a1

Site labels 0…N-1 are assigned by sorting (N·T⁻¹·v) lexicographically,
where  N·T⁻¹ = [[a−b,−b],[b,a]].

Reciprocal / momentum theory
----------------------------
Reciprocal primitive vectors (a_i·b_j = 2π δ_{ij}):
    b1 = 2π(1, -1/√3),   b2 = 2π(0, 2/√3)

Reduced coordinates q=(q1,q2) give Cartesian k = q1·b1 + q2·b2.

The dual torus matrix  TinvN^T = [[a-b, b],[-b, a]]  maps integer pairs
(m1,m2) to integer reduced coords fq = N·q.  The N allowed momenta satisfy
0 ≤ fq1, fq2 < N.

R3 acts on momenta (dual action):
    R3* : (q1,q2) → (-q2, q1-q2)  mod 1
    in integer coords: (fq1,fq2) → ((-fq2)%N, (fq1-fq2)%N)

Little group:
    C3   if the momentum is fixed by R3* — occurs at
           Γ = (0,0)          always present
           K = (2/3,1/3)      present iff N % 3 == 0
           K'= (1/3,2/3)      present iff N % 3 == 0
    C1   for all other (generic) momenta

C3 irreps  (ω = e^{2πi/3}):
    A  : χ(id)=1, χ(R3)=1,  χ(R3²)=1
    Ep : χ(id)=1, χ(R3)=ω,  χ(R3²)=ω²
    Em : χ(id)=1, χ(R3)=ω², χ(R3²)=ω

Character of symmetry (rotation index k, translation by sites[l]=(n1,n2)):
    χ = exp(2πi·(q1·n1 + q2·n2)) × χ_ρ(R3^k)

Symmetry group ordering (block structure, duplicates removed):
    block k=0 : pure translations  (N elements)
    block k=1 : R3  ∘ translations
    block k=2 : R3² ∘ translations

Allowed symmetries:
    C3 little group → all |G| indices
    C1 little group → only k=0 indices (pure translations)

Usage
-----
    python genlat.py [a [b]] [--toml FILE]

    Defaults: a=3, b=1  (7-site C3 cluster)

Examples
--------
    python genlat.py 3 0 --toml triangular.9.toml
    python genlat.py 4 2 --toml triangular.12.toml
    python genlat.py 3 1 --toml triangular.7.toml
"""

import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_a1    = np.array([1.0, 0.0])
_a2    = np.array([0.5, np.sqrt(3.0) / 2.0])
_SQRT3 = np.sqrt(3.0)

# 120° CCW rotation in *lattice* coordinates
_R3_lat = np.array([[-1, -1], [1, 0]], dtype=int)

# Cube root of unity
_OMEGA = np.exp(2j * np.pi / 3.0)

# C3 character table:  (label, [χ(id), χ(R3), χ(R3²)])
_C3_IRREPS = [
    ('A',  np.array([1.0+0j,     1.0+0j,     1.0+0j    ])),
    ('Ep', np.array([1.0+0j,     _OMEGA,     _OMEGA**2  ])),
    ('Em', np.array([1.0+0j,     _OMEGA**2,  _OMEGA     ])),
]


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def lat_to_cart(v):
    """Integer lattice coordinates (n1, n2) → Cartesian 2-vector."""
    return int(v[0]) * _a1 + int(v[1]) * _a2


def _recip_cart(q1, q2):
    """
    Reduced reciprocal coordinates (q1, q2) → Cartesian (kx, ky).

    k = q1·b1 + q2·b2  with  b1=2π(1,-1/√3),  b2=2π(0,2/√3).
    """
    kx = 2.0 * np.pi * q1
    ky = 2.0 * np.pi * (-q1 / _SQRT3 + 2.0 * q2 / _SQRT3)
    return kx, ky


# ---------------------------------------------------------------------------
# Torus construction
# ---------------------------------------------------------------------------

def build_torus(a, b):
    """
    Build a C3-compatible triangular lattice torus with parameters (a, b).

    Returns
    -------
    sites  : list of N tuples (n1, n2), sorted by N·T⁻¹·v lexicographically
    T      : (2,2) int ndarray — torus matrix (columns = torus vectors)
    TinvN  : (2,2) int ndarray — N·T⁻¹ = adjugate(T) = [[a-b,-b],[b,a]]
    N      : int — number of sites
    """
    N = a * a - a * b + b * b
    if N <= 0:
        raise ValueError(f"(a={a}, b={b}) gives N = {N} ≤ 0.")
    T     = np.array([[a, b], [-b, a - b]], dtype=int)
    TinvN = np.array([[a - b, -b], [b, a]], dtype=int)

    rng = abs(a) + abs(b) + 1
    sites = []
    for n2 in range(-rng, rng + 1):
        for n1 in range(-rng, rng + 1):
            f = TinvN @ np.array([n1, n2], dtype=int)
            if 0 <= int(f[0]) < N and 0 <= int(f[1]) < N:
                sites.append((n1, n2))

    if len(sites) != N:
        raise RuntimeError(f"Expected {N} sites, found {len(sites)}.")

    sites.sort(key=lambda v: tuple(int(x) for x in TinvN @ np.array(v, dtype=int)))
    return sites, T, TinvN, N


# ---------------------------------------------------------------------------
# Permutation helpers
# ---------------------------------------------------------------------------

def _compose(p, q):
    """Compose permutations: result[i] = q[p[i]]."""
    return [q[p[i]] for i in range(len(p))]


def _frac_mod(v, TinvN, N):
    """Fractional-coordinate tuple of v, component-wise mod N."""
    f = TinvN @ np.array(v, dtype=int)
    return (int(f[0]) % N, int(f[1]) % N)


def _build_lookup(sites, TinvN, N):
    """Dict: fractional-coord tuple → site index."""
    return {_frac_mod(s, TinvN, N): idx for idx, s in enumerate(sites)}


def compute_translation_perm(sites, TinvN, N, delta):
    """Permutation for translation by lattice vector delta = (dn1, dn2)."""
    lookup = _build_lookup(sites, TinvN, N)
    return [lookup[_frac_mod((s[0]+delta[0], s[1]+delta[1]), TinvN, N)]
            for s in sites]


def compute_rotation_perm(sites, TinvN, N):
    """Permutation for the C3 (120° CCW) rotation."""
    lookup = _build_lookup(sites, TinvN, N)
    return [lookup[_frac_mod(tuple(_R3_lat @ np.array(s, dtype=int)), TinvN, N)]
            for s in sites]


def compute_nn_bonds(sites, TinvN, N):
    """
    Compute all nearest-neighbour bonds on the triangular lattice torus.

    The six NN directions in lattice coordinates are ±(1,0), ±(0,1), ±(1,-1).
    Choosing the three "positive" half-directions avoids double-counting:
        d1 = (1,  0)   along a1
        d2 = (0,  1)   along a2
        d3 = (1, -1)   along a1 − a2  [Cartesian length = 1]

    Returns a list of 3·N ordered pairs (i, j) with i < j is NOT guaranteed
    (the ordering follows site enumeration), but each bond appears exactly once.
    """
    lookup    = _build_lookup(sites, TinvN, N)
    nn_dirs   = [(1, 0), (0, 1), (1, -1)]
    bonds     = []
    for d in nn_dirs:
        for i, s in enumerate(sites):
            j = lookup[_frac_mod((s[0] + d[0], s[1] + d[1]), TinvN, N)]
            bonds.append((i, j))
    return bonds


def generate_all_translations(sites, TinvN, N):
    """All N translation permutations (indexed by displacement = sites[l])."""
    lookup = _build_lookup(sites, TinvN, N)
    return [
        [lookup[_frac_mod((s[0]+d[0], s[1]+d[1]), TinvN, N)] for s in sites]
        for d in sites
    ]


# ---------------------------------------------------------------------------
# Symmetry group with metadata
# ---------------------------------------------------------------------------

def build_sym_group(sites, TinvN, N, perm_rot):
    """
    Build the full T_N × C3 symmetry group with per-element metadata.

    Ordering:  block k=0 (pure translations), k=1 (R3∘trans), k=2 (R3²∘trans).
    Duplicates (e.g. when R3 acts trivially on sites) are dropped;
    first occurrence wins.

    Returns
    -------
    all_syms : list of permutations (each a list of ints)
    sym_meta : list of (k, l)
        k ∈ {0,1,2}   — rotation index
        l ∈ {0,…,N-1} — translation index (displacement = sites[l])
    """
    all_trans  = generate_all_translations(sites, TinvN, N)
    r2         = _compose(perm_rot, perm_rot)
    rot_powers = [list(range(N)), perm_rot, r2]

    all_syms = []
    sym_meta = []
    seen     = set()

    for k, rot in enumerate(rot_powers):
        for l, trans in enumerate(all_trans):
            key = tuple(_compose(trans, rot))
            if key not in seen:
                seen.add(key)
                all_syms.append(list(key))
                sym_meta.append((k, l))

    return all_syms, sym_meta


# ---------------------------------------------------------------------------
# Momentum enumeration
# ---------------------------------------------------------------------------

def _r3_star_int(fq1, fq2, N):
    """
    R3* on integer reduced reciprocal coords.

    R3*: (q1,q2) → (-q2, q1-q2) mod 1
    Integer form: (fq1,fq2) → ((-fq2)%N, (fq1-fq2)%N)
    """
    return ((-fq2) % N, (fq1 - fq2) % N)


def _is_c3_fixed(fq1, fq2, N):
    """True if momentum (fq1/N, fq2/N) is a fixed point of R3*."""
    r1, r2 = _r3_star_int(fq1, fq2, N)
    return r1 == fq1 and r2 == fq2


def enumerate_momenta(a, b, N):
    """
    All N allowed momenta on the torus.

    Dual torus matrix TinvN^T = [[a-b,b],[-b,a]] maps (m1,m2) → fq = N·q.
    Condition: 0 ≤ fq1, fq2 < N.

    Returns list of (fq1, fq2, q1, q2, kx, ky) sorted by (fq1, fq2).
    """
    rng = abs(a) + abs(b) + 1
    result = []
    for m2 in range(-rng, rng + 1):
        for m1 in range(-rng, rng + 1):
            fq1 = (a - b) * m1 + b * m2
            fq2 = -b       * m1 + a * m2
            if 0 <= fq1 < N and 0 <= fq2 < N:
                q1, q2 = fq1 / N, fq2 / N
                kx, ky = _recip_cart(q1, q2)
                result.append((fq1, fq2, q1, q2, kx, ky))

    assert len(result) == N, f"Expected {N} momenta, got {len(result)}"
    result.sort(key=lambda x: (x[0], x[1]))
    return result


# ---------------------------------------------------------------------------
# Irreducible representation sections
# ---------------------------------------------------------------------------

def compute_irrep_sections(a, b, N, sites, all_syms, sym_meta):
    """
    Compute irrep sections for every momentum on the torus.

    C3-fixed momenta (Γ, K, K'):
        Three sections per point — A, Ep, Em — with all |G| symmetries allowed.
        Character formula:  χ_j = exp(2πi q·sites[l]) × χ_ρ(R3^k)
        where sym_meta[j] = (k, l).

    Generic C3-orbits {q, R3*q, R3²*q}:
        One C1.A section per orbit member.  Only the N pure-translation
        symmetries (k=0 block) are allowed.
        Character formula:  χ_j = exp(2πi q·sites[l])  (rotation char = 1).

    Returns list of dicts:
        name, characters (list of complex), allowed_symmetries, momentum [kx,ky]
    """
    momenta       = enumerate_momenta(a, b, N)
    sections      = []
    processed     = set()
    generic_count = 0

    all_allowed   = list(range(len(all_syms)))
    trans_allowed = [j for j, (k, _) in enumerate(sym_meta) if k == 0]

    for fq1, fq2, q1, q2, kx, ky in momenta:
        fq = (fq1, fq2)
        if fq in processed:
            continue

        if _is_c3_fixed(fq1, fq2, N):
            # ---- C3 fixed point: Γ, K, or K' ----
            processed.add(fq)

            if fq1 == 0 and fq2 == 0:
                pt = "Gamma"
            elif 3 * fq1 == 2 * N and 3 * fq2 == N:
                pt = "K"
            elif 3 * fq1 == N and 3 * fq2 == 2 * N:
                pt = "Kp"
            else:
                pt = f"FixedPt_{fq1}_{fq2}"

            for lbl, rot_chars in _C3_IRREPS:
                chars = [
                    np.exp(2j * np.pi * (q1 * sites[l][0] + q2 * sites[l][1]))
                    * rot_chars[k]
                    for k, l in sym_meta
                ]
                sections.append({
                    'name': f"{pt}.C3.{lbl}",
                    'characters': chars,
                    'allowed_symmetries': all_allowed,
                    'momentum': [kx, ky],
                })

        else:
            # ---- generic orbit of 3 momenta: little group = C1 ----
            orbit = []
            cur = fq
            for _ in range(3):
                orbit.append(cur)
                cur = _r3_star_int(cur[0], cur[1], N)
            assert cur == fq, "Orbit did not close after 3 steps"

            for o in orbit:
                processed.add(o)

            for midx, (ofq1, ofq2) in enumerate(orbit):
                oq1, oq2 = ofq1 / N, ofq2 / N
                okx, oky = _recip_cart(oq1, oq2)
                chars = [
                    np.exp(2j * np.pi * (oq1 * sites[sym_meta[j][1]][0]
                                       + oq2 * sites[sym_meta[j][1]][1]))
                    for j in trans_allowed
                ]
                sections.append({
                    'name': f"M{generic_count}_{midx}.C1.A",
                    'characters': chars,
                    'allowed_symmetries': trans_allowed,
                    'momentum': [okx, oky],
                })

            generic_count += 1

    return sections


# ---------------------------------------------------------------------------
# Verification helpers
# ---------------------------------------------------------------------------

def verify_perm(perm, name="perm"):
    N = len(perm)
    assert sorted(perm) == list(range(N)), \
        f"'{name}' is not a valid permutation of {{0,…,{N-1}}}!"


def verify_c3(perm_rot, N):
    """Assert R3^3 = identity."""
    r3 = _compose(_compose(perm_rot, perm_rot), perm_rot)
    if r3 != list(range(N)):
        raise AssertionError("perm_rot does not satisfy R3^3 = id.")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def print_lattice_info(sites, T, N, a, b):
    t1c = lat_to_cart(T[:, 0])
    t2c = lat_to_cart(T[:, 1])
    print("# Triangular (hexagonal) lattice — C3-compatible torus")
    print(f"# Parameters : a = {a},  b = {b}")
    print(f"# Sites      : N = {N}")
    print(f"# Torus vectors:")
    print(f"#   t1 = {T[0,0]}*a1 + ({T[1,0]})*a2  →  Cartesian {t1c}")
    print(f"#   t2 = {T[0,1]}*a1 + ({T[1,1]})*a2  →  Cartesian {t2c}")
    print("# 120° rotation R3 in lattice coords: [[-1,-1],[1,0]]")
    print()
    print(f"# {'idx':>4}   {'(n1,n2)':>10}   {'Cartesian (x, y)':>28}")
    print(f"# {'-'*4}   {'-'*10}   {'-'*28}")
    for idx, s in enumerate(sites):
        c = lat_to_cart(s)
        print(f"  {idx:4d}   ({s[0]:3d},{s[1]:3d})      "
              f"({c[0]: .10f},  {c[1]: .10f})")


def write_toml(filename, sites, T, TinvN, N, a, b, perm_ta1, perm_ta2, perm_rot,
               coupling_name="J1"):
    """
    Write a complete xdiag-compatible TOML lattice file with irrep sections.

    Sections written:
      Coordinates, Interactions (NN Heisenberg bonds),
      Translation_a1/a2, Rotation_C3, Symmetries,
      then one [Name.LittleGroup.Irrep] section per (momentum, irrep) pair.

    Parameters
    ----------
    coupling_name : str
        Label used for every nearest-neighbour bond (default 'J1').
    """
    t1, t2 = T[:, 0], T[:, 1]

    all_syms, sym_meta = build_sym_group(sites, TinvN, N, perm_rot)
    irrep_sections     = compute_irrep_sections(a, b, N, sites, all_syms, sym_meta)
    nn_bonds           = compute_nn_bonds(sites, TinvN, N)

    with open(filename, "w") as f:
        f.write("# Triangular (hexagonal) lattice with C3 (threefold) rotation symmetry\n")
        f.write(f"# Parameters: a={a}, b={b},  N={N} sites\n")
        f.write(f"# Torus vectors: t1=({t1[0]},{t1[1]}) lat,  "
                f"t2=({t2[0]},{t2[1]}) lat\n")
        f.write("# Generated by genlat.py\n\n")

        # Coordinates
        f.write("Coordinates = [\n")
        for s in sites:
            c = lat_to_cart(s)
            f.write(f"  [{c[0]:.15f}, {c[1]:.15f}],\n")
        f.write("]\n\n")

        # Nearest-neighbour Heisenberg interactions
        f.write(f"# Nearest-neighbour Heisenberg interactions ({len(nn_bonds)} bonds)\n")
        f.write("Interactions = [\n")
        for i, j in nn_bonds:
            f.write(f"  ['{coupling_name}', 'SdotS', {i}, {j}],\n")
        f.write("]\n\n")

        # Translation generators
        f.write("# Translation by primitive vector a1 = (1,0) in lattice coords\n")
        f.write(f"Translation_a1 = {perm_ta1}\n\n")
        f.write("# Translation by primitive vector a2 = (0,1) in lattice coords\n")
        f.write(f"Translation_a2 = {perm_ta2}\n\n")

        # Rotation
        f.write("# 120-degree counterclockwise rotation (C3 generator)\n")
        f.write(f"Rotation_C3 = {perm_rot}\n\n")

        # Full symmetry group
        f.write(f"# Full translation × C3 symmetry group ({len(all_syms)} elements)\n")
        f.write("Symmetries = [\n")
        for sym in all_syms:
            f.write(f"  {sym},\n")
        f.write("]\n\n")

        # Irreducible representations
        f.write("# Irreducible representations\n")
        for sec in irrep_sections:
            f.write(f"\n[{sec['name']}]\n")
            f.write("characters = [\n")
            for chi in sec['characters']:
                f.write(f"  [{chi.real:.16f}, {chi.imag:.16f}],\n")
            f.write("]\n")
            f.write(f"allowed_symmetries = {sec['allowed_symmetries']}\n")
            kx, ky = sec['momentum']
            f.write(f"momentum = [{kx:.16f}, {ky:.16f}]\n")

    print(f"\n# Wrote: {filename}  "
          f"({len(all_syms)} symmetry elements, "
          f"{len(irrep_sections)} irrep sections)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("a", nargs="?", type=int, default=3,
                        help="Torus parameter a  (default: 3)")
    parser.add_argument("b", nargs="?", type=int, default=1,
                        help="Torus parameter b  (default: 1)")
    parser.add_argument("--toml", metavar="FILE", default=None,
                        help="Write an xdiag-compatible TOML lattice file")
    parser.add_argument("--coupling-name", metavar="NAME", default="J1",
                        help="Coupling label for NN bonds in TOML (default: J1)")
    args = parser.parse_args()

    a, b = args.a, args.b

    # 1. Build the hexagonal lattice torus
    sites, T, TinvN, N = build_torus(a, b)

    # 2. Symmetry permutations
    perm_ta1 = compute_translation_perm(sites, TinvN, N, (1, 0))
    perm_ta2 = compute_translation_perm(sites, TinvN, N, (0, 1))
    perm_rot = compute_rotation_perm(sites, TinvN, N)

    # 3. Verify
    verify_perm(perm_ta1, "Translation a1")
    verify_perm(perm_ta2, "Translation a2")
    verify_perm(perm_rot,  "Rotation C3")
    verify_c3(perm_rot, N)

    # 4. Print lattice and permutations
    print_lattice_info(sites, T, N, a, b)
    print(f"\n# Translation by a1 = (1,0) [lattice]:")
    print(f"  perm_ta1 = {perm_ta1}")
    print(f"\n# Translation by a2 = (0,1) [lattice]:")
    print(f"  perm_ta2 = {perm_ta2}")
    print(f"\n# 120-degree CCW rotation (C3):")
    print(f"  perm_rot = {perm_rot}")
    print(f"\n# Verification: R3^3 = identity ✓")

    # 5. Optionally write TOML with irrep sections
    if args.toml:
        write_toml(args.toml, sites, T, TinvN, N, a, b,
                   perm_ta1, perm_ta2, perm_rot,
                   coupling_name=args.coupling_name)

    return sites, T, N, perm_ta1, perm_ta2, perm_rot


if __name__ == "__main__":
    main()


