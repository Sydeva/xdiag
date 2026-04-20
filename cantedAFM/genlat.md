# genlat.py — Irrep Extension Design Notes

## Goal

Extend `genlat.py` so that when `--toml FILE` is given it also writes, for
every distinct momentum point on the torus, the irreducible representation
(irrep) sections in the xdiag TOML format:

```toml
[Gamma.C3.A]
characters         = [[re, im], ...]
allowed_symmetries = [0, 1, 2, ...]
momentum           = [kx, ky]
```

---

## 1  Reciprocal lattice

Primitive lattice vectors (Cartesian):

```
a1 = (1, 0)
a2 = (1/2, √3/2)
```

Reciprocal primitive vectors satisfying `ai · bj = 2π δij`:

```
b1 = 2π (1, -1/√3)
b2 = 2π (0,  2/√3)
```

A momentum in **reduced coordinates** `(q1, q2)` corresponds to Cartesian
wavevector `k = q1 b1 + q2 b2`.

---

## 2  Allowed momenta on the torus

The torus matrix `T = [[a, b], [-b, a-b]]` (columns = torus vectors in
lattice coordinates) imposes periodic boundary conditions.  The N allowed
momenta are indexed by integer pairs `(m1, m2)` satisfying

```
0 ≤ [(a-b) m1 - b m2] < N
0 ≤ [  b   m1 + a m2] < N
```

i.e. `TinvN_T @ (m1, m2)` lies in `[0, N)²`, where

```
TinvN_T = transpose of (N · T⁻¹) = [[a-b,  b],
                                     [ -b,  a]]
```

The reduced coordinates of momentum `m = (m1, m2)` are:

```
q = (1/N) · TinvN_T @ m   ∈ [0,1)²
```

The Cartesian wavevector is `k = q1 b1 + q2 b2`.

---

## 3  Action of R3 on momenta

In the reciprocal-lattice basis `(q1, q2)` the 120° rotation acts as:

```
R3* : (q1, q2)  ↦  (-q2,  q1 - q2)   (mod 1)
```

This has order 3: applying it three times returns the identity.

---

## 4  Little group at each momentum

The full point group of the C3 torus is `C3 = {id, R3, R3²}`.  The
**little group** of momentum `q` is the subgroup that fixes `q` modulo the
reciprocal lattice:

| Condition | Little group | Fixed points |
|-----------|-------------|--------------|
| `R3*(q) ≡ q  (mod 1)` | **C3** | Γ = (0,0), K = (2/3,1/3), K' = (1/3,2/3) |
| otherwise | **C1 = {id}** | all generic momenta |

The K and K' points exist on the torus **only when `N ≡ 0 (mod 3)`**, i.e.
when `(a mod 3, b mod 3) ∈ {(1,2),(2,1),(0,0)}`.

---

## 5  Irreps of the little groups

### C1 (trivial little group)
One irrep labelled **A**, character = 1 for every element.

### C3
Three 1-D irreps (using `ω = e^{2πi/3}`):

| Label | χ(id) | χ(R3) | χ(R3²) |
|-------|-------|--------|--------|
| A     | 1     | 1      | 1      |
| E+    | 1     | ω      | ω²     |
| E-    | 1     | ω²     | ω      |

---

## 6  Symmetry group ordering and character formula

The full symmetry group has `|G|` elements built as:

```
sym[k*N + l]  =  compose(all_trans[l],  rot_power[k])
```

where `k ∈ {0,1,2}` indexes `{id, R3, R3²}` and `l ∈ {0,…,N-1}` indexes
the N translations (translation by `sites[l]`).  Duplicate permutations
(possible when N=3) are dropped via a `seen` set; the metadata `(k, l)` is
stored alongside each unique symmetry.

**Character of symmetry `(k, l)` in irrep `(q, ρ)`:**

```
χ = exp(2πi · q · sites[l])  ×  χ_ρ(R3^k)
```

where `q · sites[l] = q1 n1 + q2 n2` (dot product in reduced coordinates)
and `χ_ρ(R3^k)` is the C3 character from the table above (or 1 for C1).

---

## 7  Allowed symmetries

`allowed_symmetries` lists the indices (into `Symmetries`) of every symmetry
whose rotation part belongs to the little group of `q`:

- **C3 little group**: all `|G|` symmetry indices are allowed.
- **C1 little group**: only the `N` indices with `k = 0` (pure translations)
  are allowed.

---

## 8  Momentum orbit structure and naming

Momenta are grouped into C3-orbits `{q, R3*q, R3²*q}`.  Orbit sizes:

| Orbit type | Size | Members |
|-----------|------|---------|
| Γ         | 1    | (0,0) |
| K         | 1    | (2/3,1/3) (if N%3==0) |
| K'        | 1    | (1/3,2/3) (if N%3==0) |
| Generic   | 3    | {q, R3*q, R3²*q} |

**TOML key naming convention:**

| Momentum | Little group | Irrep | TOML key |
|----------|-------------|-------|----------|
| Γ        | C3          | A     | `Gamma.C3.A` |
| Γ        | C3          | E+    | `Gamma.C3.Ep` |
| Γ        | C3          | E-    | `Gamma.C3.Em` |
| K        | C3          | A/E+/E- | `K.C3.A`, `K.C3.Ep`, `K.C3.Em` |
| K'       | C3          | A/E+/E- | `Kp.C3.A`, `Kp.C3.Ep`, `Kp.C3.Em` |
| generic  | C1          | A     | `M{i}_0.C1.A`, `M{i}_1.C1.A`, `M{i}_2.C1.A` |

---

## 9  Implementation plan for `genlat.py`

1. **`_recip_cart(q1, q2)`** — convert reduced reciprocal coords to Cartesian.
2. **`enumerate_momenta(TinvN, N)`** — return list of `(m1,m2, q1,q2, kx,ky)`.
3. **`r3_star(q1, q2)`** — apply R3* in reduced coords, return `((-q2)%1, (q1-q2)%1)`.
4. **`little_group_order(q1, q2)`** — return 3 if fixed by R3*, else 1.
5. **`build_sym_group(sites, TinvN, N, perm_rot)`** — return `(all_syms, sym_meta)`
   where `sym_meta[j] = (k, l)`.
6. **`compute_irrep_sections(momenta, sites, sym_meta, N)`** — return list of dicts
   `{name, characters, allowed_symmetries, momentum}`.
7. **`write_toml(...)`** — append irrep sections after the `Symmetries` block.

---

## 10  Edge cases

- **N=3** (e.g. a=2, b=1): R3 acts as the identity permutation on the 3
  sites; `all_syms` has only 3 elements (not 9).  The K-point exists (3%3=0)
  but the E± irreps become degenerate with A since R3≡id as a permutation.
- **N not divisible by 3**: K and K' are absent; only Γ has a C3 little group.
- **Generic orbit containing only 1 or 2 elements**: impossible for a valid
  C3-compatible torus unless the orbit wraps onto a fixed point, which is
  already handled by the fixed-point detection.


