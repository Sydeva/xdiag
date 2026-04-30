# AGENTS.md — XDiag Structure & Operator Implementation
This document captures the structure of XDiag v0.4.1, with focus on the operator system architecture and performance implications.
---
## Architecture Overview
XDiag is organized around **three library targets**:
- `xdiag` — Serial + OpenMP (default)
- `xdiag_distributed` — MPI-parallelized (mutually exclusive with OpenMP)
- `xdiagjl` — Julia wrapper (shared library only)
Build system uses CMake ≥ 3.19 with **Release** as default build type. Native optimizations (`-march=native -mtune=native`) are **disabled by default** and must be enabled with `-DXDIAG_OPTIMIZE_FOR_NATIVE=On`.
---
## Operator System
### Public API (`xdiag/operators/`)
**Core Classes:**
- **`Op`** (`op.hpp/cpp`) — Represents a single operator term with optional sites and custom matrix
  - Can be built-in type (string) or custom matrix (arma::mat / arma::cx_mat)
  - Supports 1-site, 2-site, 3-site, and N-site operations via sites vector
  - No dedicated type definition in Op itself; type string is unvalidated at construction
- **`Coupling`** (`coupling.hpp/cpp`) — Encodes coupling coefficient
  - Either scalar (real/complex) or string (symbolic, resolved later)
  - Used for parametric Hamiltonians (e.g., `"J1"`, `"D"` resolved via `OpSum["J1"] = 1.0`)
- **`OpSum`** (`opsum.hpp/cpp`) — Container for operator sum
  - Holds `std::vector<std::pair<Coupling, Op>>`
  - Manages symbolic coupling constants in `std::map<string, Scalar>`
  - Supports arithmetic (+=, -=, *= with scalars)
### Operator Type System (`xdiag/operators/logic/types.hpp/cpp`)
**Type Registry** — Compile-time lists of known operator types. Built-in types and Matrix are the main categories. Classification: real types (SzSz, Exchange, SdotS, etc.) vs complex types (ScalarChirality). Site count mapping: undefined (generic), 1-site, 2-site, 3-site.

**Added operators `"S+S+"` and `"S-S-"`** (2-site, real, disjoint sites required):
- `"S+S+"` — acts as S⁺ᵢ S⁺ⱼ: fires only when both sites i and j are spin-down, flips both to spin-up. Changes nup by +2.
- `"S-S-"` — Hermitian conjugate of `"S+S+"`: fires only when both sites are spin-up, flips both to spin-down. Changes nup by −2.
- Implemented in `xdiag/basis/spinhalf/apply/apply_spsp_smsm.hpp` as a hand-optimized bit-manipulation kernel (same performance tier as `Exchange`, `SzSz`).
- `hc("S+S+")` → `"S-S-"` and vice versa (registered in `hc.cpp`).
- Site ordering: site pair is canonicalized to {min, max} by the generic 2-site path in `order.cpp` without sign change (S⁺ᵢS⁺ⱼ = S⁺ⱼS⁺ᵢ for spins on different sites).
- **Usage note**: because nup changes by ±2, these operators can only be applied to `Spinhalf` blocks without a fixed `nup` sector, or used as off-diagonal blocks in a sector-crossing context.
**Helper functions:**
- `is_known_type(type)` — Boolean check
- `is_real_type(type)` / `is_cplx_type(type)` — Classification
- `nsites_of_type(type)` — Returns site count (or `undefined` for flexible types)
### Operator Validation (`xdiag/operators/logic/valid.cpp`)
**`check_valid(Op)` dispatch:**
1. Check type is in `known_types`; throw otherwise
2. **Built-in ops** (Sz, S+, S-, Cdagup, Cup, etc.): Must have 1 site, must NOT have matrix
3. **Two-site ops** (Hop, SzSz, SdotS, Exchange, tJSzSz, tJSdotS, NupNdn, etc.): Must have exactly 2 sites, must NOT have matrix
4. **Two-site disjoint ops** (`"S+S+"`, `"S-S-"`): Must have exactly 2 **distinct** sites (i ≠ j), must NOT have matrix — enforced via `must_have_disjoint_sites`
5. **Three-site ops** (ScalarChirality): Must have exactly 3 sites, must NOT have matrix
6. **Parameterless ops** (HubbardU, Id): Must NOT have sites, must NOT have matrix
7. **Matrix type**: Must have matrix, sites must be disjoint, no site count constraint (NxN for N sites)
**Error handling:** Immediate `XDIAG_THROW` on validation failure; used at read time and during algorithm entry points.
---
## Operator Application Pipeline
### Spinhalf Block Implementation
**High-level flow** (`xdiag/basis/spinhalf/apply/apply_terms.hpp`):
- Loop over all (coupling, op) pairs
- Built-in types (SzSz, Exchange, etc.): dispatch to hand-optimized implementations
- Matrix type: call `apply_matrix` which decomposes into non-branching operators
### Built-in vs. Matrix Operators: Performance Comparison
**Built-in operators** (e.g., `apply_szsz.hpp`):
- Hand-optimized bit manipulation (e.g., `bits::popcnt(spins & mask)`)
- Single dispatch, no decomposition overhead
- Example: SzSz on 2 sites computes coefficient in ~3 operations
**Matrix operators** (`apply_matrix.hpp`):
1. **Decomposition phase**: Call `non_branching_ops<bit_t, coeff_t>(cpl, op)`
   - Scans matrix element-by-element
   - Builds lookup tables for each basis state → (new state, coefficient)
   - Returns vector of non-branching operators (up to 4 terms for generic 4×4)
2. **Application phase**: Iterate over sum of non-branching operators
   - Check if diagonal or off-diagonal
   - Incurs ~8× overhead for decomposed 4×4 matrices
**Performance impact for DM operator** (user-specified `Sx_i Sy_j - Sy_i Sx_j`):
- Two Matrix terms at cost of ~8 applications each
- Honeycomb lattice: ~72 DM bonds × 2 terms × 8 decomposition cost ≈ **significant overhead**
- Real-time complexity depends on Lanczos iterations
### Compilation & Lowering (`xdiag/operators/logic/compilation.cpp`)
**Purpose:** Transform high-level OpSum into block-specific representation before diagonalization.
**Spinhalf compilation** (`compile_spinhalf`):
- `SdotS` → `SzSz` + `Exchange`
- Same-site checks produce simplified terms
- Matrix operators pass through unchanged
**Key observation:** Matrix operators are **NOT lowered**. Users must provide block-compatible matrices directly.
---
## Operator Logic Suite (`xdiag/operators/logic/`)
### Ordering (`order.cpp`)
- Canonicalize operator ordering for efficient lookup
- Built-in 2-site ops: sort by site pair, apply conjugation rules
- Matrix ops: permute matrix by site permutation; combine matrices on same sites
- Result: lexicographic ordering by site, then by matrix content
### Quantum numbers (`qns.cpp`)
- Determine particle conservation (nup/ndn change)
- Built-in ops: pre-mapped
- Matrix ops: analyze to determine nup change or return `std::nullopt` if ambiguous
### Hermitian conjugate (`hc.cpp`)
- Built-in ops: specific mappings (S+ ↔ S-, etc.)
- Matrix ops: take Hermitian conjugate; preserve sites
### Real/complex classification (`real.cpp`)
- Built-in ops: pre-classified
- Matrix ops: check `op.matrix().isreal()`
---
## I/O System (`xdiag/io/`)
### TOML Operators (`xdiag/io/toml/operators.cpp`)
**Parsing Op from TOML array:**
- Format: `["OpType", site1, site2, ..., [[matrix_rows]]]`
- Optional matrix (array of arrays) appended
**OpSum parsing**:
- Accepts array of `[coupling, Op_array]` or table with `Interactions` + `Constants`
- Coupling can be scalar or string (symbolic)
### Read functions (`xdiag/io/read.cpp`)
**Public API:**
- `read_opsum(FileToml, tag)` — Loads OpSum from TOML file
- `read_representation(FileToml, irrep_tag, group_tag)` — Loads symmetry irrep
---
## User Extension Point: Custom Operators via Matrix
**Standard workflow:**
1. Define 4×4 complex matrices for single operator terms
2. Create Op with `Op("Matrix", {i,j}, matrix)`
3. Add to OpSum: `ops += coupling * Op(...)`
4. Pass to diagonalization
**Implementation example** (DM rewrite):
```cpp
arma::cx_mat sx_sy = arma::kron(sx, sy);
ops_out += cpl * Op("Matrix", sites, sx_sy);
ops_out += cpl * Op("Matrix", sites, minus_sy_sx);
```
**Constraints:**
- Matrix must be 2^N × 2^N where N = number of sites
- Must be "non-branching" (max 1 nonzero per row) for Lanczos stability
- No automatic lowering to built-in ops
---
## Performance Implications
### Default Release build
- `-O3` enabled
- `-march=native` **disabled** (portable binaries)
- Linear algebra: MKL (threaded) or plain BLAS/LAPACK
### To enable native tuning
```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DXDIAG_OPTIMIZE_FOR_NATIVE=On
cmake --build build && cmake --install build
```
Expected speedup: **5–20%** for Matrix-heavy codes.
---
## Summary: Operator Architecture
| Aspect | Implementation |
|--------|---|
| **Type system** | Hardcoded list; validation at construction |
| **Built-in ops** | Hand-optimized per block type |
| **Matrix ops** | Generic; decomposed into non-branching terms; 8–10× slower |
| **Compilation** | Light lowering; Matrix ops pass through |
| **Customization** | Matrix ops with 2^N × 2^N Hermitian matrices |
| **Performance** | Built-in ops optimal; Matrix ops incur decomposition; native optimizations optional |
| **`"S+S+"` / `"S-S-"`** | Built-in 2-site Spinhalf operators; S⁺ᵢS⁺ⱼ / S⁻ᵢS⁻ⱼ; i≠j enforced; nup ±2; h.c. pair; kernel in `apply_spsp_smsm.hpp` |

---
## Tests for `"S+S+"` / `"S-S-"`

**Test file:** `tests/blocks/spinhalf/test_spinhalf_spsp_smsm.cpp`
**Registered in:** `tests/CMakeLists.txt` (target `tests`)
**Run with:** `./tests/tests spinhalf_spsp_smsm`

**What is tested (124 assertions, lattices N = 2 … 5, all site pairs i ≠ j):**

| Check | Method |
|-------|--------|
| `apply(Op("S+S+",{i,j}), ψ)` matches `S⁺ᵢ(S⁺ⱼ(ψ))` | Compare result vectors to successive single-site `S+` applications |
| `apply(Op("S-S-",{i,j}), ψ)` matches `S⁻ᵢ(S⁻ⱼ(ψ))` | Compare result vectors to successive single-site `S-` applications |
| Adjoint relation `⟨φ\|S+S+\|ψ⟩ = ⟨ψ\|S-S-\|φ⟩` | Real dot-product equality |
| `nup(Op("S+S+",{0,1})) == 2` | Quantum number check |
| `nup(Op("S-S-",{0,1})) == -2` | Quantum number check |
| `hc(Op("S+S+",{0,1})) == Op("S-S-",{0,1})` | Hermitian conjugate mapping |
| `hc(Op("S-S-",{0,1})) == Op("S+S+",{0,1})` | Hermitian conjugate mapping |

The `nup`/`hc` assertions also appear inline in `tests/operators/logic/test_qns.cpp`.

**Important:** Tests use `Spinhalf(nsites)` (no fixed `nup`) because `"S+S+"` / `"S-S-"` change nup by ±2 and cannot be applied within a fixed-nup sector.
