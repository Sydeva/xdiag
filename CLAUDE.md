# CLAUDE.md — XDiag Project Guide

## Project Overview

**XDiag** (v0.4.1) is a C++17 library for exact diagonalization of quantum many-body systems. It supports local spin, t-J, and fermionic (Hubbard) models with full space group symmetry support. Parallelization is provided via OpenMP (shared memory) and MPI (distributed memory). There is also a Julia wrapper (`XDiag.jl`).

- **License**: Apache 2.0 (SPDX headers on all source files)
- **Upstream**: https://github.com/awietek/xdiag
- **Fork**: https://github.com/Sydeva/xdiag
- **Docs**: https://awietek.github.io/xdiag
- **Paper**: https://arxiv.org/abs/2505.02901

---

## Build System

CMake ≥ 3.19. Default build type is **Release**. Default install prefix is `<source_root>/install`.

### Standard build
```bash
cmake -S . -B build
cmake --build build
cmake --install build
```

### Key CMake options

| Option | Default | Description |
|---|---|---|
| `BUILD_TESTING` | Off | Build test suite |
| `BUILD_EXAMPLES` | Off | Build example programs |
| `XDIAG_DISTRIBUTED` | Off | Build MPI-distributed library (`xdiag_distributed`) |
| `XDIAG_JULIA_WRAPPER` | Off | Build Julia wrapper shared library (`xdiagjl`) |
| `XDIAG_DISABLE_OPENMP` | Off | Disable OpenMP |
| `XDIAG_DISABLE_HDF5` | Off | Disable HDF5 I/O |
| `XDIAG_DISABLE_COLOR` | Off | Disable colored terminal output |
| `XDIAG_OPTIMIZE_FOR_NATIVE` | Off | Enable `-march=native -mtune=native` |
| `XDIAG_FORCE_MKL_SEQUENTIAL` | Off | Force Intel MKL into sequential mode |
| `XDIAG_USE_SPARSE_MKL` | Off | Use sparse BLAS from Intel MKL |
| `XDIAG_SHARED_LIBS` | (unset) | Build shared instead of static library |

### Three library targets

| Target | When built | Notes |
|---|---|---|
| `xdiag` | Default | Serial + OpenMP |
| `xdiag_distributed` | `XDIAG_DISTRIBUTED=On` | OpenMP disabled; requires MPI |
| `xdiagjl` | `XDIAG_JULIA_WRAPPER=On` | Always shared; cannot combine with distributed |

### Debug build adds ASan
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```
Debug mode automatically adds `-fsanitize=address`.

### Using the installed library in a downstream project
```cmake
cmake_minimum_required(VERSION 3.19)
project(xdiag_application)

find_package(xdiag REQUIRED HINTS "/path/to/xdiag/install")
add_executable(main main.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag)
```

---

## Dependencies

### Required
- **LAPACK/BLAS** (or Intel MKL as a drop-in replacement; auto-detected)

### Optional (auto-detected)
- **OpenMP** — shared-memory parallelism (disabled for distributed build)
- **HDF5** (C++ interface) — required for `.h5` file I/O
- **MPI** — required only when `XDIAG_DISTRIBUTED=On`
- **Intel MKL** — preferred over plain LAPACK when available; supports GNU, Intel, and sequential threading backends

### Vendored (in `xdiag/extern/`)
- **Armadillo** — dense linear algebra types (`arma::vec`, `arma::mat`, `arma::cx_vec`, …)
- **fmtlib** — string formatting (`{:.12f}` style), header-only (`FMT_HEADER_ONLY`)
- **toml++** — TOML file parsing
- **flat_hash_map** — fast hash map
- **GSL** — guideline support library
- **Clara** — command line parser

---

## Repository Structure

```
xdiag/              # Library source (mirrors installed headers)
  algebra/          # dot, norm, inner, apply, matrix
  algorithms/       # Lanczos, Arnoldi, LOBPCG, time_evolve, sparse_diag
  basis/            # Internal Hilbert space bases (spinhalf, tj, electron, distributed)
  bits/             # Bit manipulation primitives
  blocks/           # Public block types: Spinhalf, tJ, Electron (+ _distributed)
  combinatorics/    # Binomial, combinations, Fermi tables, lin_table
  extern/           # Vendored libraries
  io/               # FileToml, FileH5 — read/write TOML and HDF5
  operators/        # Op, OpSum, Coupling + logic (symmetrize, hc, qns, …)
  parallel/         # MPI utilities (allreduce, alltoall, comm_pattern, …)
  random/           # Random number and hash utilities
  states/           # State, ProductState, RandomState
  symmetries/       # Permutation, PermutationGroup, Representation
  utils/            # Error, Log, Scalar, xdiag_api macros
  all.hpp           # Single-header include for the full public API
  common.hpp        # Shared types (complex, int64_t, overload, …)
tests/              # Catch2 test suite
examples/           # Standalone example programs
benchmarks/         # Performance benchmarks
cmake/              # FindMKL, sources.cmake, config templates, package configs
julia/              # Julia wrapper source
docs/               # MkDocs documentation source
```

---

## Core Concepts & Public API

### Include everything
```cpp
#include <xdiag/all.hpp>
using namespace xdiag;
```

### Blocks (Hilbert spaces)
```cpp
Spinhalf block(nsites);                     // all Sz sectors
Spinhalf block(nsites, nup);                // fixed Sz = nup - ndown
Spinhalf block(nsites, irrep);              // with space group symmetry
Spinhalf block(nsites, nup, irrep);         // Sz + symmetry

tJ       block(nsites, nup, ndn);
Electron block(nsites, nup, ndn);

// Distributed (MPI) variants:
// SpinhalfDistributed, tJDistributed, ElectronDistributed
```

### Operators
```cpp
Op op("SdotS", {i, j});         // two-site operator by type string
Op op("Sz", i);                 // single-site

OpSum H;
H += "J" * Op("SdotS", {i, j});
H["J"] = 1.0;                   // set named coupling constant
```

Built-in operator type strings include: `"SdotS"`, `"S+"`, `"S-"`, `"Sz"`, `"Hop"`,
`"Cdagup"`, `"Cup"`, `"Cdagdn"`, `"Cdn"`, `"Nup"`, `"Ndn"`, `"tJSdotS"`, etc.
Custom matrix operators are also supported via `arma::mat` / `arma::cx_mat`.

### States
```cpp
State psi(block);                                       // zero state
State psi(block, /* real= */ true, ncols);
State psi = product_state(block, {"Up","Dn","Up",...});
State psi = random_state(block);
```
Underlying storage is `arma::vec` / `arma::cx_vec`; retrieve via
`psi.vector()`, `psi.vectorC()`, `psi.matrix()`, `psi.matrixC()`.

### Algorithms
```cpp
double e0          = eigval0(ops, block);           // lowest eigenvalue
auto [e0, psi]     = eig0(ops, block);              // + eigenvector

EigvalsLanczosResult res = eigvals_lanczos(ops, block, neigvals);
EigsLanczosResult    res = eigs_lanczos(ops, block, neigvals);

State psi_t = time_evolve(H, psi, t);               // real-time evolution
State psi_b = imaginary_time_evolve(H, psi, tau);   // imaginary-time
```

All algorithms accept an optional `precision` (default `1e-12`),
`max_iterations` (default `1000`), and `random_seed` (default `42`).

### Algebra
```cpp
double   n  = norm(psi);
double   d  = dot(psi, phi);
complex  dc = dotC(psi, phi);
double   e  = inner(H, psi);        // <psi|H|psi>
complex  ec = innerC(H, psi);
State    Hv = apply(H, psi);
arma::mat M = matrix(H, block);     // dense matrix representation
```

### Symmetries
```cpp
Permutation p({1, 2, 3, 0});                              // site 0->1, 1->2, …
PermutationGroup group({p, p*p, p*p*p, Permutation(4)});  // cyclic group C4
Representation irrep(group, {1.0, -1.0, 1.0, -1.0});      // characters
Spinhalf block(nsites, nup, irrep);
```

### I/O
```cpp
FileToml fin("params.toml", "r");
FileH5   fout("data.h5", "w");
```

---

## Error Handling Convention

Every `main()` uses the try/catch pattern:
```cpp
int main() try {
  // ...
} catch (Error e) {
  error_trace(e);   // prints full stack trace
}
```

Internal code uses the macros:
```cpp
XDIAG_THROW("message");        // throw xdiag::Error with file/line info
XDIAG_RETHROW(e);              // rethrow, adding caller context
XDIAG_TRY_CATCH(statement);   // wrap a single call
```

---

## Logging & Verbosity

```cpp
set_verbosity(2);                         // 0=silent, 1=normal, 2=verbose
Log("Ground state energy: {:.12f}", e0);  // fmt-style, MPI-aware
```

`Log` resolves to `LogSerial` or `LogMPI` depending on build configuration.

---

## Testing

```bash
# Build and run serial tests
cmake -S . -B build -DBUILD_TESTING=On
cmake --build build
cd build && ctest

# Build and run distributed (MPI) tests
cmake -S . -B build -DBUILD_TESTING=On -DXDIAG_DISTRIBUTED=On
cmake --build build
cd build && mpirun -n 4 ./tests_distributed
```

- Uses **Catch2** (bundled as `tests/catch.hpp`)
- Main runners: `tests/catch_main.hpp` (serial), `tests/catch_mpi_main.hpp` (MPI)
- Tests mirror the source tree under `tests/`
- Serial binary: `tests`; distributed binary: `tests_distributed`

---

## Code Conventions

- **Standard**: C++17, no compiler extensions (`CMAKE_CXX_EXTENSIONS OFF`)
- **Namespace**: all symbols live in `namespace xdiag`
- **Header guards**: `#pragma once`
- **Public API**: all exported symbols marked `XDIAG_API`
- **Integer type**: `int64_t` used throughout; MKL uses ILP64 (`MKL_ILP64`)
- **Complex type**: `using complex = std::complex<double>` (defined in `common.hpp`)
- **Index offset**: `XDIAG_OFFSET` is `1` for the Julia wrapper (1-based indexing), `0` otherwise
- **Sentinel values**: `invalid_index = (int64_t)-1`, `undefined = std::numeric_limits<int64_t>::min()`
- **License headers**: every source file begins with an SPDX comment block
- **IntelLLVM compiler**: uses `-fp-model precise` for floating-point reproducibility
- **Armadillo** is the dense linear algebra type system throughout; avoid raw arrays in new code

---

## Common Pitfalls

- `XDIAG_DISTRIBUTED` and `XDIAG_JULIA_WRAPPER` are **mutually exclusive**.
- When `XDIAG_DISTRIBUTED=On`, OpenMP is **automatically disabled** in the distributed library.
- The Julia wrapper **forces** `BUILD_SHARED_LIBS=On`.
- HDF5 is **silently ignored** when building the Julia wrapper even if found.
- `config.hpp` is **generated** at CMake configure time from `cmake/config.hpp.in` — do not edit it manually.
- New source files must be **registered** in `cmake/sources.cmake` (in `XDIAG_SOURCES`, `XDIAG_DISTRIBUTED_SOURCES`, or `XDIAG_JULIA_SOURCES` as appropriate).

