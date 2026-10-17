# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Bertini_real is a numerical algorithm for decomposing real algebraic curves (dim 1) and
surfaces (dim 2) into cell complexes, using **Bertini 1** as the homotopy-continuation engine.
The C++ side produces the decomposition; the Python package (`python/`) parses the on-disk
output and renders/plots it. MATLAB code in `matlab_codes/` is legacy post-processing.

Two C++ executables are built: `bertini_real` (computes the decomposition) and `sampler`
(refines/adaptively samples an existing decomposition for smoother output).

## Building (C++)

As of July 2025 the build is **CMake** (previously autotools), and it requires **Bertini 1
version >= 1.7** (which itself moved to CMake — install location and headers changed).

```sh
mkdir build && cd build
cmake ../
make            # use -j for parallelism
make install    # may need sudo; installs bertini_real + sampler to bin/
```

Dependencies (all must be findable by CMake): MPFR, GMP, `bertini-parallel` (>=1.7, found via
`find_package(bertini1 1.7 CONFIG)`), Boost (>=1.53, components `filesystem` + `timer`), MPI
(openmpi/mpich), plus tools CMake, Flex, Bison. **Bertini 1 must be compiled from source
against the same GMP/MPFR/MPI libraries.** Custom CMake find-modules live in `cmake/`.
`brconfig.h.in` is configured into `build/config.h`. Source-file lists are in `files.cmake`
(not globbed) — **add new `.cpp`/`.hpp` files there**, not just to disk.

Clone recursively (`git clone --recursive`) to get the `matlab_codes/brakelab` submodule.

## Python package

```sh
cd python && pip install -e .          # installs the `bertini_real` package
```
Deps: matplotlib, trimesh, dill, algopy, sympy, scipy, networkx. `glumpy` (OpenGL rendering)
is optional and imported defensively in `__init__.py`.

Typical interactive use, run **from inside a decomposition output folder**:
```python
import bertini_real
bertini_real.gather_and_plot()   # data.gather() -> plot.plot()
```
`data.gather()` reads the raw on-disk output into a `Curve` or `Surface` object (chosen by the
dimension in the directory name); `data.gather_and_save()` also dills it to a `BRdataN.pkl`.

## Running / tests

There is **no automated test suite or test runner**. `test/curve/*` and `test/surface/*` are
example systems, each a directory containing a Bertini `input` file (and sometimes Python
plot/assemble scripts). The manual workflow for any example:

1. Run Bertini 1 on the `input` file with `tracktype: 1` to produce a numerical irreducible
   decomposition and a `witness_data` file.
2. Run `bertini_real` in that directory; it consumes `input` + `witness_data`. If there are
   multiple components it prompts for which to decompose. Optional flags:
   `-sphere <file> -pi <projectionfile>`.
3. Optionally run `sampler` to refine the decomposition.
4. Use the Python package to plot the result.

Output is plain-text files in a subfolder of cwd, written incrementally after each major stage
(so a crash still leaves the last good parsable state). Key files: `decomp`, `vertex_set`,
copies of `input` + `witness_data`; curves add `E.edge`; surfaces add `S.surf` plus curve
sub-decompositions in their own subfolders. The README and `manual/bertini_real_manual.pdf`
document these formats in full.

## C++ architecture

`src/` and `include/` mirror each other and are organized by subsystem (file lists in
`files.cmake`). `src/bertini_real.cpp` and `src/sampler/` hold the two executable `main`s;
everything else compiles into both (`common_src`).

- **`bertini1/`** — `bertini_extensions`: the C++ bridge to Bertini 1's C structures/headers
  (`bertini_headers.hpp`).
- **`nag/`** — numerical algebraic geometry core. `nid` (numerical irreducible decomposition),
  `witness_set`, `system_randomizer`, and `nag/solvers/` (the homotopy solvers: `midpoint`,
  `multilintolin`, `nullspace`, `sphere_intersection`, `postProcessing`, common `solver`).
- **`decompositions/`** — the top-level algorithms: `curve`, `surface`, base `decomposition`,
  and `checkSelfConjugate` (real-vs-complex detection).
- **`cells/`** + **`containers/`** — the cell-complex data model: `vertex`/`edge`/`face`/`cell`,
  held in `vertex_set` and `holders`.
- **`symbolics/`** — symbolic preprocessing: `derivative_systems`, `isosingular` (deflation),
  `nullspace`, `slicing`, `sphere_intersection`.
- **`io/`** — `fileops`, terminal `color`, and the **Flex** parser `partitionParse.l` (compiled
  to `partitionParse.yy.c` at build time; CMake `flex_target` with prefix `partitionParse`).
- top-level: `programConfiguration` (CLI flags / config), `parallelism` (MPI master/worker),
  `double_odometer`, `limbo`.

The program is **MPI-parallel** (head/worker model in `parallelism`). C++14.

## Python architecture

`python/bertini_real/` mirrors the C++ cell model in Python objects. `data/` does the parsing
(`gather*`), `curve`/`surface`/`edge`/`face`/`vertex`/`cell`/`decomposition` are the parsed
types, `parse/` reads the directory naming convention, `plot/`/`glumpyplotter`/`anaglypy` do
rendering, `dehomogenize/` handles projective coords, `paths/` and `util/` are helpers.

## Notes

- `documentation/` holds the Doxygen config (`bertini_real.doxy.config`) for the C++ docs at
  doc.bertinireal.com/cpp; `python/docs/` is the Sphinx source for doc.bertinireal.com/python.
- GitHub Actions only mirrors pushes to an MPI GitLab (`.github/workflows/github-gitlab-sync.yml`);
  there is no CI build/test.
