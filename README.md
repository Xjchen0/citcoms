# CitcomS (Modern C++17 Edition)

[![Version: 0.4.0](https://img.shields.io/badge/Version-citcoms__cpp--v0.4.0-brightgreen.svg)](https://github.com/Xjchen0/citcoms/releases/tag/citcoms_cpp-v0.4.0)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](COPYING)
[![C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/Build-CMake%20%E2%89%A5%203.16-brightgreen.svg)](https://cmake.org)
[![Compiler: GCC & LLVM/Clang](https://img.shields.io/badge/Compiler-GCC%20%7C%20LLVM%2FClang-orange.svg)](https://gcc.gnu.org)
[![CitcomS Manual](https://img.shields.io/badge/docs-CitcomS%20Manual-green.svg)](https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf)

A modernized, high-performance **ISO C++17** finite-element mantle convection engine for **regional (1-cap)** and **global full spherical shell (12-cap)** simulations. Built upon CIG **CitcomS** and integrating core scientific capabilities from **EarthByte/CitcomS**, this repository preserves 100% backward compatibility with legacy configurations while delivering modern Stokes/multigrid solvers, hardware-optimized RDMA communications, geological-age scheduling, and ParaView-ready VTK/HDF5 visualization.

---

## Core Highlights

- **ISO C++17 & CMake Toolchain**: Type-safe C++17 codebase verified on modern GCC (`g++` 9–16) and LLVM (`clang++` 12–22), with native C++17 post-processing (`citcoms_combine`, 50×–100× faster than Python).
- **Modern Stokes & Multigrid Solvers**: Opt-in **W-Norm FGMRES** (`uzawa=fgmres`), **Flexible CG** (`uzawa=fcg`), pure **V-cycle** (`mg_mode=vcycle`), coupled **$3\times 3$ nodal block / Chebyshev smoothers**, coarse-grid PCG, and geometric-mean viscosity coarsening — delivering **2.4×–2.7×** speedup and eliminating Uzawa stalling under $10^4\text{–}10^6$ lateral viscosity contrasts.
- **InfiniBand RDMA / UCX Optimization**: Persistent static communication buffer pools and inverted `MPI_Irecv`-before-`MPI_Isend` ordering eliminate UCX RCache faults and RNR retries at scale.
- **High-Performance I/O, Checkpointing & Exact Age Snapping**: Composite `ascii,vtk` and `ascii,hdf5` modes with automatic 3D ASCII deduplication, ParaView `.pvd` / `.xdmf` timelines, crash-safe sliding-window checkpoints (`restart=auto`), and exact geological-age step snapping (`exact_time_output=on`).
- **EarthByte 4D Plate & Slab Data Assimilation (Opt-in)**: Full support for 4D lithosphere thermal assimilation (age-dependent depth, power-law blending, deforming-zone masking), 3D slab thermal & internal velocity assimilation (synchronized *before* each Stokes solve), EBA depth-dependent $Di(z)$ scaling, 13 lower-mantle viscosity laws (`rheol=20–28, 112–118`), hybrid tracers, dynamic topography buoyancy cutoff, and 9 diagnostic outputs.

---

## Documentation Index

Detailed configuration guides, parameter tables, and tutorials are organized modularly under [`docs/`](docs/):

| Document | Contents |
|---|---|
| **[Build, Toolchain & Post-Processing](docs/build_and_install.md)** | Prerequisites, CMake build flags (GCC/Clang, HDF5, SIMD/LTO, double-precision matrix), InfiniBand RDMA/UCX notes, and `citcoms_combine` usage |
| **[I/O Formats & ParaView Visualization](docs/io_and_visualization.md)** | Single & composite formats (`ascii`, `ascii,vtk`, `ascii,hdf5`), `output_optional*` deduplication, VTK `cap`/`single`/`proc` + `.pvd`, HDF5 + `.xdmf` |
| **[Checkpointing & Geological Age Control](docs/checkpoint_and_time.md)** | `restart=auto`, sliding-window pruning, two-phase crash safety, two-stage binary continuation, `storage_interval_myr` / `end_age_myr` / `exit_at_present` & exact snapping |
| **[Modern Stokes & Multigrid Solvers](docs/stokes_solvers.md)** | Zero-drift kernel optimizations, `FGMRES` / `FCG` / `V-cycle` / `block3x3_gs` / `chebyshev` tuning guide & complete solver parameter table |
| **[EarthByte Scientific Features & Parameters](docs/earthbyte_features.md)** | Complete parameter tables for 4D assimilation, EBA $Di(z)$ scaling, 13 `rheol` viscosity structures (`20–28, 112–118`), hybrid tracers, buoyancy cutoff & diagnostic outputs |
| **[4D 板块与俯冲板片数据同化实战教程 (中文)](docs/plate_assimilation_tutorial.md)** | End-to-end 4D plate/slab data assimilation tutorial: mathematical formulation, Full/Regional file formats, pyGPlates verification script, production `.cfg` templates & FAQ |

> **Parameter Syntax Note**: Following traditional CitcomS conventions, quotes (`"..."`) are only used for path/file/prefix parameters (e.g., `datadir="./output"`, `vel_bound_file="./bvel"`). Non-path options (`output_format=ascii,vtk`, `restart=auto`, `uzawa=fgmres`, `exact_time_output=on`) do **not** require quotes (though the parser accepts both bare and quoted values as well as optional spaces around `=`).

---

## Quick Start

```bash
# Build with CMake (GCC or Clang)
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run Regional (1-cap) simulation
mpirun -np 4 ./CitcomSRegional ../examples/Regional/input.sample

# Run Global Full-Sphere (12-cap) simulation
mpirun -np 12 ./CitcomSFull ../examples/Full/input.sample

# Run automated EarthByte & regression verification tests
cd .. && CITCOMS_BUILD_DIR=./build python3 tests/test_earthbyte_features.py
```

---

## Directory Structure

```text
citcoms/
├── CMakeLists.txt          # CMake build system (C++17)
├── bin/                    # Drivers (CitcomSRegional, CitcomSFull), citcoms_combine, citcoms_h5_to_xdmf.py
├── lib/                    # Core C++17 finite-element, Stokes, MG, advection, tracer & assimilation library
├── Py2C/                   # Standalone C++17 config converter
├── visual/                 # Post-processing wrappers (combine.py, combinesurf.py)
├── docs/                   # Modular documentation & 4D plate assimilation tutorial
├── tests/                  # Automated regression & EarthByte feature test suites
└── examples/               # Regional and Full-sphere example configurations
```

---

## License, Attribution & Academic Citations

### License & Attribution
This project is free and open-source software distributed under the **GNU General Public License v2 (GPL-2.0)** — see [COPYING](COPYING).

- **Upstream Base**: Derived from **CIG CitcomS** (California Institute of Technology / Computational Infrastructure for Geodynamics, `https://github.com/geodynamics/citcoms`).
- **EarthByte Attribution**: The 4D lithosphere/slab thermal and internal velocity assimilation framework, EBA depth-dependent dissipation number scaling, lower-mantle viscosity laws (`rheol=20–28, 112–118`), hybrid tracer composition method, dynamic topography buoyancy cutoff, and extended diagnostic outputs are adapted from **EarthByte/CitcomS** (`https://github.com/EarthByte/citcoms`, licensed under GPL-2.0) and modernized for C++17 and `citcoms_cpp` solvers.

### Citations
If you use this software in your research, please cite the original **CitcomS** papers. If you additionally use the **4D plate/slab assimilation, EBA $Di(z)$ scaling, or EarthByte viscosity/tracer features**, please also cite the corresponding **EarthByte** publications:

1. **CitcomS Core Software & Methods**:
   - Moresi, L., Zhong, S., Han, L., Conrad, C., Tan, E., Gurnis, M., Choi, E., Thoutireddy, P., Manea, V., McNamara, A., Becker, T., Leng, W., & Armendariz, L. (2014). *CitcomS v3.3.1*. Zenodo. [https://doi.org/10.5281/zenodo.7271920](https://doi.org/10.5281/zenodo.7271920)
   - Zhong, S., Zuber, M. T., Moresi, L., & Gurnis, M. (2000). *Role of temperature-dependent viscosity and surface plates in spherical shell models of mantle convection*. Journal of Geophysical Research: Solid Earth, 105(B5), 11063–11082. [https://doi.org/10.1029/2000JB900003](https://doi.org/10.1029/2000JB900003)
   - Tan, E., Choi, E., Thoutireddy, P., Gurnis, M., & Aivazis, M. (2006). *GeoFramework: Coupling multiple models of mantle convection within a computational framework*. Geochemistry, Geophysics, Geosystems, 7(6), Q06001. [https://doi.org/10.1029/2005GC001155](https://doi.org/10.1029/2005GC001155)
   - Zhong, S., McNamara, A., Tan, E., Moresi, L., & Gurnis, M. (2008). *A benchmark study on mantle convection in a 3-D spherical shell using CitcomS*. Geochemistry, Geophysics, Geosystems, 9(10), Q10017. [https://doi.org/10.1029/2008GC002048](https://doi.org/10.1029/2008GC002048)
   - McNamara, A. K., & Zhong, S. (2004). *Thermochemical structures within a spherical mantle: Superplumes or piles?* Journal of Geophysical Research: Solid Earth, 109(B7), B07402. [https://doi.org/10.1029/2003JB002847](https://doi.org/10.1029/2003JB002847)

2. **EarthByte/CitcomS 4D Assimilation, EBA & Viscosity Models**:
   - Bower, D. J., Gurnis, M., & Flament, N. (2015). *Assimilating lithosphere and slab history in 4-D Earth models*. Physics of the Earth and Planetary Interiors, 238, 8–22. [https://doi.org/10.1016/j.pepi.2014.10.013](https://doi.org/10.1016/j.pepi.2014.10.013)
   - Flament, N., Gurnis, M., Williams, S., Seton, M., Skogseid, J., Heine, C., & Müller, R. D. (2014). *Topographic asymmetry of the South Atlantic from global models of mantle flow and lithospheric stretching*. Earth and Planetary Science Letters, 387, 107–119. [https://doi.org/10.1016/j.epsl.2013.11.017](https://doi.org/10.1016/j.epsl.2013.11.017)
   - Zhang, N., Zhong, S., Leng, W., & Li, Z.-X. (2010). *A model for the evolution of the Earth's mantle structure since the Early Paleozoic*. Journal of Geophysical Research: Solid Earth, 115(B6), B06401. [https://doi.org/10.1029/2009JB006896](https://doi.org/10.1029/2009JB006896)
   - Rudolph, M. L., Lekić, V., & Lithgow-Bertelloni, C. (2015). *Viscosity jump in Earth's mid-mantle*. Science, 350(6266), 1349–1352. [https://doi.org/10.1126/science.aad1929](https://doi.org/10.1126/science.aad1929)
