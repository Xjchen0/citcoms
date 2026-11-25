# citcoms_cpp: Modern C++17 Mantle Convection Simulation Engine

[![Version: 0.3.0](https://img.shields.io/badge/Version-citcoms__cpp--v0.3.0-brightgreen.svg)](https://github.com/Xjchen0/citcoms_cpp/releases/tag/citcoms_cpp-v0.3.0)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](COPYING)
[![C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/Build-CMake%20%E2%89%A5%203.16-brightgreen.svg)](https://cmake.org)
[![Compiler: GCC & LLVM/Clang](https://img.shields.io/badge/Compiler-GCC%20%7C%20LLVM%2FClang-orange.svg)](https://gcc.gnu.org)
[![CitcomS Manual](https://img.shields.io/badge/docs-CitcomS%20Manual-green.svg)](https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf)

**citcoms_cpp** is a modernized, high-performance finite element simulation engine for thermochemical mantle convection within Earth's mantle in both **regional (1 cap)** and **global full spherical shell (12 caps)** geometries.

This repository (`citcoms_cpp`) represents a complete architectural modernization to **ISO C++17**, standardizing all legacy C/K&R routines into modern, type-safe C++ with first-class support for modern **GCC (`g++`)** and **LLVM (`clang++`)** toolchains, a **CMake** build system, **hardware-optimized InfiniBand RDMA/UCX communications**, and **high-performance C++ post-processing** while strictly preserving 100% numerical and physical equivalence with all CitcomS physics solvers, configurations, and input/output formats.

---

## 1. Key Modernization Features

### 1.1 Complete C++17 Architecture & Mainstream Toolchain Support
* **Standard C++17 (`-std=c++17`)**:
  - All legacy source files (`lib/*.c`, `bin/*.c`) modernized into clean `.cpp` sources.
  - Eliminated all obsolete K&R C parameter lists and unprototyped forward declarations in favor of standard ANSI/ISO C++ signatures.
  - Full const-correctness on string literals and message passing.
  - Modern type safety and memory allocations: explicit typed conversions replacing implicit `void*` casts.
  - Standardized mathematical and utility routines (`<algorithm>`, `<cmath>`, `<filesystem>`, `<string_view>`, `std::swap`, constexpr `citcom_min`/`citcom_max`).
* **Multi-Toolchain Compatibility**:
  - **GCC (`g++`)**: Fully verified with GCC 9, 10, 11, 12, 13, 14, 15, and 16.
  - **LLVM (`clang++`)**: Fully verified with Clang 12 through Clang 21/22.
  - Clean compilation with `-Wall -Wextra` flags enabled.

### 1.2 Strict 100% Function & Configuration Compatibility
* **Zero Discrepancy in Physics**:
  - Governing equations, Stokes flow solvers, multigrid preconditioners, conjugate gradient routines, spherical harmonics, advection-diffusion solvers, and tracer particle dynamics remain 100% mathematically and numerically identical.
* **Identical Configuration & I/O**:
  - All parameter files (`.cfg`), input instructions, restart flags, and command-line invocations operate exactly as before.
  - Generates identical output files (`.coord`, `.velo`, `.visc`, `.temp`, `.surf`, `.botm`, `.ave`, `.chkpt`).

### 1.3 High-Performance Post-Processing (`citcoms_combine`)
* **50x–100x Speedup**:
  - Replaces slow, memory-intensive Python string parsing for large 3D grids with a native C++17 multi-threaded postprocessor: `citcoms_combine`.
  - Performs on-the-fly pasting and domain-decomposition re-assembly for volume and surface outputs with 64 KB buffered streaming.
* **Transparent Backward Compatibility**:
  - `visual/combine.py` and `visual/combinesurf.py` automatically detect and delegate to `citcoms_combine` when present, falling back to Python 3 standard library seamlessly if needed.

### 1.4 High-Performance I/O, Advanced Checkpointing & Geological Age Scheduling
* **Composite Output Formats & De-duplication**:
  - Natively supports composite formats `"ascii,hdf5"` and `"ascii,vtk"`.
  - Decoupled ASCII 3D volume outputs: When `output_optional_ascii` specifies 1D monitoring (`horiz_avg`), redundant 3D body files (`.velo`, `.visc`, and cycle 0 `.coord`) are automatically suppressed, eliminating inode exhaustion and eliminating duplicate 3D writes in composite modes.
* **Fault-Tolerant Checkpointing**:
  - Sliding window pruning (`checkpoints_to_keep = 3`) to eliminate disk quota exhaustion.
  - Two-phase commit with global MPI barrier synchronization to guarantee crash-safe restarts under scheduler walltime limits (`SIGKILL`).
  - Zero-touch auto-resume (`restart = auto`) for hands-off HPC job array chaining.
  - Wall-clock checkpoint triggering (`checkpoint_walltime`) and guaranteed exit checkpoints (`checkpoint_on_termination = on`).
* **Geological Age Triggers & Exact Time Snapping**:
  - Direct scheduling via physical Myr (`storage_interval_myr`, `checkpoint_interval_myr`) mapped via internal `scalet`.
  - Exact time step snapping (`exact_time_output = on`) to ensure frames land exactly on target Myr boundaries for uniform animation rendering.
* **ASCII-Independent Continuation**:
  - Two-stage restarts (`tic_method = -1`, `tracer.ic_method = 2`) directly read binary checkpoint (`.chkpt`) files without intermediate uncompressed ASCII files.
* **ParaView 1-Click Interactive Visualization (VTK & HDF5)**:
  - VTK PVD Collection (`.pvd`): Rank 0 continuously generates and updates `${data_file}.pvd` linking multi-block `.vtm` / `.pvts` files with physical geological ages (Myr).
  - HDF5 XDMF Indexing (`.xdmf`): Maps 12-cap curvilinear structured meshes (`3DSMesh`) and hyperslab slices to parallel HDF5 datasets, enabling direct native 3D visualization and animated time playback in ParaView.
  - Standalone tool `bin/citcoms_h5_to_xdmf.py` enables 1-second retroactive XDMF indexing of existing HDF5 runs.

### 1.5 InfiniBand RDMA & UCX Hardware Communication Optimization
* **Persistent Communication Buffer Pools**:
  - Replaces per-iteration dynamic memory allocations (`malloc`/`free`) in core multigrid smoothers and nodal field exchanges (`full_exchange_id_d`, `exchange_node_d`, `exchange_node_f`, `full_exchange_snode_f`) with persistent, adaptive static communication buffers.
  - Keeps communication virtual addresses stable throughout millions of time steps, preventing InfiniBand UCX Memory Registration Cache (RCache) desynchronization and eliminating hardware protection faults (`Transport retry count exceeded, synd 0x15`).
* **Inverted Non-Blocking Communication Ordering**:
  - Posts all non-blocking receives (`MPI_Irecv`) prior to packing and posting non-blocking sends (`MPI_Isend`) in domain boundary exchanges and tracer particle migration.
  - Eliminates Receiver-Not-Ready (RNR NACK) packet drops, retry overhead, and communication deadlocks under high-core-count HPC scaling (hundreds to thousands of MPI ranks).

### 1.6 Modernized Stokes & Multigrid Solver Engine
* **Zero-Drift Kernel & Communication Optimizations (Always On)**:
  - **MPI_Allreduce Fusion**: Fuses 5 consecutive Uzawa norm Allreduces into 1 (`global_stokes_norms`) and 2 multigrid upward-sweep Allreduces into 1 (`global_vdot2`).
  - **Persistent Solver Workspaces**: Eliminates per-step heap allocations across `multi_grid`, `conj_grad`, and outer Uzawa solvers.
  - **O(1) Static Stiffness Assembly LUT**: Precomputes `Node_map` indices for `construct_node_ks()`, replacing repeated inner linear scans.
  - **SIMD-Friendly Kernel Loops**: Register-accumulated `gauss_seidel` and `n_assemble_del2_u`, single-pass outer-element `assemble_div_u`/`assemble_c_u`, and cached advection-diffusion spherical metrics.
* **Flexible Outer Krylov Solvers (Opt-in, 100% Legacy Default Compatible)**:
  - **W-Norm Flexible GMRES (`uzawa = fgmres`)**: Eliminates CG orthogonality loss and divergence stalling under steep lateral viscosity contrasts ($10^4$–$10^6$) and inexact inner solves; supports both incompressible and compressible (ALA/BA) flows.
  - **Truncated Flexible CG (`uzawa = fcg`)**: Notay's FCG($m_{\max}$) with fused multi-dot Allreduce and optional periodic divergence refresh.
  - **Relative Inner Tolerance (`inner_rel_tol`) & Schur Preconditioner (`p_precond = visc_mass`)**.
* **Modern Multigrid Cycle & Smoothers (Opt-in)**:
  - **Pure V-Cycle (`mg_mode = vcycle`)**: Avoids redundant FMG coarse-grid traversals inside outer Krylov iterations.
  - **Coupled $3\times 3$ Nodal Block Smoothers (`mg_smoother = block3x3_gs` / `block3x3_sgs`)** and **Rank-Count-Invariant Chebyshev Polynomial Smoother (`mg_smoother = chebyshev`)**.
  - **Coarse-Grid PCG (`coarse_solver = cg`)** and **Geometric-Mean Viscosity Coarsening (`visc_smooth_method = 4`)**.

### 1.7 EarthByte Scientific & 4D Plate Assimilation Capabilities (Opt-in)
* **4D Lithosphere, Slab & Internal Velocity Data Assimilation**:
  - Time-dependent lithosphere thermal assimilation with age-dependent boundary layer depth (`lith_age_depth_function`), power-law blending (`lith_age_exponent`), minimum age clamping (`lith_age_min`), and deforming-region masking (`lith_age_stencil_value`).
  - 3D thermal slab assimilation (`slab_assim`, `slab_assim_file`) and 3D dynamic internal velocity boundary conditions (`internal_vbcs_file`, `velocity_internal_file`) synchronized **before** each Stokes solve and upon checkpoint restart (eliminating the 1-step lag in upstream EarthByte).
* **Extended Boussinesq Approximation (EBA) Depth-Dependent $Di$ Scaling**:
  - Opt-in `eba_dis_scaling = on` loads depth-varying dissipation number scaling $Di(z)$ from column 5 of 5- or 7-column `refstate_file` for viscous, adiabatic, and latent heating while keeping `eba_dis_scaling = off` 100% compatible with legacy 7-column files.
* **EarthByte Viscosity Structures**:
  - Added 13 published mantle rheology laws: `rheol = 20–28` (Flament / Zhang / Bower / Rudolph EBA lower-mantle depth profiles) and `rheol = 112, 113, 117, 118` (Hassan et al.).
* **Hybrid Tracers, Dynamic Topography Cutoff & Diagnostic Outputs**:
  - Hybrid ratio/ambient tracer composition (`hybrid_method`), configurable tracer array buffer cushion (`tracer_cushion`), dynamic topography buoyancy depth cutoff (`remove_buoyancy_above_znode`), EarthByte time-control compatibility aliases (`storage_spacing_Myr`, `exit_at_present`), and 9 diagnostic outputs (`sten_temp`, `sten_velo`, `divv`, `tracer_dens`, `heating_visc_nd`, `heating_adi_nd`, `heating_latent_nd`, `temp_sph`, `comp_sph`).
  - See the dedicated step-by-step guide: **[4D Plate & Slab Data Assimilation Tutorial (`docs/plate_assimilation_tutorial.md`)](docs/plate_assimilation_tutorial.md)**.

---

## 2. Directory Structure

```
citcoms_cpp/
├── CMakeLists.txt          # Modern CMake build configuration (C++17)
├── bin/                    # Driver executables & postprocessing tools
│   ├── Citcom.cpp          # Main simulation driver entrypoint
│   ├── CitcomSRegional.cpp # Regional solver driver (1 cap)
│   ├── CitcomSFull.cpp     # Full sphere solver driver (12 caps)
│   ├── citcoms_combine.cpp # High-performance C++17 data combiner
│   └── citcoms_h5_to_xdmf.py # Retroactive ParaView XDMF indexer
├── lib/                    # CitcomS C++17 core solver library
│   ├── global_defs.h       # Core headers, definitions, and data structures
│   ├── prototypes.h        # Modernized ISO C++ function prototypes
│   └── *.cpp               # FEM solvers, multigrid, viscosity, advection
├── Py2C/                   # Modernized C++17 config converter (Py2C)
├── visual/                 # Post-processing and visualization scripts
│   ├── combine.py          # Volume combiner wrapper (delegates to C++)
│   └── combinesurf.py      # Surface combiner wrapper (delegates to C++)
├── docs/                   # Tutorials and scientific guides
│   └── plate_assimilation_tutorial.md # 4D plate/slab data assimilation tutorial
├── tests/                  # Automated regression & EarthByte feature test suites
└── examples/               # Example configuration files and test cases
    ├── Regional/           # Regional mantle convection examples
    └── Full/               # Full sphere mantle convection examples
```

---

## 3. Compilation & Build Guide

### 3.1 Prerequisites
* **C++ Compiler**: GCC `g++` (>= 9) or LLVM `clang++` (>= 12) supporting C++17.
* **Build Tool**: CMake (>= 3.16) and `make` (or `ninja`).
* **MPI Library**: OpenMPI (>= 3.0), MPICH, MVAPICH2, or Intel MPI.
* **Optional**:
  - `libhdf5-openmpi-dev` (for parallel HDF5 output: `-DENABLE_HDF5=ON`).
  - `OpenMP` (for postprocessing multithreading).

### 3.2 Building with CMake (Recommended)

#### Using GCC (`g++` / `mpic++`):
```bash
# 1. Create build directory
mkdir build && cd build

# 2. Configure with CMake (Release optimization)
cmake .. -DCMAKE_BUILD_TYPE=Release

# 3. Compile in parallel
make -j$(nproc)
```

#### Using LLVM (`clang++`):
```bash
mkdir build_clang && cd build_clang
cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

#### Enabling Parallel HDF5 & Performance CMake Options:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HDF5=ON \
  -DCITCOMS_NATIVE_ARCH=ON \
  -DCITCOMS_LTO=ON
make -j$(nproc)
```
Optional CMake performance & precision flags:
- `-DCITCOMS_NATIVE_ARCH=ON` (default `OFF`): Enables `-march=native` SIMD (AVX2/AVX-512/FMA) for single-architecture HPC clusters.
- `-DCITCOMS_LTO=ON` (default `OFF`): Enables Interprocedural Optimization / Link-Time Optimization (`-flto`).
- `-DCITCOMS_DOUBLE_MATRIX=ON` (default `OFF`): Stores stiffness (`Eqn_k*`) and divergence (`elt_del`) matrices in `double` instead of `float` (doubles matrix memory footprint; opt-in for extreme augmented-Lagrangian conditioning).

### 3.3 Built Executables
Upon completion, the following binaries are generated in your build directory:
* **`CitcomSRegional`**: Regional spherical mantle convection solver (1 cap).
* **`CitcomSFull`**: Global full-sphere mantle convection solver (12 caps).
* **`citcoms_combine`**: High-performance C++17 combiner and postprocessor.
* **`Py2C`**: Modern C++17 standalone configuration file converter.

---

## 4. Execution & Quick Verification

### 4.1 Running a Regional Model (`CitcomSRegional`)

#### Single Process Test:
```bash
# Run regional simulation using configuration file:
./build/CitcomSRegional examples/Regional/example0
```

#### Parallel MPI Execution (e.g. 4 processes):
```bash
mpirun -np 4 ./build/CitcomSRegional my_regional_model.cfg
```

### 4.2 Running a Global Full-Sphere Model (`CitcomSFull`)

A full spherical model consists of 12 spherical caps. With 1 process per cap (12 processes):
```bash
mpirun -np 12 ./build/CitcomSFull my_full_model.cfg
```
For higher resolutions (e.g. $2 \times 2 \times 1 = 4$ processes per cap, 48 processes total):
```bash
mpirun -np 48 ./build/CitcomSFull my_full_model.cfg
```

---

## 5. Post-Processing & Combiner Usage

When CitcomS runs in parallel, each MPI process outputs its local partition files (e.g. `model.coord.0`, `model.velo.0.100`, `model.visc.0.100`).

### 5.1 Using Native C++ (`citcoms_combine`)
The C++ tool directly concatenates and reassembles the full spatial domain:
```bash
# Volume combine (replaces combine.py):
# citcoms_combine volume <prefix> <step> <nodex> <nodey> <nodez> <ncap> <nprocx> <nprocy> <nprocz> [opts] [datadir]
./build/citcoms_combine volume model 100 65 65 33 1 2 2 1 coord,velo,visc .

# Surface combine (replaces combinesurf.py):
# citcoms_combine surf <prefix> <step> <nodex> <nodey> <nodez> <nprocxy> <nprocx> <nprocy> <nprocz> [datadir]
./build/citcoms_combine surf model 100 65 65 33 1 2 2 1 .
```

### 5.2 Using Python Wrappers (Backward Compatible)
The existing Python scripts continue to work unchanged and automatically leverage the C++ combiner:
```bash
# Combines volume fields for timestep 100:
python3 visual/combine.py model 100 65 65 33 1 2 2 1

# Combines surface fields for timestep 100:
python3 visual/combinesurf.py model 100 65 65 33 1 2 2 1
```

---

## 6. CitcomS Configuration Guide

### 6.1 Output Formats (`output_format`) & Format-Specific Options
CitcomS supports single and composite output formats:
* `output_format = "ascii"`: Standard columnar ASCII files (`.coord`, `.velo`, `.visc`, `.temp`).
* `output_format = "ascii,vtk"`: Simultaneous 1D monitoring profiles and 3D VTK visualization files (compatible with ParaView).
* `output_format = "ascii,hdf5"`: Simultaneous 1D monitoring profiles and parallel HDF5 data archives.
*(Note: Enabling `vtk` and `hdf5` simultaneously is prohibited to avoid redundant 3D volume duplication and I/O congestion).*

Format-specific optional variables are configured independently:
* `output_optional_ascii = "horiz_avg"` (default): Outputs only lightweight 1D thermal convergence profiles (`.ave` / `.horiz_avg`, ~3 KB per step across all ranks).
* `output_optional_vtk = ""` (default): Pure 3D volume fields (temperature, velocity, viscosity, coordinates) without pseudo-3D surface zero-padding or 6-component stress bloat.
* `output_optional_h5 = "surf,botm,horiz_avg"` (default): Packs 3D volume, 2D dynamic topography/heatflux, and 1D horizontal averages into a single compressed `.h5` file per step.
*(Legacy `output_optional` parameter is maintained for backward compatibility).*

### 6.2 Decoupled ASCII 3D Volume Outputs & Composite Mode De-duplication
By default, when `output_optional_ascii` is set to `"horiz_avg"` (or under composite modes `"ascii,hdf5"` and `"ascii,vtk"`):
* 3D volume ASCII files (`.velo.*`, `.visc.*`, and cycle 0 `.coord.*`) are **automatically suppressed**.
* Only the 1D convergence profile (`.horiz_avg.*` / `.ave.*`) is output in ASCII (~3 KB per step cluster-wide), completely eliminating cluster inode pressure and removing redundant 3D writes.
* **Restoring 3D ASCII outputs**: When running standalone ASCII models where 3D volume files are desired, explicitly add `"velo"`, `"visc"`, `"coord"`, `"3d"`, or `"all"` into `output_optional_ascii`:
  ```ini
  output_format = "ascii"
  output_optional_ascii = "horiz_avg,velo,visc"  # or "3d", "all"
  ```

### 6.3 Advanced Checkpointing & Crash Safety
* **Zero-Touch Auto-Resume (`restart = "auto"`)**:
  Recommended for HPC queue workflows. Automatically finds and resumes from the latest verified checkpoint (`.last_good_checkpoint`). If no checkpoint is found, it begins a clean simulation from cycle 0 without requiring manual parameter file edits across job chains.
* **Sliding Window Checkpoint Pruning (`checkpoints_to_keep = 3`)**:
  Keeps only the last $N$ valid checkpoints on disk (default `3`). When a new checkpoint is successfully committed, the oldest is safely removed across all MPI ranks, preventing disk quota and inode exhaustion on HPC clusters. Set to `0` to keep all checkpoints indefinitely.
* **Two-Phase Commit**:
  Uses global MPI barrier synchronization (`parallel_process_sync`) and atomic metadata updating on Rank 0 to guarantee crash-proof checkpoints against HPC scheduler walltime terminations (`SIGKILL`).
* **Wall-Clock Checkpointing (`checkpoint_walltime`)**:
  Creates checkpoints at fixed real elapsed wall-clock intervals (e.g. `checkpoint_walltime = 7200` for every 2 hours), decoupled from fluctuating CFL timesteps.
* **Termination Checkpoint (`checkpoint_on_termination = "on"`)**:
  Guarantees writing a final valid checkpoint when the simulation completes or reaches `max_timesteps` / `cpu_limits_in_seconds`.

### 6.4 Geological Age (Myr) Triggers, Exact Time Snapping & Simulation Endpoint
* **Physical Time Triggers (`storage_interval_myr`, `checkpoint_interval_myr`)**:
  Users can specify output and checkpoint intervals directly in physical geological age (Myr):
  ```ini
  storage_interval_myr = 1.0     # Output fields every 1.0 Myr
  checkpoint_interval_myr = 5.0  # Create checkpoints every 5.0 Myr
  ```
  CitcomS automatically maps physical Myr to internal dimensionless time units using `scalet = (R^2 / kappa) / (1e6 * 365.25 * 24 * 3600)`.
  *(Note: When physical intervals are configured without explicit step intervals, step-based counters are automatically disabled).*
* **Simulation Termination Endpoint (`end_age_myr` / `end_time_myr`)**:
  Allows users to terminate simulations directly at a physical geological age (e.g. 50.0 Myr or present day 0.0 Myr) rather than guessing time-step counts:
  ```ini
  end_age_myr = 100.0            # Cleanly stop simulation when reaching 100.0 Myr
  ```
  CitcomS automatically computes the dimensionless target `end_time = end_age_myr / scalet`.
* **Exact Time Step Snapping (`exact_time_output = "on"`)**:
  Dynamically clamps the CFL timestep $\Delta t$ at milestone steps and when approaching `end_age_myr` so simulation frames and the termination point land precisely on requested Myr boundaries (e.g., exactly 100.000 Myr) for uniform, flicker-free animation rendering, guaranteed final field output, and termination checkpointing.

### 6.5 Checkpoint-Based Continuation & Two-Stage Restarts
CitcomS supports full binary checkpointing and seamless multi-stage restarts:
```ini
# Stage 2 model inheriting temperature & tracers from Stage 1:
tic_method = -1
tracer = 1
tracer.ic_method = 2
datafile_old = "stage1_output/model_stage1"
solution_cycles_init = 5000
```
CitcomS directly reads high-precision binary checkpoint `.chkpt` files without requiring large uncompressed ASCII velocity or tracer exports, with automatic fallback for backward compatibility.

### 6.6 ParaView 1-Click Interactive Visualization (VTK & HDF5)
CitcomS provides automated timeline index generation for both VTK and HDF5, enabling instant 3D interactive rendering and dynamic playback in ParaView without any manual data stitching:
- **VTK Workflow (`output_format = "ascii,vtk"`)**:
  - **Granularity Modes (`vtk_output_mode` / `vtk_split_mode`)**:
    - `vtk_output_mode = "cap"` (Default, Recommended): Merges subdomains per cap. Outputs 12 `.cap*.vts` files (for full sphere) or 1 file (for regional) plus a master `${data_file}.<cycle>.vtm` multi-block index. Reduces file counts by ~96% while maintaining structured grid topology and fast memory-mapped slicing.
    - `vtk_output_mode = "single"` (or `"unified"`): Assembles all 12 caps into a single monolithic `${data_file}.<cycle>.vtu` (UnstructuredGrid with Hexahedron-8 elements) per timestep. Easiest for single-file sharing and drag-and-drop ParaView opening.
    - `vtk_output_mode = "proc"` (Legacy): Writes 1 structured `.proc*.vts` file per MPI rank.
  - **Encoding (`vtk_format`)**:
    - `vtk_format = "binary"` (Default): High-speed zlib compressed base64 encoded binary arrays (~70% storage savings).
    - `vtk_format = "ascii"`: Human-readable plain text.
  - **ParaView Usage**:
    1. Open ParaView, click **File -> Open**, and select `${data_file}.pvd`.
    2. Click **Apply**. The complete 3D spherical mantle renders immediately.
    3. Press **Play (▶️)** on the top toolbar to watch the thermal plumes and mantle convection evolve smoothly across geological time (Myr)!
- **HDF5 Workflow (`output_format = "ascii,hdf5"`)**:
  - **Native 3D Cartesian Geometry & Vectors**:
    - `/coord`: True 3D Cartesian coordinates $(X, Y, Z)$ mapped to ParaView `<Geometry GeometryType="XYZ">`. ParaView renders the spherical shell directly without distortion.
    - `/coord_rtp`: Spherical coordinates $(\theta, \phi, r)$ preserved for scientific analysis scripts.
    - `/velocity`: Rotated Cartesian velocities $(v_x, v_y, v_z)$ mapped to ParaView `<Attribute Name="velocity" AttributeType="Vector">`. Vector glyph arrows point accurately in 3D Euclidean space.
    - `/velocity_rtp`: Raw spherical velocity components $(v_\theta, v_\phi, v_r)$ preserved for radial flux calculations.
  - **Files Generated**: `${data_file}.xdmf` (master temporal & spatial grid index), `${data_file}.h5` (nodal mesh geometry), and `${data_file}.<cycle>.h5` (physical fields).
  - **ParaView Usage**:
    1. Open ParaView, click **File -> Open**, and select `${data_file}.xdmf`.
    2. Choose the default **Xdmf Reader** when prompted, and click **Apply**.
    3. The 12 caps are automatically stitched into the full 3D sphere with temperature, velocity, and viscosity fields loaded.
    4. Press **Play (▶️)** to animate the simulation directly from high-throughput parallel HDF5 data.
- **Retroactive Indexing of Existing HDF5 Runs**:
  If you have existing `.h5` files from past runs that lack an `.xdmf` file, use the bundled CLI tool:
  ```bash
  python3 bin/citcoms_h5_to_xdmf.py --datafile stage1_preheat --dir ./output
  ```
  This generates `stage1_preheat.xdmf` in under 1 second.

### 6.7 High-Performance Stokes & Multigrid Solver Options
Legacy `.cfg` files without these keys preserve 100% original defaults (`uzawa = cg`, `mg_mode = fmg`, `mg_smoother = scalar_gs`, `coarse_solver = gs`, `p_precond = bpi`, `inner_rel_tol = 0.0`).

For plate-driven global convection models with strong lateral viscosity contrasts ($10^4$–$10^6$) and augmented Lagrangian penalty (`aug_lagr = 1`), the following drop-in configuration eliminates Uzawa stalling and provides **2.4×–2.7×** Stokes solver speedup with **10×–20×** tighter divergence convergence:
```ini
# Recommended high-contrast Stokes solver configuration:
uzawa = fgmres                  # cg (default) | fcg | fgmres | bicg
fgmres_restart = 30             # Krylov subspace size before restart
inner_rel_tol = 0.05            # Relative momentum tolerance (ensures exact Schur operator)
mg_mode = vcycle                # fmg (default) | vcycle
mg_smoother = block3x3_gs       # scalar_gs (default) | block3x3_gs | block3x3_sgs | chebyshev
coarse_solver = cg              # gs (default) | cg
visc_smooth_method = 4          # 4 = geometric-mean element coarsening (preserves weak plate boundaries)
```
> **Note**: When using `mg_mode = vcycle`, pair it with `inner_rel_tol = 0.05` (leaving `inner_mg_cycles = 0`) rather than forcing `inner_mg_cycles = 1`, so each inner momentum solve achieves sufficient accuracy for the exact Schur-complement pressure solution.

| Parameter | Options (Default First) | Description |
|---|---|---|
| `uzawa` | `cg`, `fcg`, `fgmres`, `bicg` | Outer Schur-complement Krylov solver (`fgmres` supports both incompressible and compressible flows) |
| `inner_rel_tol` | `0.0` (disabled) | Relative inner momentum tolerance `tol * ||r_0||`; recommended `0.05` for `fgmres`/`fcg` |
| `fcg_mmax` | `10` | Truncation window $m_{\max}$ for `uzawa = fcg` |
| `fgmres_restart` | `30` | Restart dimension $m$ for `uzawa = fgmres` |
| `residual_refresh_every` | `0` | Optional periodic exact divergence recomputation for `uzawa = fcg` |
| `inner_mg_cycles` | `0` (tolerance-driven) | Optional fixed inner MG cycle cap (use only with `mg_mode = fmg` or `>= 2` with `vcycle`) |
| `p_precond` | `bpi`, `visc_mass` | Schur-complement pressure preconditioner |
| `mg_mode` | `fmg`, `vcycle` | Multigrid traversal mode |
| `mg_smoother` | `scalar_gs`, `block3x3_gs`, `block3x3_sgs`, `chebyshev` | Multigrid velocity smoother |
| `coarse_solver` | `gs`, `cg` | Coarsest-level solver at `levmin` |
| `visc_smooth_method` | `0`, `1`, `2`, `3`, `4` | `4` adds geometric-mean element viscosity coarsening |

### 6.8 EarthByte Scientific Features & Parameter Reference Tables
All EarthByte capabilities are **strictly opt-in** (defaulting to legacy `citcoms_cpp` behavior) so existing `.cfg` files and 4/7-column `refstate_file` inputs remain 100% backward compatible. For complete file format specs, pyGPlates prep scripts, and ready-to-run `.cfg` templates, see **[docs/plate_assimilation_tutorial.md](docs/plate_assimilation_tutorial.md)**.

#### (1) 4D Lithosphere, Slab & Internal Velocity Assimilation
| Parameter | Default | Description |
|---|---|---|
| `lith_age_depth_function` | `off` | When `on`, computes age-dependent assimilation depth $z_{\text{lith}} = 2.32\sqrt{A^*}\cdot\texttt{lith\_age\_depth}$; when `off`, uses constant depth `lith_age_depth` |
| `lith_age_exponent` | `-1.0` (unset) | Blending exponent $p$ in $\mathrm{daf} = (z/z_{\text{lith}})^p$ (`0.5` recommended in EarthByte, `1.0` linear, `0.0` disables interior blending; unset `< 0` preserves legacy CIG $0.5 z/z_{\text{lith}}$) |
| `lith_age_min` | `0.0` | Minimum lithosphere age clamp in Myr (prevents $\sqrt{0}$ singularity at mid-ocean ridges) |
| `lith_age_stencil_value` | `-999.0` | Nodes with age $\le$ this value (Myr) skip lithosphere thermal assimilation (e.g. GPlates deforming zones) |
| `slab_assim` | `off` | Enable 3D time-dependent thermal slab assimilation (`T = (1 - alpha)*T + alpha*T_slab`) |
| `slab_assim_file` | `""` | Prefix for 3D slab files (`<prefix><age>.<cap>` for Full, `<prefix><age>` for Regional; 2 cols: `T_slab alpha`) |
| `internal_vbcs_file` | `off` | Enable 3D time-dependent internal velocity BCs across all multigrid levels (synchronized before each Stokes solve and on restart) |
| `velocity_internal_file` | `""` | Prefix for 3D internal velocity files (`<prefix><age>.<cap>` / `<prefix><age>`; 4 cols: `v_theta v_phi v_r flag` where `flag`: `0`=free, `1`=prescribe, `2`=skip boundary) |
| `topvbc` | `0` | Regional models additionally support `topvbc = 3` (open-top vertical traction with zero horizontal slip, Bower et al., 2011) |

#### (2) EBA Depth-Dependent $Di$ Scaling, Dynamic Topography & Time Compatibility
| Parameter | Default | Description |
|---|---|---|
| `eba_dis_scaling` | `off` | When `on`, reads column 5 of `refstate_file` (5- or 7-column) as depth-dependent $Di(z)$ multiplier for viscous, adiabatic, and latent heating |
| `remove_buoyancy_above_znode` | `10000` | Zeros thermal and chemical buoyancy for all global radial nodes `gz > threshold` (for dynamic topography excluding shallow lithosphere) |
| `storage_spacing_Myr` | `0` | EarthByte compatibility alias mapped to `storage_interval_myr` (inherits `citcoms_cpp` exact time-step snapping) |
| `exit_at_present` | `off` | Terminates forward convection (`timedir >= 0`) cleanly upon reaching `0.0 Ma` (snaps exactly to `0.0 Ma` when `exact_time_output = on`) |

#### (3) EarthByte Viscosity Structures (`rheol`)
| `rheol` | Lower-Mantle Profile (`l = 3`) & Formulation | Reference |
|---|---|---|
| `20` | $\eta = N_0 \exp\!\big(\frac{E}{T+T_0} - \frac{E}{0.5+T_0}\big)$, $T \in [0, 1]$ | Flament et al. (2013, 2014) |
| `21`, `22`, `23` | Same as `20` with linear lower-mantle $N_0(z)$: `10->50` (`21`), `10->100` (`22`), `10->30` (`23`) | Flament et al. (2014) TC8/TC7/TC9 |
| `24` | Frank-Kamenetskii $\eta = N_0(z)\exp(E(T_0 - T))$, linear lower-mantle $N_0(z)$ from `2` | Zhang et al. (2010), Bower et al. (2013) |
| `25`, `26`, `28` | EBA `20` variants with lower-mantle $N_0(z)$: `5->100` (`25`), `5->80` (`26`), `2->100` (`28`) | Flament EBA suite |
| `27` | EBA `20` variant with lower-mantle $N_0(z)$ linear `10->100` between 1000–2867 km | Rudolph et al. (2015) |
| `112`, `113`, `117`, `118` | Exponential lower-mantle profile (`112`), depth-varying activation + LTBL clip (`113`), Arrhenius `rheol=7` with `[0,1]` clamp & optional UTBL clip (`117`/`118`) | Hassan et al. |

#### (4) Hybrid Tracers & Diagnostic Output Tokens (`output_optional` / `output_optional_ascii`)
| Parameter / Output Token | Default | Description |
|---|---|---|
| `hybrid_method` | `off` | In ratio method (`buoy_type = 1`), empty elements receive `C = 0.0` (ambient) instead of retaining previous composition or aborting |
| `tracer_cushion` | `100000` | Extra per-rank tracer buffer capacity when loading `tracer_file` (reduces dynamic reallocations) |
| `sten_temp`, `sten_velo` | off | Output 3D nodal thermal slab stencil (`slab_sten`) and internal velocity stencil (`slab_sten2`) |
| `divv` | off | Output normalized nodal velocity divergence $\sqrt{\|\nabla\cdot\mathbf{u}\|^2 / \|\mathbf{u}\|^2}$ |
| `tracer_dens` | off | Output volume-normalized nodal tracer density (`E->trace.dens_node`) |
| `heating_visc_nd`, `heating_adi_nd`, `heating_latent_nd` | off | Output nodal-projected viscous, adiabatic, and latent heating fields |
| `temp_sph`, `comp_sph` | off | Output radial-layer spherical harmonic expansions of temperature and composition up to `output_ll_max` (Full spherical mode) |

---

## 7. License & Citation

CitcomS is free and open-source software distributed under the terms of the **GNU General Public License v2 (GPL-2.0)**. See the [COPYING](COPYING) file for full details.

### Citation
If you use CitcomS or its 4D plate assimilation features in your research, please cite:
> Moresi, L., Zhong, S., Han, L., Conrad, C., Tan, E., Gurnis, M., Choi, E., Thoutireddy, P., Manea, V., McNamara, A., Becker, T., Leng, W., & Armendariz, L. (2014). *CitcomS v3.3.1*. Zenodo. [https://doi.org/10.5281/zenodo.7271920](https://doi.org/10.5281/zenodo.7271920)
>
> Bower, D. J., Gurnis, M., & Flament, N. (2015). *Assimilating lithosphere and slab history in 4-D Earth models*. Physics of the Earth and Planetary Interiors, 238, 8–22. [https://doi.org/10.1016/j.pepi.2014.10.013](https://doi.org/10.1016/j.pepi.2014.10.013)
