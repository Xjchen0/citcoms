# CitcomS: High-Performance Mantle Convection Simulation

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](COPYING)
[![CitcomS Manual](https://img.shields.io/badge/docs-CitcomS%20Manual-green.svg)](https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf)

CitcomS is a finite element code designed to solve thermochemical mantle convection problems within the Earth's mantle in both regional (1 cap) and full spherical shell (12 caps) geometries.

---

## CitcomS Release Notes: High-Performance I/O & Continuation Modernization

### 1. Purpose & Motivation
This release addresses long-standing I/O bottlenecks, storage inefficiencies, and runtime stability issues in CitcomS:

* **Dual-Format Output**: Eliminates the 20-year mutual exclusivity between ASCII, HDF5, and VTK. CitcomS now natively supports composite formats `"ascii,hdf5"` and `"ascii,vtk"`.
* **ASCII-Independent Continuation & Two-Stage Restart**: Completely decouples simulation continuation from fragmented ASCII files. Two-stage modeling (e.g., thermal spin-up in Stage 1 passing steady-state temperature to Stage 2 via `tic_method = -1`, and particle continuation via `tracer.ic_method = 2`) now directly reads high-precision binary checkpoint (`.chkpt`) files, with automatic backwards compatibility fallback to legacy `.velo`/`.tracer`.
* **Inode Governance & Format-Specific Defaults**: Establishes minimal sensible defaults tailored for each format to eliminate HPC inode exhaustion while preserving essential 1D thermal convergence profiles (`.ave`).
* **Mutual Exclusion between VTK and HDF5**: Strictly prohibits enabling `vtk` and `hdf5` simultaneously to avoid redundant 3D volume duplication and I/O congestion.
* **HDF5 Stability & Array Bounds Fixes**: Resolves critical memory out-of-bounds indexing and longitudinal coordinate misalignment in HDF5 slice exports (surface/bottom heat flux and dynamic topography).
* **0-Dependency Python 2.6/2.7 & 3.x Post-Processing**: Upgrades legacy visual tools and introduces a pure standard-library archive/governance utility (`citcoms_archive.py`) requiring zero third-party packages (no numpy/scipy/pip required).

---

### 2. Modified Files & Components
* **C Solver & Checkpoint Subsystem (`lib/`)**:
  - `lib/checkpoints.h` & `lib/Checkpoints.c`: Implemented direct checkpoint extraction routines (`read_tic_from_checkpoint`, `get_checkpoint_elapsed_time`, and `read_tracers_from_checkpoint`) with mesh verification, seek offsets, and sentinel checks.
  - `lib/Initial_temperature.c`: Updated `read_tic_from_file` (`tic_method = -1`) to load initial temperature directly from `.chkpt` with fallback to `.velo`.
  - `lib/Problem_related.c`: Updated `get_initial_elapsed_time` to load initial elapsed time from `.chkpt` with fallback to `.velo`.
  - `lib/Tracer_setup.c`: Updated `read_old_tracer_file` (`tracer.ic_method = 2`) to restore particles and composition directly from `.chkpt` with fallback to `.tracer`.
  - `lib/global_defs.h` & `lib/Instructions.c`: Added composite output routing (`dual_output_ascii_h5`, `dual_output_ascii_vtk`), format-specific option buffers, re-entrant option parsing (`output_parse_optional`), and strict VTK/HDF5 mutual exclusion.
  - `lib/Output.c`: Defined format-specific sensible defaults and enabled VTK format parsing under composite modes.
  - `lib/Output_h5.c`: Fixed 2D nodal indexing (`p = 1 + i + j*nx`), corrected bottom topography pointer (`topo = E->slice.tpgb[1]`), corrected surface topography indexing, and added `#ifndef USE_HDF5` stubs.
* **Post-Processing & Inode Management (`visual/`)**:
  - `visual/parser.py`, `visual/combine.py`, `visual/combinesurf.py`: Updated syntax, replaced deprecated builtins (`file()`), fixed integer division (`//`), and ensured 100% Python 2.6/2.7 and Python 3.x compatibility.
  - `visual/citcoms_archive.py`: New standalone zero-dependency archive and inode management utility using Python standard library only.

---

### 3. Configuration & Parameter Guide
* **Output Format (`output_format`)**:
  - `output_format = "ascii,hdf5"` : Recommended for complete scientific workflows. Outputs lightweight 1D ASCII profile (`.ave`) for runtime monitoring alongside full-featured HDF5 3D/2D archive.
  - `output_format = "ascii,vtk"`  : Recommended for fast ParaView 3D visualization. Pair with `vtk_format = "binary"`.
  - *Note*: `"hdf5,vtk"` is prohibited and will abort with an error.

* **Format-Specific Optional Variables**:
  Mandatory baseline fields (grid coordinates, 3D temperature $T$, velocity vector $\vec{v}$, and viscosity $\eta$) are always generated automatically for all formats.
  Optional variables are configured independently per format with the following defaults:
  - `output_optional_ascii = "horiz_avg"` (default): Only outputs one lightweight (~3 KB) `.ave` file per timestep across all MPI ranks. Zero inode overhead.
  - `output_optional_vtk = ""` (default): Pure 3D volume output without pseudo-3D surface zero-padding or 6-component stress bloat.
  - `output_optional_h5 = "surf,botm,horiz_avg"` (default): Packs 3D volume, 2D dynamic topography/heatflux, and 1D horizontal averages into a single `.h5` file per step.
  *(Legacy `output_optional` parameter is maintained for backward compatibility).*

* **Two-Stage Simulation Restart (`tic_method = -1`, `tracer.ic_method = 2`)**:
  Stage 2 models inheriting temperature or tracers from Stage 1:
  ```ini
  tic_method = -1
  tracer = 1
  tracer.ic_method = 2
  datafile_old = "stage1_output/model_stage1"
  solution_cycles_init = 5000
  ```
  CitcomS will automatically read `stage1_output/model_stage1.chkpt.*.5000`. Stage 1 no longer needs to export large ASCII `.velo` or `.tracer` files!

* **Standalone Inode Governance Tool (`citcoms_archive.py`)**:
  Run on HPC clusters using system Python (no pip or virtualenv required):
  ```bash
  # Status check (display unarchived steps and file counts):
  python visual/citcoms_archive.py -c model.cfg --status

  # Batch archive all timesteps into .tar.gz and clean raw files:
  python visual/citcoms_archive.py -c model.cfg --all --clean

  # Lossless restore of an archived step when raw ASCII is needed:
  python visual/citcoms_archive.py -c model.cfg -s 5000 --extract
  ```

---

## 4. Compilation & Installation Guide

### 4.1 System Requirements & Prerequisites
CitcomS is written in C and parallelized using MPI. Ensure the following prerequisites are installed on your system:
* **Mandatory**:
  - **C Compiler**: GCC (>= 4.x), Intel oneAPI/ICC, or Clang.
  - **MPI Library**: OpenMPI, MPICH, MVAPICH2, or Intel MPI (providing `mpicc`).
  - **GNU Autotools**: `autoconf` (>= 2.59), `automake` (>= 1.9), and `libtool` (required if building directly from a Git checkout).
* **Optional**:
  - **Parallel HDF5 (PHDF5)**: Required only if compiling with HDF5 output.
    - Debian/Ubuntu: `apt-get install libhdf5-openmpi-dev` (or `libhdf5-mpich-dev`).
    - RHEL/CentOS/Rocky: `yum install hdf5-openmpi-devel` (or `hdf5-mpich-devel`).
  - **zlib (`libz`)**: Enables transparent gzip file compression (`USE_GZDIR`).
  - **GGRD / NetCDF / GMT**: Required only if `--with-ggrd` is specified.

### 4.2 Step-by-Step Build Procedure
1. **Generate the GNU build scripts** (for Git checkouts):
   ```bash
   autoreconf -i
   ```

2. **Configure the build environment**:
   * *Basic build* (MPI wrapper `mpicc` detected automatically):
     ```bash
     ./configure CC=mpicc
     ```
   * *Production optimization* (recommended for compute clusters):
     ```bash
     ./configure CC=mpicc CFLAGS="-O3 -march=native"
     ```
   * *Enable Parallel HDF5*:
     ```bash
     # If HDF5 is installed in a custom prefix:
     export PHDF5_HOME=/opt/phdf5
     ./configure CC=mpicc --with-hdf5
     ```
   * *Specifying custom include and library directories*:
     ```bash
     ./configure CC=mpicc \
       CPPFLAGS="-I/opt/custom/include" \
       LDFLAGS="-L/opt/custom/lib"
     ```

3. **Compile**:
   ```bash
   make -j$(nproc)
   ```

4. **Executable Binaries**:
   Upon successful compilation, two pure-C solver binaries are generated in `bin/`:
   - `bin/CitcomSRegional` : Regional spherical mantle convection solver (1 cap).
   - `bin/CitcomSFull`     : Global full spherical mantle convection solver (12 caps).

5. *(Optional) Install to a Destination Directory*:
   ```bash
   make install prefix=$HOME/citcoms
   ```

### 4.3 Running CitcomS & Quick Verification
Run a test simulation using the bundled examples:
* **Uniprocessor test**:
  ```bash
  ./bin/CitcomSRegional examples/Regional/example0
  ```
* **MPI Parallel Execution**:
  ```bash
  # Run a regional model across 4 MPI processes:
  mpirun -np 4 ./bin/CitcomSRegional my_regional_model.cfg

  # Run a full spherical model across 12 MPI processes (1 rank per cap):
  mpirun -np 12 ./bin/CitcomSFull my_global_model.cfg
  ```

---

## 5. Documentation, Help & Citations

See the file `INSTALL` and the Manual for building and installation instructions.

See the file `NEWS` for new features and bug fixes of this release.

The Manual can be downloaded at:
[https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf](https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf)

For bug reports, please open an issue.

For help, please post to the mantle convection forum:
[https://community.geodynamics.org/c/mantle-convection/12](https://community.geodynamics.org/c/mantle-convection/12)

### Citation
CitcomS v3.3.1 can be cited as:
> Moresi, Louis, Zhong, Shijie, Han, Lijie, Conrad, Clint, Tan, Eh, Gurnis, Michael, Choi, Eunseo, Thoutireddy, Pururav, Manea, Vlad, McNamara, Allen, Becker, Thorsten, Leng, Wei, & Armendariz, Luis. (2014). CitcomS v3.3.1 (v3.3.1). Zenodo. [https://doi.org/10.5281/zenodo.7271920](https://doi.org/10.5281/zenodo.7271920)

Full citation information can be found on the CIG website:
[https://geodynamics.org/resources/citcoms/howtocite](https://geodynamics.org/resources/citcoms/howtocite)

### License
CitcomS is free software published under the GPL2.0 license. See the file `COPYING` for copying conditions.
