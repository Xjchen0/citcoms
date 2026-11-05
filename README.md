# CitcomS: High-Performance Mantle Convection Simulation

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](COPYING)
[![CitcomS Manual](https://img.shields.io/badge/docs-CitcomS%20Manual-green.svg)](https://github.com/geodynamics/citcoms/blob/master/doc/citcoms-manual.pdf)

CitcomS is a finite element code designed to solve thermochemical mantle convection problems within the Earth's mantle in both regional (1 cap) and full spherical shell (12 caps) geometries.

---

## CitcomS Release Notes: High-Performance I/O, Advanced Checkpointing & Continuation Modernization

### 1. Purpose & Motivation
This release addresses long-standing I/O bottlenecks, storage inefficiencies, and runtime stability issues in CitcomS:

* **Dual-Format Output**: Eliminates the 20-year mutual exclusivity between ASCII, HDF5, and VTK. CitcomS now natively supports composite formats `"ascii,hdf5"` and `"ascii,vtk"`.
* **ASCII-Independent Continuation & Two-Stage Restart**: Completely decouples simulation continuation from fragmented ASCII files. Two-stage modeling (e.g., thermal spin-up in Stage 1 passing steady-state temperature to Stage 2 via `tic_method = -1`, and particle continuation via `tracer.ic_method = 2`) now directly reads high-precision binary checkpoint (`.chkpt`) files, with automatic backwards compatibility fallback to legacy `.velo`/`.tracer`.
* **Sliding Window Checkpoint Pruning (`checkpoints_to_keep`)**: Caps disk footprint by retaining only the last $N$ valid checkpoints (default `3`). Older checkpoint files across all MPI ranks are automatically removed, preventing disk quota and inode exhaustion on HPC clusters during long-running simulations.
* **Two-Phase Commit & Crash Safety (`last_good_checkpoint`)**: Uses a global MPI barrier (`parallel_process_sync`) to guarantee that all ranks finish writing before updating the valid checkpoint record. Protects against corrupted restarts when HPC schedulers terminate jobs via `SIGKILL` at walltime expiration.
* **Zero-Touch Auto-Resume (`restart = auto`)**: Enables hands-off Slurm/PBS job array chaining. Automatically resumes from the latest verified checkpoint if one exists, or starts a fresh simulation from cycle 0 without manually editing `.cfg` parameter files.
* **Geological Age & Physical Time Triggers (`storage_interval_myr`, `checkpoint_interval_myr`)**: Users can specify output and checkpoint intervals directly in physical geological age (Myr). CitcomS automatically converts them to internal dimensionless units via `scalet = (R^2 / kappa) / (1e6 * 365.25 * 24 * 3600)`.
* **Exact Time Step Snapping (`exact_time_output = on`)**: Dynamically clamps $\Delta t$ at milestone steps so simulation frames land precisely on requested Myr boundaries (e.g., 1.000 Myr, 2.000 Myr) for uniform, flicker-free animation rendering.
* **Decoupled Wall-Clock Checkpointing (`checkpoint_walltime`) & Termination Checkpoint (`checkpoint_on_termination`)**: Checkpoint based on real execution wall-clock time (e.g. every 2 hours) regardless of fluctuating CFL time steps, and guarantee a final clean checkpoint before exiting upon reaching `max_timesteps` or `cpu_limits_in_seconds`.
* **Inode Governance & Format-Specific Defaults**: Establishes minimal sensible defaults tailored for each format to eliminate HPC inode exhaustion while preserving essential 1D thermal convergence profiles (`.ave`).
* **Mutual Exclusion between VTK and HDF5**: Strictly prohibits enabling `vtk` and `hdf5` simultaneously to avoid redundant 3D volume duplication and I/O congestion.
* **HDF5 Stability & Array Bounds Fixes**: Resolves critical memory out-of-bounds indexing and longitudinal coordinate misalignment in HDF5 slice exports (surface/bottom heat flux and dynamic topography).
* **ParaView PVD Collection Indexing (`.pvd`)**: Master Rank 0 automatically generates and updates `${data_file}.pvd` referencing each output's multi-block `.vtm` or `.pvts` with accurate geological age timestamps (Myr), enabling seamless 1-click loading and animated timeline playback in ParaView.
* **HDF5 XDMF Curvilinear Grid Indexing (`.xdmf`)**: Bridges the gap between CitcomS HDF5 outputs and ParaView. CitcomS generates an accompanying `${data_file}.xdmf` descriptor that maps global 12-cap structured meshes (`3DSMesh`) and hyperslab slices to HDF5 arrays, unlocking direct 1-click 3D visualization and time playback in ParaView without any Python conversion scripts.
* **0-Dependency Python 2.6/2.7 & 3.x Post-Processing**: Upgrades legacy visual tools and introduces a pure standard-library archive/governance utility (`citcoms_archive.py`) requiring zero third-party packages (no numpy/scipy/pip required).

---

### 2. Modified Files & Components
* **C Solver & Checkpoint Subsystem (`lib/`, `bin/`)**:
  - `lib/checkpoints.h` & `lib/Checkpoints.c`: Implemented direct checkpoint extraction routines (`read_tic_from_checkpoint`, `get_checkpoint_elapsed_time`, and `read_tracers_from_checkpoint`), sliding window pruning (`prune_old_checkpoints`), two-phase commit with MPI barrier, auto-restart checkpoint detection (`determine_auto_restart_checkpoint`), and next output scheduling (`update_next_output_times`).
  - `lib/global_defs.h`: Added variables for composite output, sliding window pruning, crash safety, auto-restart, geological age intervals (Myr), dimensionless time intervals, wall-clock intervals, and exact time snapping.
  - `lib/Instructions.c`: Added composite output routing (`dual_output_ascii_h5`, `dual_output_ascii_vtk`), format-specific option buffers, re-entrant option parsing (`output_parse_optional`), strict VTK/HDF5 mutual exclusion, and parser support for `restart = auto`, `storage_interval_myr`, `checkpoint_interval_myr`, `checkpoint_walltime`, `checkpoints_to_keep`, `checkpoint_on_termination`, and `exact_time_output`.
  - `lib/Advection_diffusion.c`: Implemented exact time step snapping in `std_timestep` to dynamically clamp $\Delta t$ to milestone output and checkpoint times.
  - `bin/Citcom.c`: Integrated multi-criteria output and checkpoint triggering (step-based, model time/geological age, wall-clock time), tracking of last saved checkpoint, and final checkpoint on termination.
  - `lib/Initial_temperature.c`: Updated `read_tic_from_file` (`tic_method = -1`) to load initial temperature directly from `.chkpt` with fallback to `.velo`.
  - `lib/Problem_related.c`: Updated `get_initial_elapsed_time` to load initial elapsed time from `.chkpt` with fallback to `.velo`.
  - `lib/Tracer_setup.c`: Updated `read_old_tracer_file` (`tracer.ic_method = 2`) to restore particles and composition directly from `.chkpt` with fallback to `.tracer`.
  - `lib/Output.c`: Defined format-specific sensible defaults and enabled VTK format parsing under composite modes.
  - `lib/Output_vtk.c`: Implemented `write_pvd` for continuous ParaView collection generation, supporting true geological physical time (Myr) and automatic history truncation on simulation rewind/restart.
  - `lib/Output_h5.c`: Fixed 2D nodal indexing, corrected bottom topography pointer, unified conditional compilation guards (`#ifdef USE_HDF5`), and implemented `h5output_write_xdmf` mapping multi-cap curvilinear coordinates and fields directly to XDMF 2.0 temporal grids.
* **Post-Processing, Indexing & Inode Management (`visual/`, `bin/`, `tests/`)**:
  - `bin/citcoms_h5_to_xdmf.py`: Standalone CLI tool to retroactively inspect and index past CitcomS HDF5 simulations into ParaView-ready `.xdmf` files.
  - `visual/parser.py`, `visual/combine.py`, `visual/combinesurf.py`: Updated syntax, replaced deprecated builtins (`file()`), fixed integer division (`//`), and ensured 100% Python 2.6/2.7 and Python 3.x compatibility.
  - `visual/citcoms_archive.py`: New standalone zero-dependency archive and inode management utility using Python standard library only.
  - `tests/`: Added automated test suite (`test_pvd_and_xdmf.py`, `test_checkpoint_features.py`, `test_ascii_dedup.py`) covering PVD/XDMF XML validity, restart pruning, and 1D monitoring de-duplication.

---

### 3. Configuration & Parameter Guide
* **Output Format (`output_format`)**:
  - `output_format = "ascii,hdf5"` : Recommended for complete scientific workflows. Outputs lightweight 1D ASCII profile (`.ave`) for runtime monitoring alongside full-featured HDF5 3D/2D archive.
  - `output_format = "ascii,vtk"`  : Recommended for fast ParaView 3D visualization. Pair with `vtk_format = "binary"`.
  - *Note*: `"hdf5,vtk"` is prohibited and will abort with an error.

* **Format-Specific Optional Variables & De-Duplication**:
  - `output_optional_ascii = "horiz_avg"` (default): Strictly outputs only the lightweight 1D profile (`.ave` / `.horiz_avg`, ~3 KB per step across the cluster) for real-time monitoring. 3D volume ASCII files (`.velo`, `.visc`, `.coord`) are suppressed by default to prevent inode exhaustion, eliminating 3D volume duplication in composite modes (`"ascii,hdf5"` and `"ascii,vtk"`). To output 3D ASCII files, explicitly include `"velo"`, `"visc"`, or `"3d"` in `output_optional_ascii`.
  - `output_optional_vtk = ""` (default): Pure 3D volume output (temperature, velocity, viscosity, coordinates) without pseudo-3D surface zero-padding or 6-component stress bloat.
  - `output_optional_h5 = "surf,botm,horiz_avg"` (default): Packs 3D volume (temperature, velocity, viscosity), 2D dynamic topography/heatflux, and 1D horizontal averages into a single compressed `.h5` archive per step.
  *(Legacy `output_optional` parameter is maintained for backward compatibility).*

* **Auto-Resume & Checkpoint Management**:
  - `restart = "auto"` : Recommended for HPC queue workflows. Automatically finds and resumes from the latest verified checkpoint (`.last_good_checkpoint`). If no checkpoint is found, it begins a clean simulation from cycle 0.
  - `checkpoints_to_keep = 3` (default): Keeps only the last 3 valid checkpoints on disk. When a new checkpoint is committed, the oldest is safely removed across all ranks. Set to `0` to keep all checkpoints indefinitely.
  - `checkpoint_walltime = 7200` : Create checkpoints every 7200 seconds (2 hours) of wall-clock time, independent of timestep size.
  - `checkpoint_on_termination = "on"` (default): Guarantees saving a final checkpoint when the run exits normally or reaches `max_timesteps` / `cpu_limits_in_seconds`.

* **Geological Age (Myr) & Physical Time Triggers**:
  - `storage_interval_myr = 1.0` : Output simulation fields every 1.0 Myr. Automatically converted to dimensionless units using internal `scalet`.
  - `checkpoint_interval_myr = 5.0` : Create checkpoints every 5.0 Myr.
  - `storage_interval_time = 0.0` : Output fields at dimensionless time intervals (default `0.0` = off).
  - `exact_time_output = "on"` (default): Automatically clamps $\Delta t$ at milestone steps so field outputs and checkpoints land exactly on the requested Myr boundaries for uniform, flicker-free animation rendering.
  *(Note: If `storage_interval_myr` is specified without explicit `storage_spacing`, step-based output is automatically disabled).*

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

* **ParaView 1-Click Interactive Visualization (VTK & HDF5)**:
  CitcomS provides automated timeline index generation for both VTK and HDF5, enabling instant 3D interactive rendering and dynamic playback in ParaView without any manual data stitching:
  - **VTK Workflow (`output_format = "ascii,vtk"`)**:
    - **Files Generated**: `${data_file}.pvd` (master temporal collection), `${data_file}.<cycle>.vtm` (multi-block index), and `${data_file}.proc*.<cycle>.vts` (per-rank structured grids).
    - **ParaView Usage**:
      1. Open ParaView, click **File -> Open**, and select `${data_file}.pvd`.
      2. Click **Apply**. The complete 3D spherical mantle renders immediately.
      3. Press **Play (▶️)** on the top toolbar to watch the thermal plumes and mantle convection evolve smoothly across geological time (Myr)!
  - **HDF5 Workflow (`output_format = "ascii,hdf5"`)**:
    - **Files Generated**: `${data_file}.xdmf` (master temporal & spatial grid index), `${data_file}.h5` (nodal mesh geometry), and `${data_file}.<cycle>.h5` (physical fields).
    - **ParaView Usage**:
      1. Open ParaView, click **File -> Open**, and select `${data_file}.xdmf`.
      2. Choose the default **Xdmf Reader** when prompted, and click **Apply**.
      3. The 12 caps are automatically stitched into the full 3D sphere with temperature, velocity, and viscosity fields loaded.
      4. Press **Play (▶️)** to animate the simulation directly from high-throughput parallel HDF5 data.
  - **Retroactive Indexing of Existing HDF5 Runs**:
    If you have existing `.h5` files from past runs that lack an `.xdmf` file, use the bundled CLI tool:
    ```bash
    python bin/citcoms_h5_to_xdmf.py --datafile stage1_preheat --dir ./output
    ```
    This generates `stage1_preheat.xdmf` in under 1 second.

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
