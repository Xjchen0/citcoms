# Build, Toolchain, RDMA Optimization & Post-Processing

← Back to [README](../README.md)

---

## 1. C++17 Architecture & Compiler Support

- **ISO C++17 (`-std=c++17`)**: All legacy C/K&R routines (`lib/*.cpp`, `bin/*.cpp`, `Py2C/*.cpp`) are standardized with full ANSI/ISO prototypes, const-correctness, explicit typed allocations, and standard library helpers.
- **Verified Compilers**:
  - **GCC (`g++`)**: GCC 9 – 16
  - **LLVM (`clang++`)**: Clang 12 – 22
- **Prerequisites**:
  - CMake >= 3.16 and `make` (or `ninja`)
  - MPI library (OpenMPI >= 3.0, MPICH, MVAPICH2, or Intel MPI)
  - Optional: `libhdf5-openmpi-dev` (for `-DENABLE_HDF5=ON`), `zlib` (for VTK binary compression), and OpenMP (for multithreaded `citcoms_combine`)

---

## 2. Building with CMake

### Standard GCC Build
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### LLVM / Clang Build
```bash
mkdir -p build_clang && cd build_clang
cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Enabling HDF5 & Performance Options
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HDF5=ON \
  -DCITCOMS_NATIVE_ARCH=ON \
  -DCITCOMS_LTO=ON
make -j$(nproc)
```

| CMake Option | Default | Description |
|---|---|---|
| `ENABLE_HDF5` | `OFF` | Enable parallel HDF5 output (`output_format=hdf5` / `ascii,hdf5`) |
| `CITCOMS_NATIVE_ARCH` | `OFF` | Enable `-march=native` SIMD (AVX2 / AVX-512 / FMA) on homogeneous HPC nodes |
| `CITCOMS_LTO` | `OFF` | Enable Link-Time / Interprocedural Optimization (`-flto`) |
| `CITCOMS_DOUBLE_MATRIX` | `OFF` | Store stiffness (`Eqn_k*`) and divergence (`elt_del`) matrices in `double` instead of `float` (2× matrix memory; opt-in for extreme augmented-Lagrangian conditioning) |

### Built Executables
- `CitcomSRegional`: Regional spherical-cap solver (1 cap)
- `CitcomSFull`: Global full-sphere solver (12 caps)
- `citcoms_combine`: High-performance C++17 domain-reassembly combiner
- `Py2C`: Standalone C++17 Pyre-to-flat `.cfg` converter

---

## 3. Running Simulations

### Regional (1 Cap)
```bash
mpirun -np 4 ./build/CitcomSRegional my_regional_model.cfg
```

### Global Full Sphere (12 Caps)
```bash
# 1 rank per cap (12 ranks total)
mpirun -np 12 ./build/CitcomSFull my_full_model.cfg

# 2x2x1 ranks per cap (48 ranks total)
mpirun -np 48 ./build/CitcomSFull my_full_model.cfg
```

---

## 4. InfiniBand RDMA & UCX Communication Optimizations

- **Persistent Communication Buffer Pools**: Replaces per-iteration `malloc`/`free` inside multigrid smoothers and nodal field exchanges (`full_exchange_id_d`, `exchange_node_d`, `exchange_node_f`, `full_exchange_snode_f`) with persistent adaptive buffers. Stable virtual addresses prevent InfiniBand UCX Registration Cache (RCache) desynchronization (`Transport retry count exceeded, synd 0x15`).
- **Inverted Non-Blocking Ordering**: Posts all `MPI_Irecv` calls before packing and posting `MPI_Isend` calls across boundary exchanges and tracer migration, eliminating Receiver-Not-Ready (RNR NACK) retries and deadlocks on hundreds to thousands of MPI ranks.

---

## 5. High-Performance Post-Processing (`citcoms_combine`)

`citcoms_combine` delivers a **50×–100×** speedup over legacy Python scripts via buffered 64 KB C++17 streaming and OpenMP parallelism:

```bash
# Volume combine:
# citcoms_combine volume <prefix> <step> <nodex> <nodey> <nodez> <ncap> <nprocx> <nprocy> <nprocz> [opts] [datadir]
./build/citcoms_combine volume model 100 65 65 33 1 2 2 1 coord,velo,visc .

# Surface combine:
# citcoms_combine surf <prefix> <step> <nodex> <nodey> <nodez> <nprocxy> <nprocx> <nprocy> <nprocz> [datadir]
./build/citcoms_combine surf model 100 65 65 33 1 2 2 1 .
```

The existing Python wrappers (`visual/combine.py` and `visual/combinesurf.py`) remain 100% backward compatible and automatically delegate to `citcoms_combine` when available.
