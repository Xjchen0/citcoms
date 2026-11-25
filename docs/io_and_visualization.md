# I/O Formats & ParaView Interactive Visualization

← Back to [README](../README.md)

---

## 1. Output Formats (`output_format`) & Format-Specific Controls

`CitcomS` supports both standalone and composite output formats:

| `output_format` Value | Behavior |
|---|---|
| `ascii` | Standard columnar ASCII files (`.coord`, `.velo`, `.visc`, `.horiz_avg`, etc.) |
| `vtk` | VTK 3D outputs + ParaView `.pvd` timeline |
| `hdf5` | Parallel HDF5 archives + ParaView `.xdmf` timeline (requires `-DENABLE_HDF5=ON`) |
| `ascii,vtk` | Simultaneous lightweight 1D/diagnostic ASCII + 3D VTK visualization |
| `ascii,hdf5` | Simultaneous lightweight 1D/diagnostic ASCII + parallel HDF5 archives |

> Enabling `vtk` and `hdf5` simultaneously is disallowed to avoid duplicate 3D volume I/O congestion.

### Format-Specific Optional Variables
In composite modes (`ascii,vtk` or `ascii,hdf5`), optional fields can be configured independently per format:

```ini
output_format=ascii,vtk
output_optional_ascii=horiz_avg          # Lightweight 1D monitoring (~3 KB/step cluster-wide)
output_optional_vtk=comp_nd,surf         # Export composition1..N and surface topo into VTK
```

- In **composite modes**, `output_optional_ascii` defaults to `horiz_avg` so redundant 3D ASCII volume files (`.velo.*`, `.visc.*`, `.coord.*`) are automatically suppressed, avoiding inode exhaustion.
- In **standalone `ascii` or `hdf5` mode**, `output_optional` is honored directly (unless explicitly overridden by `output_optional_ascii` / `output_optional_h5`), ensuring 100% compatibility with legacy and EarthByte `.cfg` files.
- To restore full 3D ASCII output in standalone ASCII mode, include `velo`, `visc`, `coord`, `3d`, or `all`:
  ```ini
  output_format=ascii
  output_optional=horiz_avg,velo,visc
  ```

---

## 2. VTK Workflow & ParaView 1-Click Playback

Configure VTK granularity and encoding without quotes:

```ini
output_format=ascii,vtk
output_optional_ascii=horiz_avg
output_optional_vtk=comp_nd,surf
vtk_output_mode=cap        # cap (default) | single | proc
vtk_format=binary          # binary (default, zlib compressed) | ascii
```

### Granularity Modes (`vtk_output_mode` / `vtk_split_mode`)
- **`cap` (Default, Recommended)**: Merges MPI subdomains per spherical cap — writes 12 `.cap*.vts` files (Full) or 1 `.vts` file (Regional) plus a `${datafile}.<cycle>.vtm` multi-block index. Reduces file count by ~96% while preserving structured grid topology. Exports `temperature`, `velocity`, `viscosity`, plus optional `stress`, `composition1..N` (`comp_nd`), and `surf`.
- **`single` (aliases: `unified`, `monolithic`, `all`)**: Assembles the entire domain into one monolithic `${datafile}.<cycle>.vtu` per timestep for single-file sharing.
- **`proc` (Legacy)**: Writes one `.proc*.vts` file per MPI rank.

### Opening in ParaView
1. Rank 0 automatically maintains `${datafile}.pvd` mapping every output cycle to physical geological time (Myr) and auto-pruning superseded entries on checkpoint rewind.
2. Open `${datafile}.pvd` in ParaView → click **Apply** → press **Play (▶)** to animate 3D mantle evolution.

---

## 3. Parallel HDF5 & XDMF Workflow

```ini
output_format=ascii,hdf5
output_optional_ascii=horiz_avg
output_optional_h5=surf,botm,horiz_avg
```

### Native 3D Cartesian Datasets
- `/coord` & `/velocity`: True 3D Cartesian coordinates $(X, Y, Z)$ and rotated 3D vectors $(v_x, v_y, v_z)$ mapped directly to ParaView `<Geometry GeometryType="XYZ">` and vector glyphs.
- `/coord_rtp` & `/velocity_rtp`: Spherical coordinates $(\theta, \phi, r)$ and components $(v_\theta, v_\phi, v_r)$ preserved for radial analysis scripts.
- Automatically maintains `${datafile}.xdmf` linking `${datafile}.h5` (mesh) and `${datafile}.<cycle>.h5` (fields).

### Opening in ParaView & Retroactive Indexing
1. Open `${datafile}.xdmf` in ParaView with the default **Xdmf Reader** → click **Apply** → press **Play (▶)**.
2. To generate an `.xdmf` index retroactively for existing `.h5` runs:
   ```bash
   python3 bin/citcoms_h5_to_xdmf.py --datafile stage1_preheat --dir ./output
   ```
