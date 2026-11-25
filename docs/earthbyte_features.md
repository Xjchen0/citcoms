# EarthByte Scientific Features & Complete Parameter Reference

← Back to [README](../README.md) | Hands-on Tutorial: **[4D 板块与俯冲板片数据同化实战教程](plate_assimilation_tutorial.md)**

---

## Overview & Backward Compatibility

This engine absorbs the core scientific capabilities of **EarthByte/CitcomS** (`https://github.com/EarthByte/citcoms`, GPL-2.0) directly into the modernized `citcoms_cpp` C++17 codebase:
- All new features are **strictly opt-in** (defaulting to legacy `citcoms_cpp` behavior) so legacy `.cfg` files and 4/7-column `refstate_file` profiles produce identical numerics.
- All EarthByte features integrate seamlessly with `citcoms_cpp`'s **FGMRES + V-cycle + Block 3×3** Stokes solvers, **exact geological-age step snapping**, **zero-lag velocity BC synchronization** (updating internal velocity BCs *before* each Stokes solve and upon checkpoint restart), and **VTK/HDF5 ParaView output**.

---

## 1. 4D Lithosphere, Slab & Internal Velocity Assimilation Parameters

| Parameter | Default | Description |
|---|---|---|
| `lith_age_depth_function` | `off` | `on`: age-dependent assimilation depth $z_{\text{lith}} = 2.32\sqrt{A_{\text{nd}}}\cdot d_0$ ($d_0$ = `lith_age_depth`, $A_{\text{nd}} = A/\text{scalet}$); `off`: constant depth `lith_age_depth` |
| `lith_age_exponent` | `-1.0` (unset) | Blending exponent $p$ in $\mathrm{daf}=(z/z_{\text{lith}})^p$ (`0.5` EarthByte recommended, `1.0` linear, `0.0` disables interior blending; `< 0` preserves legacy CIG $0.5\,z/z_{\text{lith}}$) |
| `lith_age_min` | `0.0` | Minimum lithosphere age clamp in Myr (prevents $\sqrt{0}$ singularity at mid-ocean ridges) |
| `lith_age_stencil_value` | `-999.0` | Nodes with age ≤ this value (Myr) skip lithosphere thermal assimilation (e.g., GPlates deforming zones) |
| `slab_assim` | `off` | Enable 3D time-dependent thermal slab assimilation: $T = (1-\alpha)T + \alpha T_{\text{slab}}$ |
| `slab_assim_file` | `""` | Path prefix for 3D slab files (`<prefix><age>.<cap>` Full, `<prefix><age>` Regional; 2 columns: `T_slab alpha`) |
| `internal_vbcs_file` | `off` | Enable 3D time-dependent internal velocity BCs across all multigrid levels (synchronized before each Stokes solve and on restart) |
| `velocity_internal_file` | `""` | Path prefix for 3D internal velocity files (`<prefix><age>.<cap>` / `<prefix><age>`; 4 columns: `v_theta v_phi v_r flag`, where `flag`: `0`=free, `1`=prescribe, `2`=skip boundary) |
| `topvbc` | `0` | Regional models additionally support `topvbc=3` (open-top vertical traction with zero horizontal slip; Bower et al., 2011) |

---

## 2. EBA Depth-Dependent $Di(z)$ Scaling, Dynamic Topography & Time Aliases

| Parameter | Default | Description |
|---|---|---|
| `eba_dis_scaling` | `off` | When `on`, reads column 5 of `refstate_file` (5- or 7-column format) as depth-dependent $Di(z)$ multiplier for viscous, adiabatic, and latent heating. When `off`, sets $Di(z)=1.0$ for 100% compatibility with legacy 7-column files whose column 5 is `0.0` |
| `remove_buoyancy_above_znode` | `10000` | Zeros thermal and chemical buoyancy for all global radial nodes `gz > threshold` (for dynamic topography excluding shallow lithospheric buoyancy) |
| `storage_spacing_Myr` | `0` | EarthByte compatibility alias mapped to `storage_interval_myr` (inherits `citcoms_cpp` exact time-step snapping) |
| `exit_at_present` | `off` | Terminates forward convection (`timedir >= 0`) cleanly upon reaching `0.0 Ma` (snaps directly onto `0.0 Ma` when `exact_time_output=on`) |
| `reset_startage` | `0` | Resets `start_age` to `find_age_in_MY(E)` (fixed EarthByte forward-time calculation) |

---

## 3. EarthByte Mantle Viscosity Structures (`rheol`)

| `rheol` | Formulation & Lower-Mantle (`l = 3`) Profile | Reference |
|---|---|---|
| `20` | $\eta = N_0 \exp\big(\frac{E}{T+T_0} - \frac{E}{0.5+T_0}\big)$, $T \in [0,1]$ (reference $T=0.5$) | Flament et al. (2013, 2014) |
| `21` | Same as `20`, linear lower-mantle $N_0(z)$ from `10` to `50` (670–2867 km) | Flament et al. (2014) TC8 |
| `22` | Same as `20`, linear lower-mantle $N_0(z)$ from `10` to `100` (670–2867 km) | Flament et al. (2014) TC7 / Zhang et al. (2010) |
| `23` | Same as `20`, linear lower-mantle $N_0(z)$ from `10` to `30` (670–2867 km) | Flament et al. (2014) TC9 |
| `24` | Frank-Kamenetskii $\eta = N_0(z)\exp(E(T_0 - T))$, linear lower-mantle $N_0(z)$ from `2` | Zhang et al. (2010) / Bower et al. (2013) |
| `25` | EBA `20` variant, linear lower-mantle $N_0(z)$ from `5` to `100` | Flament EBA suite |
| `26` | EBA `20` variant, linear lower-mantle $N_0(z)$ from `5` to `80` | Flament EBA suite |
| `27` | EBA `20` variant, linear lower-mantle $N_0(z)$ from `10` to `100` between **1000–2867 km** | Rudolph et al. (2015) / Flament EBA |
| `28` | EBA `20` variant, linear lower-mantle $N_0(z)$ from `2` to `100` | Flament EBA suite |
| `112` | Exponential lower-mantle prefactor $N_0(z) = 0.5 \times 10^{4z/(r_o - r_i) - 0.5}$ | Hassan et al. |
| `113` | Depth-dependent lower-mantle activation term + LTBL viscosity clipping `[20, 100]` ($z > 0.43$) | Hassan et al. |
| `117` | Arrhenius (`rheol=7` variant) with $T\in[0,1]$ clamp and UTBL clip ($z < 0.06 \Rightarrow \eta \le 1$) | Hassan et al. |
| `118` | Same as `117` without UTBL viscosity clipping | Hassan et al. |

---

## 4. Hybrid Tracers & Extended Diagnostic Outputs

### Tracer Parameters
| Parameter | Default | Description |
|---|---|---|
| `hybrid_method` | `off` | Under ratio method (`buoy_type=1`), empty elements receive `C = 0.0` (ambient mantle) instead of retaining previous composition or triggering low-tracer warnings |
| `tracer_cushion` | `100000` | Extra per-rank tracer buffer capacity when loading `tracer_file` |

### Diagnostic Tokens for `output_optional` / `output_optional_ascii`
| Token | Output File Pattern | Description |
|---|---|---|
| `sten_temp` | `<datafile>.sten_temp.<rank>.<cycle>` | 3D nodal thermal slab assimilation stencil (`slab_sten`) |
| `sten_velo` | `<datafile>.sten_velo.<rank>.<cycle>` | 3D nodal internal velocity assimilation stencil (`slab_sten2`) |
| `divv` | `<datafile>.divv.<rank>.<cycle>` | Normalized nodal velocity divergence $\lVert \nabla \cdot \mathbf{u} \rVert_2 / \lVert \mathbf{u} \rVert_2$ |
| `tracer_dens` | `<datafile>.tracer_dens.<rank>.<cycle>` | Volume-normalized nodal tracer density (`E->trace.dens_node`) |
| `heating_visc_nd` | `<datafile>.heating_visc_nd.<rank>.<cycle>` | Nodal-projected viscous dissipation heating |
| `heating_adi_nd` | `<datafile>.heating_adi_nd.<rank>.<cycle>` | Nodal-projected adiabatic heating |
| `heating_latent_nd` | `<datafile>.heating_latent_nd.<rank>.<cycle>` | Nodal-projected phase-transition latent heating |
| `temp_sph` | `<datafile>.temp_sph.<rank>.<cycle>` | Radial-layer spherical harmonic expansion of temperature up to `output_ll_max` (Full sphere) |
| `comp_sph` | `<datafile>.comp_sph.<rank>.<cycle>` | Radial-layer spherical harmonic expansion of composition up to `output_ll_max` (Full sphere) |

---

## 5. Attribution & Required Citations

The features documented on this page are adapted from **EarthByte/CitcomS** (`https://github.com/EarthByte/citcoms`, GPL-2.0). Publications using these capabilities should cite:
- Bower, D. J., Gurnis, M., & Flament, N. (2015). *Assimilating lithosphere and slab history in 4-D Earth models*. Physics of the Earth and Planetary Interiors, 238, 8–22. [https://doi.org/10.1016/j.pepi.2014.10.013](https://doi.org/10.1016/j.pepi.2014.10.013)
- Flament, N., Gurnis, M., Williams, S., Seton, M., Skogseid, J., Heine, C., & Müller, R. D. (2014). *Topographic asymmetry of the South Atlantic from global models of mantle flow and lithospheric stretching*. Earth and Planetary Science Letters, 387, 107–119. [https://doi.org/10.1016/j.epsl.2013.11.017](https://doi.org/10.1016/j.epsl.2013.11.017)
- Zhang, N., Zhong, S., Leng, W., & Li, Z.-X. (2010). *A model for the evolution of the Earth's mantle structure since the Early Paleozoic*. JGR: Solid Earth, 115, B06401. [https://doi.org/10.1029/2009JB006896](https://doi.org/10.1029/2009JB006896)
- Rudolph, M. L., Lekić, V., & Lithgow-Bertelloni, C. (2015). *Viscosity jump in Earth's mid-mantle*. Science, 350(6266), 1349–1352. [https://doi.org/10.1126/science.aad1929](https://doi.org/10.1126/science.aad1929)
