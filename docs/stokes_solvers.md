# Modernized Stokes & Multigrid Solver Engine

← Back to [README](../README.md)

---

## 1. Zero-Drift Optimizations (Always On)

These structural improvements preserve 100% numerical equivalence with legacy CitcomS while reducing communication and memory overhead:
- **MPI Allreduce Fusion**: Fuses 5 consecutive Uzawa norm reductions into 1 (`global_stokes_norms`) and 2 upward multigrid reductions into 1 (`global_vdot2`).
- **Persistent Solver Workspaces**: Eliminates per-step heap allocations across `multi_grid`, `conj_grad`, and outer Uzawa solvers.
- **O(1) Static Stiffness Assembly LUT**: Precomputes `Node_map` indices for `construct_node_ks()`, replacing inner linear scans.
- **SIMD-Friendly Kernels**: Register-accumulated `gauss_seidel` and `n_assemble_del2_u`, single-pass outer-element `assemble_div_u`/`assemble_c_u`, and cached advection-diffusion spherical metrics.

---

## 2. Flexible Krylov Solvers & Modern Multigrid (Opt-in)

Legacy `.cfg` files without these keys preserve 100% original defaults (`uzawa=cg`, `mg_mode=fmg`, `mg_smoother=scalar_gs`, `coarse_solver=gs`, `p_precond=bpi`, `inner_rel_tol=0.0`).

For plate-driven global convection models with strong lateral viscosity contrasts ($10^4\text{–}10^6$), prescribed plate/slab velocity BCs, and augmented Lagrangian penalty (`aug_lagr=on`), the following drop-in configuration eliminates Uzawa stalling and delivers **2.4×–2.7×** Stokes solver speedup with **10×–20×** tighter divergence convergence:

```ini
Solver=multigrid
aug_lagr=on
aug_number=2.0e3
uzawa=fgmres                  # cg (default) | fcg | fgmres | bicg
fgmres_restart=30             # Krylov subspace dimension before restart
inner_rel_tol=0.05            # Relative inner momentum tolerance (ensures exact Schur operator)
mg_mode=vcycle                # fmg (default) | vcycle
mg_smoother=block3x3_gs       # scalar_gs (default) | block3x3_gs | block3x3_sgs | chebyshev
coarse_solver=cg              # gs (default) | cg
visc_smooth_method=4          # 4 = geometric-mean element coarsening (preserves weak plate boundaries)
```

> **Tip**: When using `mg_mode=vcycle`, pair it with `inner_rel_tol=0.05` (leaving `inner_mg_cycles=0`) rather than forcing `inner_mg_cycles=1`, so each inner momentum solve achieves adequate accuracy for the outer Schur-complement pressure solve.

---

## 3. Solver Parameter Reference Table

| Parameter | Options (Default First) | Description |
|---|---|---|
| `uzawa` | `cg`, `fcg`, `fgmres`, `bicg` | Outer Schur-complement Krylov solver (`fgmres` supports both incompressible and compressible ALA/BA flows) |
| `inner_rel_tol` | `0.0` (disabled) | Relative inner momentum tolerance $\text{tol} \cdot \lVert r_0 \rVert_2$; recommended `0.05` with `fgmres`/`fcg` |
| `fcg_mmax` | `10` | Truncation window $m_{\max}$ for Notay's Flexible CG (`uzawa=fcg`) |
| `fgmres_restart` | `30` | Restart dimension $m$ for W-Norm Flexible GMRES (`uzawa=fgmres`) |
| `residual_refresh_every` | `0` (disabled) | Optional periodic exact divergence recomputation interval for `uzawa=fcg` |
| `inner_mg_cycles` | `0` (tolerance-driven) | Optional fixed inner MG cycle cap (use only with `mg_mode=fmg` or `>= 2` with `vcycle`) |
| `p_precond` | `bpi`, `visc_mass` | Schur-complement pressure preconditioner |
| `mg_mode` | `fmg`, `vcycle` | Multigrid cycle mode (`vcycle` avoids redundant coarse traversals inside outer Krylov iterations) |
| `mg_smoother` | `scalar_gs`, `block3x3_gs`, `block3x3_sgs`, `chebyshev` | Multigrid velocity smoother (`block3x3_*` couples $(v_\theta, v_\phi, v_r)$ per node; `chebyshev` is MPI-rank-count invariant) |
| `coarse_solver` | `gs`, `cg` | Coarsest-level solver at `levmin` |
| `visc_smooth_method` | `0`, `1`, `2`, `3`, `4` | `4` adds geometric-mean element viscosity coarsening across multigrid levels |
