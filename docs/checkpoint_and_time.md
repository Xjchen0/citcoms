# Fault-Tolerant Checkpointing & Geological Age Scheduling

← Back to [README](../README.md)

---

## 1. Crash-Safe Checkpointing & Auto-Resume

Designed for unattended HPC batch-queue job chaining under strict walltime limits (`SIGKILL`):

```ini
restart=auto                       # Auto-resume from latest verified checkpoint, or start at cycle 0
checkpoints_to_keep=3              # Sliding window: retain last 3 checkpoints (0 = keep all)
checkpoint_interval_myr=5.0        # Physical geological age checkpoint spacing (Myr)
checkpoint_walltime=7200           # Optional wall-clock checkpoint interval (seconds)
checkpoint_on_termination=on       # Always save a final checkpoint before clean exit
```

### Key Mechanisms
- **Zero-Touch Auto-Resume (`restart=auto`)**: Reads `.last_good_checkpoint` to resume seamlessly from the latest committed cycle; starts cleanly from cycle 0 if none exists — no `.cfg` edits needed between queued jobs.
- **Two-Phase Commit**: Synchronizes all MPI ranks via `parallel_process_sync` before Rank 0 atomically updates `.last_good_checkpoint`, eliminating half-written/corrupted restarts.
- **Sliding-Window Pruning (`checkpoints_to_keep`)**: Automatically deletes checkpoints older than the last $N$ verified commits across all ranks to prevent inode and disk quota exhaustion.

---

## 2. Two-Stage Binary Continuation (ASCII-Independent)

Stage 2 models can inherit temperature and tracer fields directly from Stage 1 binary `.chkpt` files without exporting massive uncompressed ASCII velocity/tracer dumps:

```ini
tic_method=-1
tracer=on
tracer_ic_method=2
datadir_old="./stage1_output"
datafile_old="model_stage1"
solution_cycles_init=5000
```

*(Falls back automatically to legacy ASCII `.velo`/`.tracer` files if binary `.chkpt` files are absent.)*

---

## 3. Geological Age Triggers, Exact Snapping & Endpoints

Schedule outputs, checkpoints, and simulation termination directly in physical **Myr** (automatically converted via internal `scalet = (R^2 / kappa) / (1e6 * 365.25 * 24 * 3600)`):

```ini
start_age=200.0                    # Initial geological age (Ma)
storage_interval_myr=1.0           # Output every 1.0 Myr (EarthByte alias: storage_spacing_Myr=1)
checkpoint_interval_myr=5.0        # Checkpoint every 5.0 Myr
exact_time_output=on               # Clamp CFL dt to land exactly on target Myr milestones
end_age_myr=100.0                  # Stop cleanly after 100.0 Myr elapsed (or use exit_at_present=on)
exit_at_present=off                # Set on to terminate forward models cleanly at 0.0 Ma
```

### Parameter Summary

| Parameter | Default | Description |
|---|---|---|
| `restart` | `off` | `off` / `on` / `auto` (zero-touch auto-resume from `.last_good_checkpoint`) |
| `checkpoints_to_keep` | `3` | Number of recent checkpoints retained (`0` = unlimited) |
| `checkpoint_on_termination` | `on` | Force saving a valid checkpoint on clean completion / walltime / step limit |
| `checkpoint_walltime` | `0.0` | Real wall-clock seconds between automatic checkpoints (`0` = disabled) |
| `checkpoint_interval_myr` | `0.0` | Checkpoint spacing in physical Myr (`0` = disabled) |
| `checkpoint_interval_time` | `0.0` | Checkpoint spacing in dimensionless time (`0` = disabled) |
| `storage_interval_myr` | `0.0` | Field output spacing in physical Myr (disables default step counter if `storage_spacing` is omitted) |
| `storage_spacing_Myr` | `0` | EarthByte compatibility alias mapped to `storage_interval_myr` (inherits exact snapping) |
| `storage_interval_time` | `0.0` | Field output spacing in dimensionless time |
| `exact_time_output` | `on` | Dynamically snaps CFL $\Delta t$ so outputs and termination land precisely on target Myr boundaries |
| `end_age_myr` / `end_time_myr` | `0.0` | Terminate simulation after given elapsed physical Myr |
| `end_time` | `0.0` | Terminate simulation at given dimensionless elapsed time |
| `exit_at_present` | `off` | Terminate forward convection (`timedir >= 0`) cleanly when geological age reaches `0.0 Ma` |
