#!/usr/bin/env python3
"""Comprehensive end-to-end verification of EarthByte features integrated into citcoms_cpp."""
import math
import os
import re
import shutil
import subprocess
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD_DIR = os.environ.get("CITCOMS_BUILD_DIR", os.path.join(REPO_ROOT, "build_linux"))
REG_BIN = os.path.join(BUILD_DIR, "CitcomSRegional")
FULL_BIN = os.path.join(BUILD_DIR, "CitcomSFull")


def run_cmd(cmd, cwd=None, timeout=180):
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    if res.returncode != 0:
        print("FAILED:", " ".join(cmd))
        print("STDOUT:", res.stdout[-2000:])
        print("STDERR:", res.stderr[-2000:])
        raise RuntimeError(f"Command failed with rc={res.returncode}")
    return res


def make_cfg(params):
    return "\n".join(f"{k}={v}" for k, v in params.items()) + "\n"


def base_regional_params(datadir):
    return {
        "datadir": datadir,
        "datafile": "test",
        "datadir_old": datadir,
        "datafile_old": "test",
        "output_format": "ascii",
        "output_optional": "horiz_avg,velo,visc",
        "nproc_surf": 1,
        "nprocx": 1,
        "nprocy": 1,
        "nprocz": 1,
        "nodex": 9,
        "nodey": 9,
        "nodez": 9,
        "mgunitx": 2,
        "mgunity": 2,
        "mgunitz": 2,
        "levels": 3,
        "theta_min": 1.0,
        "theta_max": 1.5,
        "fi_min": 0.0,
        "fi_max": 0.5,
        "radius_inner": 0.55,
        "radius_outer": 1.0,
        "minstep": 1,
        "maxstep": 3,
        "maxtotstep": 100,
        "storage_spacing": 1,
        "checkpointFrequency": 100,
        "cpu_limits_in_seconds": 3600,
        "rayleigh": "1e6",
        "dissipation_number": 0.0,
        "gruneisen": 0.0,
        "Solver": "multigrid",
        "uzawa": "fgmres",
        "mg_mode": "vcycle",
        "mg_smoother": "block3x3_gs",
        "accuracy": "1e-4",
        "piterations": 50,
        "vhighstep": 3,
        "vlowstep": 10,
        "topvbc": 0,
        "botvbc": 0,
        "toptbc": 1,
        "toptbcval": 0.0,
        "bottbc": 1,
        "bottbcval": 1.0,
        "tic_method": 0,
        "num_perturbations": 1,
        "perturbl": 2,
        "perturbm": 1,
        "perturblayer": 5,
        "perturbmag": 0.05,
        "mantle_temp": 0.5,
        "ADVECTION": "on",
        "monitor_max_T": "off",
        "finetunedt": 0.5,
        "Viscosity": "system",
        "VISC_UPDATE": "on",
        "num_mat": 4,
        "TDEPV": "on",
        "rheol": 3,
        "visc0": "1,1,1,10",
        "viscE": "2.0,2.0,2.0,2.0",
        "viscT": "0.5,0.5,0.5,0.5",
        "viscZ": "0.1,0.1,0.1,0.1",
        "VMIN": "on",
        "visc_min": 0.01,
        "VMAX": "on",
        "visc_max": 1000.0,
        "checkpoint_on_termination": "off",
    }


def parse_total_heating(stderr_text):
    matches = re.findall(r"Total_heating\(visc, adi\):\s+([-\deE.+]+)\s+([-\deE.+]+)", stderr_text)
    assert matches, f"No Total_heating found in stderr:\n{stderr_text[-1000:]}"
    return float(matches[-1][0]), float(matches[-1][1])


def test_eba_dis_scaling(root):
    print("[1/5] Testing EBA Di scaling & 100% legacy refstate_file compatibility...")
    ref7 = os.path.join(root, "ref7.dat")
    with open(ref7, "w") as f:
        for _ in range(9):
            f.write("1.0 1.0 1.0 1.0 0.0 0.0 0.0\n")
    ref5 = os.path.join(root, "ref5.dat")
    with open(ref5, "w") as f:
        for _ in range(9):
            f.write("1.0 1.0 1.0 1.0 2.5\n")

    def run_eba(tag, ref_file, eba_flag):
        dd = os.path.join(root, f"eba_{tag}")
        os.makedirs(dd, exist_ok=True)
        cfg = os.path.join(root, f"eba_{tag}.cfg")
        p = base_regional_params(dd)
        p.update({
            "reference_state": 0,
            "refstate_file": ref_file,
            "eba_dis_scaling": eba_flag,
            "dissipation_number": 0.5,
            "output_optional": "horiz_avg,heating,heating_visc_nd,heating_adi_nd,heating_latent_nd",
        })
        with open(cfg, "w") as f:
            f.write(make_cfg(p))
        res = run_cmd(["mpirun", "-np", "1", REG_BIN, cfg])
        for kind in ["heating_visc_nd", "heating_adi_nd", "heating_latent_nd"]:
            assert os.path.exists(os.path.join(dd, f"test.{kind}.0.1")), f"Missing {kind}"
        return parse_total_heating(res.stderr)

    v_legacy, a_legacy = run_eba("legacy7_off", ref7, "off")
    v_zero, a_zero = run_eba("legacy7_on", ref7, "on")
    v_scaled, a_scaled = run_eba("eb5_on", ref5, "on")

    print(f"    legacy7_off: visc={v_legacy:.6e}, adi={a_legacy:.6e}")
    print(f"    legacy7_on : visc={v_zero:.6e}, adi={a_zero:.6e}")
    print(f"    eb5_on(2.5): visc={v_scaled:.6e}, adi={a_scaled:.6e}")
    assert abs(v_legacy) > 1e-12, "Legacy col5=0 with eba_dis_scaling=off should have non-zero visc heating!"
    assert abs(v_zero) < 1e-15 and abs(a_zero) < 1e-15, "col5=0 with eba_dis_scaling=on should zero heating!"
    ratio = v_scaled / v_legacy
    assert 2.2 < ratio < 2.8, f"Expected ~2.5x viscous heating scaling, got {ratio:.3f}"
    print("    PASS: EBA Di scaling & legacy compatibility verified!")


def test_viscosity_laws(root):
    print("[2/5] Testing EarthByte viscosity laws (rheol=20..28, 112, 113, 117, 118)...")
    for r in [20, 21, 22, 23, 24, 25, 26, 27, 28, 112, 113, 117, 118]:
        dd = os.path.join(root, f"visc_{r}")
        os.makedirs(dd, exist_ok=True)
        cfg = os.path.join(root, f"visc_{r}.cfg")
        p = base_regional_params(dd)
        p.update({"maxstep": 1, "rheol": r})
        with open(cfg, "w") as f:
            f.write(make_cfg(p))
        run_cmd(["mpirun", "-np", "1", REG_BIN, cfg])
        vfile = os.path.join(dd, "test.visc.0.1")
        assert os.path.exists(vfile)
        with open(vfile) as f:
            vals = [float(line.strip()) for line in f.readlines()[1:] if line.strip()]
        assert all(math.isfinite(x) and x > 0 for x in vals), f"Invalid viscosity for rheol={r}"
    print("    PASS: All 13 EarthByte viscosity laws ran cleanly with finite positive viscosity!")


def test_4d_assimilation_and_restart(root):
    print("[3/5] Testing 4D lith/slab/internal-velocity assimilation + storage_spacing_Myr + exit_at_present + restart...")
    assim_dir = os.path.join(root, "assim_inputs")
    os.makedirs(assim_dir, exist_ok=True)
    nox, noy, noz = 9, 9, 9
    # Create files for ages 0..4 Ma (NO age 5 file -> verifies exact integer age single-file read!)
    for age in range(0, 5):
        with open(os.path.join(assim_dir, f"bvel{age}"), "w") as f:
            for _ in range(nox * noy):
                f.write("2.0 -1.0\n")
        with open(os.path.join(assim_dir, f"age{age}"), "w") as f:
            for j in range(1, noy + 1):
                for i in range(1, nox + 1):
                    if i <= 4:
                        f.write("-999.0\n")  # deforming / masked region
                    else:
                        f.write(f"{40.0 + age:.1f}\n")
        with open(os.path.join(assim_dir, f"slab{age}"), "w") as f:
            for j in range(1, noy + 1):
                for i in range(1, nox + 1):
                    for k in range(1, noz + 1):
                        if i == 5 and j == 5 and k == 6:
                            f.write("0.15 1.0\n")
                        else:
                            f.write("0.5 0.0\n")
        with open(os.path.join(assim_dir, f"ivel{age}"), "w") as f:
            for j in range(1, noy + 1):
                for i in range(1, nox + 1):
                    for k in range(1, noz + 1):
                        if k in (1, noz) or i in (1, nox) or j in (1, noy):
                            f.write("0.0 0.0 0.0 2\n")
                        elif i == 5 and j == 5 and k == 6:
                            f.write("1.5 -0.5 -2.5 1\n")
                        else:
                            f.write("0.0 0.0 0.0 0\n")

    dd = os.path.join(root, "assim_out")
    os.makedirs(dd, exist_ok=True)
    cfg = os.path.join(root, "assim.cfg")
    p = base_regional_params(dd)
    del p["storage_spacing"]  # let storage_spacing_Myr drive output spacing!
    p.update({
        "maxstep": 500,
        "maxtotstep": 500,
        "storage_spacing_Myr": 2,
        "checkpointFrequency": 2,
        "exit_at_present": "on",
        "exact_time_output": "on",
        "topvbc": 1,
        "botvbc": 0,
        "file_vbcs": 1,
        "vel_bound_file": f"{assim_dir}/bvel",
        "start_age": 4.0,
        "lith_age": 1,
        "lith_age_time": 1,
        "lith_age_file": f"{assim_dir}/age",
        "lith_age_depth_function": 1,
        "lith_age_depth": 1.0,
        "lith_age_exponent": 0.5,
        "lith_age_min": 1.0,
        "lith_age_stencil_value": -100.0,
        "slab_assim": 1,
        "slab_assim_file": f"{assim_dir}/slab",
        "internal_vbcs_file": 1,
        "velocity_internal_file": f"{assim_dir}/ivel",
        "output_optional": "velo,visc,horiz_avg,sten_temp,sten_velo,divv",
        "rayleigh": "1e7",
        "finetunedt": 0.5,
    })
    with open(cfg, "w") as f:
        f.write(make_cfg(p))

    res = run_cmd(["mpirun", "-np", "1", REG_BIN, cfg])
    assert "Reached present day" in res.stderr, f"exit_at_present did not trigger!\n{res.stderr[-1500:]}"

    velo_files = sorted(
        [f for f in os.listdir(dd) if f.startswith("test.velo.0.")],
        key=lambda s: int(s.split(".")[-1]),
    )
    assert len(velo_files) >= 3, f"Expected multiple storage_spacing_Myr outputs (0, 2, 4 Ma), got {velo_files}"
    last_cycle = int(velo_files[-1].split(".")[-1])
    for kind in ["sten_temp", "sten_velo", "divv"]:
        fp = os.path.join(dd, f"test.{kind}.0.{last_cycle}")
        assert os.path.exists(fp), f"Missing diagnostic output {fp}"

    node_idx = 6 + (5 - 1) * noz + (5 - 1) * nox * noz
    with open(os.path.join(dd, f"test.velo.0.{last_cycle}")) as f:
        lines = f.readlines()[2:]
    v1, v2, v3, Tval = [float(x) for x in lines[node_idx - 1].split()]
    scalev = (6371e3 / 1e-6) / (100.0 * 365.25 * 24.0 * 3600.0)
    assert abs(Tval - 0.15) < 1e-4, f"Expected slab T=0.15 at assimilated node, got {Tval}"
    assert abs(v1 - 1.5 * scalev) < 1e-2 * scalev, f"Expected v1={1.5*scalev}, got {v1}"
    assert abs(v3 - (-2.5 * scalev)) < 1e-2 * scalev, f"Expected v3={-2.5*scalev}, got {v3}"

    with open(os.path.join(dd, f"test.sten_temp.0.{last_cycle}")) as f:
        st_val = float(f.readlines()[2 + node_idx - 1].strip())
    with open(os.path.join(dd, f"test.sten_velo.0.{last_cycle}")) as f:
        sv_val = int(f.readlines()[2 + node_idx - 1].strip())
    assert abs(st_val - 1.0) < 1e-5 and sv_val == 1, f"Stencil mismatch: {st_val}, {sv_val}"

    chkpts = sorted(
        [int(f.split(".")[-1]) for f in os.listdir(dd) if f.startswith("test.chkpt.0.")]
    )
    assert chkpts, "No checkpoints created"
    restart_cycle = chkpts[0]
    cfg_r = os.path.join(root, "assim_restart.cfg")
    p_r = dict(p)
    p_r.update({
        "restart": "on",
        "solution_cycles_init": restart_cycle,
        "maxstep": 2,
    })
    with open(cfg_r, "w") as f:
        f.write(make_cfg(p_r))
    run_cmd(["mpirun", "-np", "1", REG_BIN, cfg_r])
    print(f"    PASS: 4D assimilation (T_slab={Tval:.4f}, v3_match=OK, outputs={len(velo_files)}, exit@0Ma, restart@cycle={restart_cycle}) verified!")


def test_buoyancy_cutoff_and_hybrid_tracers(root):
    print("[4/5] Testing remove_buoyancy_above_znode & hybrid_method + tracer_dens...")
    dd_full = os.path.join(root, "buoy_full")
    dd_cut = os.path.join(root, "buoy_cut")
    os.makedirs(dd_full, exist_ok=True)
    os.makedirs(dd_cut, exist_ok=True)
    for tag, dd, zcut in [("full", dd_full, 10000), ("cut", dd_cut, 1)]:
        cfg = os.path.join(root, f"buoy_{tag}.cfg")
        p = base_regional_params(dd)
        p.update({"maxstep": 1, "remove_buoyancy_above_znode": zcut})
        with open(cfg, "w") as f:
            f.write(make_cfg(p))
        run_cmd(["mpirun", "-np", "1", REG_BIN, cfg])

    def max_vel(dd):
        with open(os.path.join(dd, "test.velo.0.1")) as f:
            lines = f.readlines()[2:]
        return max(
            math.sqrt(sum(float(x) ** 2 for x in l.split()[:3]))
            for l in lines
            if l.strip()
        )

    vf, vc = max_vel(dd_full), max_vel(dd_cut)
    print(f"    max|V| full={vf:.4e}, cutoff(znode=1)={vc:.4e}")
    assert vf > 1e-2 and vc < 1e-6 * vf, f"Expected buoyancy cutoff to suppress flow! vf={vf}, vc={vc}"

    # Hybrid method with sparse tracer file (20 tracers in 512 elements -> >95% empty!)
    tfile = os.path.join(root, "sparse_tracers.dat")
    with open(tfile, "w") as f:
        f.write("20 4\n")
        for k in range(20):
            f.write(f"{1.15 + 0.01*k:.6f} {0.15 + 0.01*k:.6f} {0.65 + 0.005*k:.6f} 1.0\n")
    dd_hyb = os.path.join(root, "hybrid_out")
    os.makedirs(dd_hyb, exist_ok=True)
    cfg_h = os.path.join(root, "hybrid.cfg")
    p_h = base_regional_params(dd_hyb)
    p_h.update({
        "maxstep": 2,
        "tracer": "on",
        "tracer_ic_method": 1,
        "tracer_file": tfile,
        "tracer_flavors": 2,
        "ic_method_for_flavors": 0,
        "z_interface": 0.7,
        "chemical_buoyancy": "on",
        "buoy_type": 1,
        "hybrid_method": "on",
        "buoyancy_ratio": 0.5,
        "itracer_warnings": "on",
        "output_optional": "horiz_avg,comp_nd,tracer_dens",
    })
    with open(cfg_h, "w") as f:
        f.write(make_cfg(p_h))
    run_cmd(["mpirun", "-np", "1", REG_BIN, cfg_h])
    assert os.path.exists(os.path.join(dd_hyb, "test.tracer_dens.0.1"))
    with open(os.path.join(dd_hyb, "test.comp_nd.0.1")) as f:
        clines = f.readlines()[2:]
    cvals = [float(l.strip()) for l in clines if l.strip()]
    assert any(v == 0.0 for v in cvals) and any(v > 0.0 for v in cvals), "Expected both 0 (empty) and >0 (tracer) comp_nd"
    print("    PASS: Buoyancy cutoff, hybrid_method (sparse tracers), and tracer_dens verified!")


def test_full_sphere_sph_outputs(root):
    print("[5/5] Testing Full spherical temp_sph & comp_sph outputs (12 MPI ranks)...")
    dd = os.path.join(root, "full_sph")
    os.makedirs(dd, exist_ok=True)
    cfg = os.path.join(root, "full_sph.cfg")
    p = base_regional_params(dd)
    for k in ("theta_min", "theta_max", "fi_min", "fi_max"):
        del p[k]
    p.update({
        "datafile": "fulltest",
        "datafile_old": "fulltest",
        "nproc_surf": 12,
        "maxstep": 1,
        "rayleigh": "1e5",
        "accuracy": "1e-3",
        "piterations": 30,
        "rheol": 22,
        "tracer": "on",
        "tracer_ic_method": 0,
        "tracers_per_element": 5,
        "tracer_flavors": 2,
        "ic_method_for_flavors": 0,
        "z_interface": 0.7,
        "chemical_buoyancy": "on",
        "buoy_type": 1,
        "buoyancy_ratio": 0.3,
        "output_ll_max": 4,
        "output_optional": "horiz_avg,temp_sph,comp_sph,divv,tracer_dens",
    })
    with open(cfg, "w") as f:
        f.write(make_cfg(p))
    run_cmd(["mpirun", "--oversubscribe", "-np", "12", FULL_BIN, cfg])
    for kind in ["temp_sph", "comp_sph", "divv", "tracer_dens"]:
        fp = os.path.join(dd, f"fulltest.{kind}.0.1")
        assert os.path.exists(fp) and os.path.getsize(fp) > 0, f"Missing or empty {fp}"
    print("    PASS: Full spherical temp_sph, comp_sph, divv, tracer_dens & rheol=22 verified on 12 ranks!")


def test_parser_unquoted_and_spacing(root):
    print("[6/6] Testing parser: unquoted non-path params, spaces around '=', inline comments, CRLF...")
    dd = os.path.join(root, "parse_test")
    os.makedirs(dd, exist_ok=True)
    cfg = os.path.join(root, "parse_test.cfg")
    # Mix quoted paths, unquoted non-path params, spaces around '=', spaces after commas, comments, CRLF
    lines = [
        f'datadir = "{dd}"',
        'datafile="parsetest"',
        "output_format = ascii, vtk   # composite without quotes",
        "output_optional_ascii = horiz_avg, divv  # comma list with spaces",
        "output_optional_vtk = surf",
        "vtk_output_mode = cap",
        "vtk_format = ascii",
        "nproc_surf=1 nprocx=1 nprocy=1 nprocz=1",
        "nodex = 9",
        "nodey = 9",
        "nodez = 9",
        "mgunitx=2",
        "mgunity=2",
        "mgunitz=2",
        "levels=3",
        "theta_min=1.0",
        "theta_max=1.5",
        "fi_min=0.0",
        "fi_max=0.5",
        "radius_inner=0.55",
        "radius_outer=1.0",
        "minstep=1",
        "maxstep=1",
        "maxtotstep=10",
        "storage_spacing=1",
        "checkpointFrequency=100",
        "cpu_limits_in_seconds=3600",
        "rayleigh = 1e6",
        "Solver = multigrid",
        "uzawa = fgmres # inline comment",
        "mg_mode = vcycle",
        "mg_smoother = block3x3_gs",
        "accuracy = 1e-4",
        "piterations = 50",
        "topvbc=0",
        "botvbc=0",
        "toptbc=1",
        "toptbcval=0.0",
        "bottbc=1",
        "bottbcval=1.0",
        "tic_method=0",
        "num_perturbations=1",
        "perturbl=2",
        "perturbm=1",
        "perturblayer=5",
        "perturbmag=0.05",
        "mantle_temp=0.5",
        "Viscosity = system",
        "VISC_UPDATE = true",
        "num_mat = 4",
        "TDEPV = yes",
        "rheol = 20",
        "visc0 = 1, 1, 1, 10",
        "viscE = 2.0, 2.0, 2.0, 2.0",
        "viscT = 0.5, 0.5, 0.5, 0.5",
        "checkpoint_on_termination = false",
    ]
    with open(cfg, "wb") as f:
        f.write(("\r\n".join(lines) + "\r\n").encode("utf-8"))
    run_cmd(["mpirun", "-np", "1", REG_BIN, cfg])
    assert os.path.exists(os.path.join(dd, "parsetest.horiz_avg.0.1"))
    assert os.path.exists(os.path.join(dd, "parsetest.divv.0.1"))
    assert os.path.exists(os.path.join(dd, "parsetest.pvd"))
    assert not os.path.exists(os.path.join(dd, "parsetest.velo.0.1")), "Composite ascii,vtk should suppress 3D ASCII velo!"
    print("    PASS: Parser handles unquoted non-path params, spaces around '=', comments, and CRLF cleanly!")


def main():
    root = tempfile.mkdtemp(prefix="verify_earthbyte_")
    try:
        test_eba_dis_scaling(root)
        test_viscosity_laws(root)
        test_4d_assimilation_and_restart(root)
        test_buoyancy_cutoff_and_hybrid_tracers(root)
        test_full_sphere_sph_outputs(root)
        test_parser_unquoted_and_spacing(root)
        print("\n=== ALL 6 EARTHBYTE & PARSER INTEGRATION TEST SUITES PASSED! ===")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
