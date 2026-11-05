#!/usr/bin/env python3
"""
Unit tests for ParaView PVD collection indexing and HDF5 XDMF temporal grid indexing.
"""

import os
import sys
import shutil
import tempfile
import subprocess
import xml.etree.ElementTree as ET

# Ensure bin is on path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
BIN_DIR = os.path.join(PROJECT_ROOT, "bin")
sys.path.insert(0, BIN_DIR)

from citcoms_h5_to_xdmf import generate_xdmf


def test_pvd_xml_generation_and_restart():
    """
    Test VTK PVD generation, time tracking, and restart pruning logic.
    """
    print("Testing VTK PVD generation and restart pruning...")
    tmpdir = tempfile.mkdtemp(prefix="test_pvd_")
    try:
        prefix = "test_model"
        hist_file = os.path.join(tmpdir, f"{prefix}.vtk_history.txt")
        pvd_file = os.path.join(tmpdir, f"{prefix}.pvd")

        # Step 1: Simulate write_pvd in Python mirroring the C implementation
        def simulate_write_pvd(cycles, time_val, caps=12):
            target_file = f"{prefix}.{cycles}.vtm" if caps == 12 else f"{prefix}.{cycles}.pvts"
            entries = []
            if os.path.isfile(hist_file):
                with open(hist_file, "r") as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            c = int(parts[0])
                            t = float(parts[1])
                            fn = parts[2]
                            if c < cycles:
                                entries.append((c, t, fn))
            entries.append((cycles, time_val, target_file))

            with open(hist_file, "w") as f:
                for c, t, fn in entries:
                    f.write(f"{c} {t:.8e} {fn}\n")

            with open(pvd_file, "w") as f:
                f.write('<?xml version="1.0"?>\n')
                f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
                f.write('  <Collection>\n')
                for c, t, fn in entries:
                    f.write(f'    <DataSet timestep="{t:.6g}" group="" part="0" file="{fn}"/>\n')
                f.write('  </Collection>\n')
                f.write('</VTKFile>\n')

        # Run 3 consecutive steps: 0, 1000, 2000
        simulate_write_pvd(0, 0.0)
        simulate_write_pvd(1000, 10.5)
        simulate_write_pvd(2000, 21.0)

        # Parse PVD
        tree = ET.parse(pvd_file)
        root = tree.getroot()
        assert root.tag == "VTKFile"
        assert root.attrib["type"] == "Collection"
        collection = root.find("Collection")
        datasets = collection.findall("DataSet")
        assert len(datasets) == 3
        assert datasets[0].attrib["timestep"] == "0"
        assert datasets[0].attrib["file"] == f"{prefix}.0.vtm"
        assert datasets[1].attrib["timestep"] == "10.5"
        assert datasets[1].attrib["file"] == f"{prefix}.1000.vtm"
        assert datasets[2].attrib["timestep"] == "21"
        assert datasets[2].attrib["file"] == f"{prefix}.2000.vtm"
        print("  [PASS] Initial 3-step PVD verified.")

        # Simulate restart rewind from cycle 1000 to cycle 2500
        # Restart occurs at 1000, then steps 1500 and 2500 are written
        simulate_write_pvd(1000, 10.5) # Re-writing 1000 prunes step 2000!
        simulate_write_pvd(1500, 15.75)
        simulate_write_pvd(2500, 26.25)

        tree = ET.parse(pvd_file)
        root = tree.getroot()
        datasets = root.find("Collection").findall("DataSet")
        assert len(datasets) == 4, f"Expected 4 datasets after rewind, got {len(datasets)}"
        steps = [d.attrib["file"] for d in datasets]
        assert steps == [
            f"{prefix}.0.vtm",
            f"{prefix}.1000.vtm",
            f"{prefix}.1500.vtm",
            f"{prefix}.2500.vtm"
        ], f"Unexpected steps after rewind: {steps}"
        print("  [PASS] Restart rewind and auto-pruning verified.")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_xdmf_generation_and_structure():
    """
    Test XDMF generation for both 12-cap full sphere and 1-cap regional models.
    """
    print("Testing XDMF generation and XML schema...")
    prefix = "global_preheat"
    caps = 12
    nox, noy, noz = 65, 65, 65
    timesteps = [
        (0, 0.0, f"{prefix}.0.h5"),
        (1000, 10.5, f"{prefix}.1000.h5"),
        (2000, 21.0, f"{prefix}.2000.h5")
    ]

    xml_str = generate_xdmf(prefix, caps, nox, noy, noz, timesteps, has_pressure=True, has_stress=True)

    # Validate XML syntax
    root = ET.fromstring(xml_str)
    assert root.tag == "Xdmf"
    assert root.attrib["Version"] == "2.0"

    domain = root.find("Domain")
    time_series = domain.find("Grid")
    assert time_series.attrib["Name"] == "TimeSeries"
    assert time_series.attrib["CollectionType"] == "Temporal"

    time_steps = time_series.findall("Grid")
    assert len(time_steps) == 3

    for idx, step_grid in enumerate(time_steps):
        assert step_grid.attrib["CollectionType"] == "Spatial"
        time_elem = step_grid.find("Time")
        expected_t = f"{timesteps[idx][1]:.6g}"
        assert time_elem.attrib["Value"] == expected_t

        cap_grids = step_grid.findall("Grid")
        assert len(cap_grids) == 12
        for c, cap in enumerate(cap_grids):
            assert cap.attrib["Name"] == f"Cap_{c}"
            topo = cap.find("Topology")
            assert topo.attrib["TopologyType"] == "3DSMesh"
            assert topo.attrib["Dimensions"] == f"{nox} {noy} {noz}"

            geom = cap.find("Geometry")
            assert geom.attrib["GeometryType"] == "XYZ"
            data_items = geom.findall("DataItem")
            assert len(data_items) >= 1
            hs = data_items[0]
            assert hs.attrib["ItemType"] == "HyperSlab"
            assert hs.attrib["Dimensions"] == f"{nox} {noy} {noz} 3"

            # Check attributes: temperature, velocity, viscosity, pressure, stress
            attrs = {a.attrib["Name"]: a for a in cap.findall("Attribute")}
            assert "temperature" in attrs
            assert attrs["temperature"].attrib["AttributeType"] == "Scalar"
            assert "velocity" in attrs
            assert attrs["velocity"].attrib["AttributeType"] == "Vector"
            assert "viscosity" in attrs
            assert attrs["viscosity"].attrib["AttributeType"] == "Scalar"
            assert "pressure" in attrs
            assert "stress" in attrs
            assert attrs["stress"].attrib["AttributeType"] == "Tensor6"

    print("  [PASS] 12-cap XDMF structure and attributes fully verified.")

    # Also test regional 1-cap
    xml_reg = generate_xdmf("regional", 1, 33, 33, 33, [(0, 0.0, "regional.0.h5")])
    root_reg = ET.fromstring(xml_reg)
    step0 = root_reg.find("Domain").find("Grid").find("Grid")
    assert step0.attrib["GridType"] == "Uniform"
    topo = step0.find("Topology")
    assert topo.attrib["TopologyType"] == "3DSMesh"
    assert topo.attrib["Dimensions"] == "33 33 33"
    print("  [PASS] Regional 1-cap XDMF verified.")


def test_offline_cli_tool():
    """
    Test bin/citcoms_h5_to_xdmf.py CLI interface.
    """
    print("Testing offline citcoms_h5_to_xdmf.py CLI...")
    tmpdir = tempfile.mkdtemp(prefix="test_cli_")
    try:
        prefix = "stage1"
        # Create dummy h5 files
        with open(os.path.join(tmpdir, f"{prefix}.h5"), "w") as f:
            f.write("dummy")
        with open(os.path.join(tmpdir, f"{prefix}.0.h5"), "w") as f:
            f.write("dummy")
        with open(os.path.join(tmpdir, f"{prefix}.1000.h5"), "w") as f:
            f.write("dummy")

        cmd = [
            sys.executable,
            os.path.join(BIN_DIR, "citcoms_h5_to_xdmf.py"),
            "--datafile", prefix,
            "--dir", tmpdir,
            "--caps", "12",
            "--nox", "65",
            "--noy", "65",
            "--noz", "65",
            "--scalet", "10.0"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        assert res.returncode == 0, f"CLI failed: {res.stderr}"

        xdmf_out = os.path.join(tmpdir, f"{prefix}.xdmf")
        assert os.path.isfile(xdmf_out), "Output xdmf file was not created"

        tree = ET.parse(xdmf_out)
        root = tree.getroot()
        assert root.tag == "Xdmf"
        print("  [PASS] Offline CLI tool created valid XDMF.")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    test_pvd_xml_generation_and_restart()
    test_xdmf_generation_and_structure()
    test_offline_cli_tool()
    print("\nALL PVD AND XDMF TESTS PASSED SUCCESSFULLY!")
