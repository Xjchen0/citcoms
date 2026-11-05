#!/usr/bin/env python3
"""
citcoms_h5_to_xdmf.py
Generate a ParaView-compatible XDMF (eXtensible Data Model and Format) XML index
for CitcomS HDF5 output files.

Usage:
    python citcoms_h5_to_xdmf.py --datafile stage1_preheat
    python citcoms_h5_to_xdmf.py --datafile stage1_preheat --dir ./output --scalet 1.0
"""

import os
import sys
import re
import argparse

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


def scan_h5_timesteps(data_dir, prefix):
    """
    Scan data_dir for files matching <prefix>.<cycles>.h5 and sort by cycles.
    """
    pattern = re.compile(rf'^{re.escape(prefix)}\.(\d+)\.h5$')
    steps = []

    if not os.path.isdir(data_dir):
        print(f"Error: Directory '{data_dir}' not found.", file=sys.stderr)
        return steps

    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if m:
            cycle = int(m.group(1))
            steps.append((cycle, fname))

    steps.sort(key=lambda x: x[0])
    return steps


def inspect_metadata(const_h5_path, steps, data_dir, user_scalet=None):
    """
    Inspect dimensions and time attributes from HDF5 files if h5py is available.
    """
    caps = 12
    nox = 65
    noy = 65
    noz = 65
    has_pressure = False
    has_stress = False
    timesteps_with_time = []
    scalet = user_scalet

    if HAS_H5PY and os.path.isfile(const_h5_path):
        try:
            with h5py.File(const_h5_path, 'r') as f:
                if 'coord' in f:
                    shape = f['coord'].shape
                    # shape is (caps, nox, noy, noz, 3)
                    if len(shape) == 5:
                        caps = shape[0]
                        nox = shape[1]
                        noy = shape[2]
                        noz = shape[3]
                if scalet is None and 'input' in f:
                    input_grp = f['input']
                    if 'scalet' in input_grp.attrs:
                        scalet = float(input_grp.attrs['scalet'])
        except Exception as e:
            print(f"Warning: Could not read metadata from '{const_h5_path}': {e}", file=sys.stderr)

    if scalet is None or scalet <= 0.0:
        scalet = 1.0

    for cycle, fname in steps:
        step_path = os.path.join(data_dir, fname)
        t_val = float(cycle)
        if HAS_H5PY:
            try:
                with h5py.File(step_path, 'r') as f:
                    if 'time' in f.attrs:
                        raw_t = float(f.attrs['time'])
                        t_val = raw_t * scalet if scalet > 0.0 else raw_t
                    elif 'time' in f:
                        raw_t = float(f['time'][...])
                        t_val = raw_t * scalet if scalet > 0.0 else raw_t
                    if 'pressure' in f:
                        has_pressure = True
                    if 'stress' in f:
                        has_stress = True
            except Exception:
                pass
        timesteps_with_time.append((cycle, t_val, fname))

    return caps, nox, noy, noz, has_pressure, has_stress, timesteps_with_time


def generate_xdmf(prefix, caps, nox, noy, noz, timesteps, has_pressure=False, has_stress=False):
    """
    Generate XDMF 2.0 XML string mapping CitcomS HDF5 datasets to 3DSMesh topologies.
    """
    const_h5 = f"{prefix}.h5"

    lines = [
        '<?xml version="1.0" ?>',
        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
        '<Xdmf Version="2.0">',
        '  <Domain>',
        '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    ]

    for cyc, t_val, fname in timesteps:
        step_h5 = fname
        if caps > 1:
            lines.append(f'      <Grid Name="Step_{cyc}" GridType="Collection" CollectionType="Spatial">')
        else:
            lines.append(f'      <Grid Name="Step_{cyc}" GridType="Uniform">')

        lines.append(f'        <Time Value="{t_val:.6g}"/>')

        for c in range(caps):
            if caps > 1:
                lines.append(f'        <Grid Name="Cap_{c}" GridType="Uniform">')

            lines.append(f'          <Topology TopologyType="3DSMesh" Dimensions="{nox} {noy} {noz}"/>')
            lines.append('          <Geometry GeometryType="XYZ">')
            lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz} 3" Type="HyperSlab">')
            lines.append('              <DataItem Dimensions="3 5" Format="XML">')
            lines.append(f'                {c} 0 0 0 0')
            lines.append('                1 1 1 1 1')
            lines.append(f'                1 {nox} {noy} {noz} 3')
            lines.append('              </DataItem>')
            lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz} 3" NumberType="Float" Precision="4" Format="HDF">')
            lines.append(f'                {const_h5}:/coord')
            lines.append('              </DataItem>')
            lines.append('            </DataItem>')
            lines.append('          </Geometry>')

            # Temperature
            lines.append('          <Attribute Name="temperature" AttributeType="Scalar" Center="Node">')
            lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz}" Type="HyperSlab">')
            lines.append('              <DataItem Dimensions="3 4" Format="XML">')
            lines.append(f'                {c} 0 0 0')
            lines.append('                1 1 1 1')
            lines.append(f'                1 {nox} {noy} {noz}')
            lines.append('              </DataItem>')
            lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz}" NumberType="Float" Precision="4" Format="HDF">')
            lines.append(f'                {step_h5}:/temperature')
            lines.append('              </DataItem>')
            lines.append('            </DataItem>')
            lines.append('          </Attribute>')

            # Velocity
            lines.append('          <Attribute Name="velocity" AttributeType="Vector" Center="Node">')
            lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz} 3" Type="HyperSlab">')
            lines.append('              <DataItem Dimensions="3 5" Format="XML">')
            lines.append(f'                {c} 0 0 0 0')
            lines.append('                1 1 1 1 1')
            lines.append(f'                1 {nox} {noy} {noz} 3')
            lines.append('              </DataItem>')
            lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz} 3" NumberType="Float" Precision="4" Format="HDF">')
            lines.append(f'                {step_h5}:/velocity')
            lines.append('              </DataItem>')
            lines.append('            </DataItem>')
            lines.append('          </Attribute>')

            # Viscosity
            lines.append('          <Attribute Name="viscosity" AttributeType="Scalar" Center="Node">')
            lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz}" Type="HyperSlab">')
            lines.append('              <DataItem Dimensions="3 4" Format="XML">')
            lines.append(f'                {c} 0 0 0')
            lines.append('                1 1 1 1')
            lines.append(f'                1 {nox} {noy} {noz}')
            lines.append('              </DataItem>')
            lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz}" NumberType="Float" Precision="4" Format="HDF">')
            lines.append(f'                {step_h5}:/viscosity')
            lines.append('              </DataItem>')
            lines.append('            </DataItem>')
            lines.append('          </Attribute>')

            # Pressure (optional)
            if has_pressure:
                lines.append('          <Attribute Name="pressure" AttributeType="Scalar" Center="Node">')
                lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz}" Type="HyperSlab">')
                lines.append('              <DataItem Dimensions="3 4" Format="XML">')
                lines.append(f'                {c} 0 0 0')
                lines.append('                1 1 1 1')
                lines.append(f'                1 {nox} {noy} {noz}')
                lines.append('              </DataItem>')
                lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz}" NumberType="Float" Precision="4" Format="HDF">')
                lines.append(f'                {step_h5}:/pressure')
                lines.append('              </DataItem>')
                lines.append('            </DataItem>')
                lines.append('          </Attribute>')

            # Stress (optional)
            if has_stress:
                lines.append('          <Attribute Name="stress" AttributeType="Tensor6" Center="Node">')
                lines.append(f'            <DataItem ItemType="HyperSlab" Dimensions="{nox} {noy} {noz} 6" Type="HyperSlab">')
                lines.append('              <DataItem Dimensions="3 5" Format="XML">')
                lines.append(f'                {c} 0 0 0 0')
                lines.append('                1 1 1 1 1')
                lines.append(f'                1 {nox} {noy} {noz} 6')
                lines.append('              </DataItem>')
                lines.append(f'              <DataItem Dimensions="{caps} {nox} {noy} {noz} 6" NumberType="Float" Precision="4" Format="HDF">')
                lines.append(f'                {step_h5}:/stress')
                lines.append('              </DataItem>')
                lines.append('            </DataItem>')
                lines.append('          </Attribute>')

            if caps > 1:
                lines.append('        </Grid>')

        lines.append('      </Grid>')

    lines.append('    </Grid>')
    lines.append('  </Domain>')
    lines.append('</Xdmf>')
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(
        description="Generate ParaView-compatible XDMF collection for CitcomS HDF5 outputs."
    )
    parser.add_argument("--datafile", "-f", required=True, help="Datafile prefix (e.g. stage1_preheat)")
    parser.add_argument("--dir", "-d", default=".", help="Directory containing HDF5 files (default: .)")
    parser.add_argument("--caps", type=int, default=None, help="Number of caps (default: auto or 12)")
    parser.add_argument("--nox", type=int, default=None, help="Number of nodes in X")
    parser.add_argument("--noy", type=int, default=None, help="Number of nodes in Y")
    parser.add_argument("--noz", type=int, default=None, help="Number of nodes in Z")
    parser.add_argument("--scalet", type=float, default=None, help="Time scaling factor to convert to Myr")
    parser.add_argument("--output", "-o", default=None, help="Output .xdmf filename (default: <dir>/<prefix>.xdmf)")

    args = parser.parse_args()

    data_dir = args.dir
    prefix = args.datafile
    const_h5 = os.path.join(data_dir, f"{prefix}.h5")

    steps = scan_h5_timesteps(data_dir, prefix)
    if not steps:
        print(f"No files matching '{prefix}.<cycle>.h5' found in '{data_dir}'.", file=sys.stderr)
        sys.exit(1)

    caps, nox, noy, noz, has_p, has_s, timesteps = inspect_metadata(const_h5, steps, data_dir, args.scalet)

    if args.caps is not None:
        caps = args.caps
    if args.nox is not None:
        nox = args.nox
    if args.noy is not None:
        noy = args.noy
    if args.noz is not None:
        noz = args.noz

    print(f"Generating XDMF for '{prefix}':")
    print(f"  Caps: {caps}")
    print(f"  Mesh: {nox} x {noy} x {noz}")
    print(f"  Timesteps: {len(timesteps)} found (from cycle {timesteps[0][0]} to {timesteps[-1][0]})")

    xdmf_content = generate_xdmf(prefix, caps, nox, noy, noz, timesteps, has_p, has_s)

    out_file = args.output if args.output else os.path.join(data_dir, f"{prefix}.xdmf")
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(xdmf_content)

    print(f"Successfully wrote XDMF index to '{out_file}'.")
    print(f"In ParaView, open '{out_file}' to visualize and animate in 3D.")


if __name__ == "__main__":
    main()
