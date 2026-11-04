#!/usr/bin/env python3
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

'''
Cut-paste and combine Citcom surf data (Accelerated via Modern C++ citcoms_combine)

usage: combinesurf.py modelname timestep nodex nodey nodez n_surf_proc nprocx nprocy nprocz
'''

import os
import sys
import subprocess
import shutil

def find_cpp_combiner():
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "build", "citcoms_combine"),
        os.path.join(os.path.dirname(__file__), "..", "build_clang", "citcoms_combine"),
        os.path.join(os.path.dirname(__file__), "..", "bin", "citcoms_combine"),
        shutil.which("citcoms_combine") or ""
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None

class CombineSurf(object):

    def __init__(self, grid):
        self.saved = [''] * (grid['nox'] * grid['noy'])

    def readData(self, crdfilename, surffilename, grid, cap):
        nprocx = int(cap['nprocx'])
        nprocy = int(cap['nprocy'])
        nprocz = int(cap['nprocz'])

        nox = grid['nox']
        noy = grid['noy']
        noz = grid['noz']

        mynox = 1 + (nox - 1) // nprocx
        mynoy = 1 + (noy - 1) // nprocy
        mynoz = 1 + (noz - 1) // nprocz
        snodes = mynox * mynoy

        with open(crdfilename, 'r') as fp1, open(surffilename, 'r') as fp2:
            fp1.readline()
            fp2.readline()
            surf = fp2.readlines()
            crd = fp1.readlines()

        data = [''] * snodes
        for i in range(len(data)):
            parts = crd[(i + 1) * mynoz - 1].split()
            x, y = parts[0], parts[1]
            data[i] = f'{x} {y} {surf[i]}'

        return data

    def join(self, data, me, grid, cap):
        nprocx = cap['nprocx']
        nprocy = cap['nprocy']
        nprocz = cap['nprocz']

        mylocz = me % nprocz
        mylocx = ((me - mylocz) // nprocz) % nprocx
        mylocy = (((me - mylocz) // nprocz - mylocx) // nprocx) % nprocy

        nox = int(grid['nox'])
        noy = int(grid['noy'])

        mynox = 1 + (nox - 1) // nprocx
        mynoy = 1 + (noy - 1) // nprocy

        mynxs = (mynox - 1) * mylocx
        mynys = (mynoy - 1) * mylocy

        n = 0
        for i in range(mynys, mynys + mynoy):
            for j in range(mynxs, mynxs + mynox):
                m = j + i * nox
                self.saved[m] = data[n]
                n += 1

    def write(self, filename, grid, data):
        with open(filename, 'w') as fp:
            fp.write(f"{grid['nox']} x {grid['noy']}\n")
            fp.writelines(data)


if __name__ == '__main__':
    if len(sys.argv) != 10:
        print(__doc__)
        sys.exit(1)

    prefix = sys.argv[1]
    step = int(sys.argv[2])

    grid = {
        'nox': int(sys.argv[3]),
        'noy': int(sys.argv[4]),
        'noz': int(sys.argv[5])
    }

    nprocxy = int(sys.argv[6])
    cap = {
        'nprocx': int(sys.argv[7]),
        'nprocy': int(sys.argv[8]),
        'nprocz': int(sys.argv[9])
    }

    # Check if C++ combiner is available for massive speedup
    cpp_bin = find_cpp_combiner()
    if cpp_bin:
        cmd = [cpp_bin, "surf", prefix, str(step), str(grid['nox']), str(grid['noy']), str(grid['noz']),
               str(nprocxy), str(cap['nprocx']), str(cap['nprocy']), str(cap['nprocz'])]
        res = subprocess.run(cmd)
        if res.returncode == 0:
            sys.exit(0)

    # Fallback to pure Python
    nproc_per_cap = cap['nprocx'] * cap['nprocy'] * cap['nprocz']
    for i in range(nprocxy):
        cb = CombineSurf(grid)
        for n in range(i * nproc_per_cap + (cap['nprocz'] - 1), (i + 1) * nproc_per_cap, cap['nprocz']):
            crdfilename = f'{prefix}.coord.{n}'
            surffilename = f'{prefix}.surf.{n}.{step}'
            print('reading ' + surffilename)
            data = cb.readData(crdfilename, surffilename, grid, cap)
            cb.join(data, n, grid, cap)

        filename = f'{prefix}.surf{i}.{step}'
        print('writing ' + filename)
        cb.write(filename, grid, cb.saved)
