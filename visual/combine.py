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
Combine the pasted Citcom Data (Accelerated via Modern C++ citcoms_combine)

usage: combine.py datafile timestep nodex nodey nodez ncap nprocx nprocy nprocz [opts]
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

class Combine(object):

    def __init__(self, nodex, nodey, nodez, nprocx, nprocy, nprocz):
        self.nodex = nodex
        self.nodey = nodey
        self.nodez = nodez

        self.nprocx = nprocx
        self.nprocy = nprocy
        self.nprocz = nprocz

        # data storage
        self.saved = [''] * (nodex * nodey * nodez)
        return

    def readData(self, filename, headlines=1):
        with open(filename, 'r') as fp:
            for _ in range(headlines):
                fp.readline()
            return fp.readlines()

    def join(self, data, me):
        nprocx, nprocy, nprocz = self.nprocx, self.nprocy, self.nprocz
        mylocz = me % nprocz
        mylocx = ((me - mylocz) // nprocz) % nprocx
        mylocy = (((me - mylocz) // nprocz - mylocx) // nprocx) % nprocy

        nodex, nodey, nodez = self.nodex, self.nodey, self.nodez
        mynodex = 1 + (nodex - 1) // nprocx
        mynodey = 1 + (nodey - 1) // nprocy
        mynodez = 1 + (nodez - 1) // nprocz

        if not len(data) == mynodex * mynodey * mynodez:
            raise ValueError("incorrect data size")

        mynxs = (mynodex - 1) * mylocx
        mynys = (mynodey - 1) * mylocy
        mynzs = (mynodez - 1) * mylocz

        n = 0
        for i in range(mynys, mynys + mynodey):
            for j in range(mynxs, mynxs + mynodex):
                for k in range(mynzs, mynzs + mynodez):
                    m = k + j * nodez + i * nodex * nodez
                    self.saved[m] = data[n]
                    n += 1

    def join_EL(self, data, me):
        nprocx, nprocy, nprocz = self.nprocx, self.nprocy, self.nprocz
        mylocz = me % nprocz
        mylocx = ((me - mylocz) // nprocz) % nprocx
        mylocy = (((me - mylocz) // nprocz - mylocx) // nprocx) % nprocy

        nodex, nodey, nodez = self.nodex, self.nodey, self.nodez
        mynodex = 1 + (nodex - nprocx) // nprocx
        mynodey = 1 + (nodey - nprocy) // nprocy
        mynodez = 1 + (nodez - nprocz) // nprocz

        if not len(data) == mynodex * mynodey * mynodez:
            raise ValueError("incorrect data size")

        mynxs = mynodex * mylocx
        mynys = mynodey * mylocy
        mynzs = mynodez * mylocz

        n = 0
        for i in range(mynys, mynys + mynodey):
            for j in range(mynxs, mynxs + mynodex):
                for k in range(mynzs, mynzs + mynodez):
                    m = k + j * nodez + i * nodex * nodez
                    self.saved[m] = data[n]
                    n += 1

    def write(self, filename):
        with open(filename, 'w') as fp:
            fp.write(f'{self.nodex} x {self.nodey} x {self.nodez}\n')
            fp.writelines(self.saved)

##############################################

def combine(prefix, *args):
    # Support both (prefix, step, ...) and (prefix, opts, step, ...)
    if len(args) == 9:
        step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz, opts = args
    elif len(args) == 8:
        step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz = args
        opts = 'coord,velo,visc'
    else:
        # legacy positional signature with opts second
        opts, step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz = args

    # Check if high-performance C++ citcoms_combine binary is available
    cpp_bin = find_cpp_combiner()
    if cpp_bin:
        cmd = [cpp_bin, "volume", str(prefix), str(step), str(nodex), str(nodey), str(nodez),
               str(ncap), str(nprocx), str(nprocy), str(nprocz), str(opts)]
        res = subprocess.run(cmd)
        if res.returncode == 0:
            return

    # Fallback to pure Python implementation
    combined_files = []
    nproc_per_cap = nprocx * nprocy * nprocz
    for i in range(ncap):
        cb = Combine(nodex, nodey, nodez, nprocx, nprocy, nprocz)
        for n in range(i * nproc_per_cap, (i+1) * nproc_per_cap):
            filename = f'{prefix}.{opts}.{n}.{step}.pasted'
            print('reading ' + filename)
            data = cb.readData(filename, 0)
            if opts == 'heating':
                cb.join_EL(data, n)
            else:
                cb.join(data, n)

        if opts == 'coord,velo,visc':
            filename = f'{prefix}.cap{i:02d}.{step}'
        elif opts == 'heating':
            filename = f'{prefix}.ele{i:02d}.{step}'
        else:
            filename = f'{prefix}.opt{i:02d}.{step}'

        print('writing ' + filename)
        cb.write(filename)
        combined_files.append(filename)

    return combined_files


if __name__ == '__main__':
    if len(sys.argv) not in (10, 11):
        print(__doc__)
        sys.exit(1)

    prefix = sys.argv[1]
    step = int(sys.argv[2])
    nodex = int(sys.argv[3])
    nodey = int(sys.argv[4])
    nodez = int(sys.argv[5])
    ncap = int(sys.argv[6])
    nprocx = int(sys.argv[7])
    nprocy = int(sys.argv[8])
    nprocz = int(sys.argv[9])
    opts = sys.argv[10] if len(sys.argv) == 11 else 'coord,velo,visc'

    combine(prefix, step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz, opts)
