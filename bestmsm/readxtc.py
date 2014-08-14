#!/usr/bin/env python

import argparse
import os
import string
import sys
import numpy as np
import PDB
import xtc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read XTC trajectory files')
    parser.add_argument('-s', metavar='file_tpr', required=True, help='Gromacs run input file',
            dest='file_tpr')
    parser.add_argument('-p', metavar='file_pdb', required=True, help='PDB file',
            dest='file_pdb')
    parser.add_argument('-f', metavar='file_xtc', required=True, nargs='*',
            help='Gromacs XTC trajectory files', dest='file_xtc')
    args = parser.parse_args()

    A = PDB.LoadPDB(args.file_pdb)

    A["XYZList"] = []
    for c in xtc.XTCReader(args.file_xtc):
        A["XYZList"].append(np.array(c.coords).copy())
    i = 0
    for c in xtc.XTCReader(args.file_xtc):
        if i==0:
            ConfShape=np.shape(c.coords)
        i=i+1
    Shape=np.array((i,ConfShape[0],ConfShape[1]))
    print A
    print Shape
