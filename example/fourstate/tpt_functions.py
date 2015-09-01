#!/bin/env python

import numpy as np

def gen_path_lengths(keys, J, pfold, flux, FF, UU):
    nkeys = len(keys)
    I = [x for x in range(nkeys) if x not in FF+UU]
    Jnode = []
    # calculate flux going through nodes
    for i in range(nkeys):
        Jnode.append(np.sum([J[i,x] for x in range(nkeys) \
                             if pfold[x] < pfold[i]]))
    # define matrix with edge lengths
    Jpath = np.zeros((nkeys, nkeys), float)
    for i in UU:
        for j in I + FF:
            if J[j,i] > 0:
                Jpath[j,i] = np.log(flux/J[j,i]) + 1e-9 # I add 1e-9
    for i in I:
        for j in [x for x in FF+I if pfold[x] > pfold[i]]:
            if J[j,i] > 0:
                Jpath[j,i] = np.log(Jnode[i]/J[j,i]) + 1e-9 # I add 1e-9
    return Jnode, Jpath
