#!/usr/bin/env python3
import os
from pathlib import Path
import sys

__projectdir__ = Path(os.path.dirname(os.path.realpath(__file__)) + '/')

import numpy as np

def getinputdict(loglineareqs = True):
    inputdict = {}

    inputdict['paramssdict'] = {'GAMMA': 1, 'RHO': 0.05, 'ETA': 2, 'MU': 100, 'SIGMA': 6, 'PHI_pi': 1.5, 'RHO_A': 0.9, 'Abar': 1, 'Pistar': 1}

    inputdict['states'] = ['A']
    inputdict['controls'] = ['C', 'Pi', 'I', 'W', 'L', 'MC', 'Y']
    inputdict['irfshocks'] = ['A']

    # equations:{{{
    inputdict['equations'] = []

    if loglineareqs is True:
        inputdict['equations'].append('C_dot = 1 / GAMMA * (I - Pi)')
    else:
        inputdict['equations'].append('C_dot / C = 1 / GAMMA * (log(I) - log(Pi) - RHO)')

    if loglineareqs is True:
        inputdict['equations'].append('W - GAMMA * C = ETA * L')
    else:
        inputdict['equations'].append('W * C ** (-GAMMA) = L ** ETA')
    
    if loglineareqs is True:
        inputdict['equations'].append('Pi_dot = RHO * Pi - SIGMA * MC_ss / MU * MC')
    else:
        inputdict['equations'].append('Pi_dot / Pi = RHO * log(Pi) - 1 / MU * (SIGMA * MC - (SIGMA - 1))')
    
    if loglineareqs is True:
        inputdict['equations'].append('MC = W - A')
    else:
        inputdict['equations'].append('MC = W / A')
    
    if loglineareqs is True:
        inputdict['equations'].append('Y = A + L')
    else:
        inputdict['equations'].append('Y = A * L')
    
    if loglineareqs is True:
        inputdict['equations'].append('I = PHI_pi * Pi')
    else:
        inputdict['equations'].append('I = exp(RHO) * Pi_ss * (Pi / Pi_ss) ** PHI_pi')
    
    if loglineareqs is True:
        inputdict['equations'].append('C = Y')
    else:
        inputdict['equations'].append('C = Y')
    
    if loglineareqs is True:
        inputdict['equations'].append('A_dot = (RHO_A - 1) * A')
    else:
        inputdict['equations'].append('A_dot = (RHO_A - 1) * (A - 1)')
    
    # equations:}}}

    # steady state:{{{
    p = inputdict['paramssdict']

    p['A'] = p['Abar']
    p['Pi'] = p['Pistar']
    p['I'] = np.exp(p['RHO'] + np.log(p['Pi']))
    p['MC'] = (p['SIGMA'] - 1) / p['SIGMA'] + p['RHO'] * p['MU'] / p['SIGMA'] * np.log(p['Pi'])
    p['W'] = p['MC'] * p['A']
    p['L'] = p['W'] ** (1 / (p['GAMMA'] + p['ETA']))
    p['C'] = p['L']
    p['Y'] = p['C']
    # steady state:}}}

    if loglineareqs is True:
        inputdict['loglineareqs'] = True
    else:
        inputdict['logvars'] = inputdict['states'] + inputdict['controls']

    return(inputdict)


def check():
    inputdict_loglin = getinputdict(loglineareqs = True)
    inputdict_log = getinputdict(loglineareqs = False)
    sys.path.append(str(__projectdir__ / Path('submodules/dsge-perturbation/')))
    from dsge_continuous_func import checksame_inputdict_cont
    checksame_inputdict_cont(inputdict_loglin, inputdict_log)
    

def dsgefull():
    inputdict = getinputdict()

    sys.path.append(str(__projectdir__ / Path('submodules/dsge-perturbation/')))
    from dsge_continuous_func import continuouslineardsgefull
    continuouslineardsgefull(inputdict)


# Run:{{{1
check()
dsgefull()
