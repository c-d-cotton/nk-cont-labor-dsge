#!/usr/bin/env python3
# PYTHON_PREAMBLE_START_STANDARD:{{{

# Christopher David Cotton (c)
# http://www.cdcotton.com

# modules needed for preamble
import importlib
import os
from pathlib import Path
import sys

# Get full real filename
__fullrealfile__ = os.path.abspath(__file__)

# Function to get git directory containing this file
def getprojectdir(filename):
    curlevel = filename
    while curlevel is not '/':
        curlevel = os.path.dirname(curlevel)
        if os.path.exists(curlevel + '/.git/'):
            return(curlevel + '/')
    return(None)

# Directory of project
__projectdir__ = Path(getprojectdir(__fullrealfile__))

# Function to call functions from files by their absolute path.
# Imports modules if they've not already been imported
# First argument is filename, second is function name, third is dictionary containing loaded modules.
modulesdict = {}
def importattr(modulefilename, func, modulesdict = modulesdict):
    # get modulefilename as string to prevent problems in <= python3.5 with pathlib -> os
    modulefilename = str(modulefilename)
    # if function in this file
    if modulefilename == __fullrealfile__:
        return(eval(func))
    else:
        # add file to moduledict if not there already
        if modulefilename not in modulesdict:
            # check filename exists
            if not os.path.isfile(modulefilename):
                raise Exception('Module not exists: ' + modulefilename + '. Function: ' + func + '. Filename called from: ' + __fullrealfile__ + '.')
            # add directory to path
            sys.path.append(os.path.dirname(modulefilename))
            # actually add module to moduledict
            modulesdict[modulefilename] = importlib.import_module(''.join(os.path.basename(modulefilename).split('.')[: -1]))

        # get the actual function from the file and return it
        return(getattr(modulesdict[modulefilename], func))

# PYTHON_PREAMBLE_END:}}}

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
    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsge_continuous_func.py'), 'checksame_inputdict_cont')(inputdict_loglin, inputdict_log)
    

def dsgefull():
    inputdict = getinputdict()

    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsge_continuous_func.py'), 'continuouslineardsgefull')(inputdict)


# Run:{{{1
check()
dsgefull()
