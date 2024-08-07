#!/usr/bin/env python                                                                                                                                                               

"""                                                                                                                                                                                 
                                                                                                                                                                                    
Generates a value of star formation efficiency                                                                                                                            
                                                                                                                                                                                    
Usage: e_ff_value.py <file1> <file2> ... [options]                                                                                                                                        
                                                                                                                                                                                    
Options:                                                                                                                                                                            
                                                                                                                                                                                    
-h --help                  Show this screen>                                                                                                                                        
--BoxSize=<BX>             Box size of the simulation (defaults to 100)                                                                                                             
--Cutoff=<cutoff>          Cutoff density where star formation efficiency is calculated (defaults to 100)                                                                                                                                                                                              
--output_path              output path for data (defaults to /E_FFvalue in file 1 parent directory)                                                                                       
"""

from scipy.spatial import KDTree
import pandas as pd
from astropy.constants import G
import astropy.units as u
import os
import h5py
import numpy as np
from docopt import docopt
from joblib import Parallel, delayed
import pathlib
from os import mkdir
from os.path import isdir
import argparse

np.random.seed(42)

def SFE_values(path1,path2,CUTOFF, Box_Size, OUTPATH):
    #Getting YSO Coordinates                                                                                                                                                       
    with h5py.File(path1) as A:
        YSO_X = A["X_pc"][:]
        YSO_Y = A["Y_pc"][:]
        
    with h5py.File(path2) as F:
        X = F["X_pc"][:]
        Y = F["Y_pc"][:]
        Surface_Den = F["SurfaceDensity_Msun_pc2"][:]
        
    #Flattening out Density and Dust Coordinates                                                                                                                              
    new_den = (Surface_Den.flatten()) * (u.solMass / u.pc**2)
    dust_cords = np.c_[X.flatten(), Y.flatten()]

    # Finding Surface Density of YSO From Nearest Pixel                                                                                                                             
    e_ff = []
    targets = np.c_[YSO_X.flatten(), YSO_Y.flatten()]
    T = KDTree(dust_cords)
    dist, idx = T.query(targets)
    df = pd.DataFrame({ "YSOs" : targets[:,0], "YSO Surface Den" : new_den[idx]})

    #Calculating Star Formation Rate                                                                                                                                               
    NumYSOs = df["YSOs"]
    yso_surface_den = df["YSO Surface Den"]
    StarRate = (((yso_surface_den > CUTOFF).sum() * 0.5) /  0.5) * (u.solMass / u.Myr)

    # Calculating Area                                                                                                                                                          
    surface_unit = 1 * (u.solMass / u.pc**2)
    L = Box_Size / 5
    Mask = new_den > CUTOFF * surface_unit
    A_i= (L / 1024)**2 * u.pc**2
    A = np.sum(Mask) *  A_i

    #Calculating Mass_Gas                                                                                                                                                      
    M_Gas = np.sum(A_i * new_den[Mask] )

    #Calculating e_ff                                                                                                                                                           
    new_G = G.to(u.pc**3/(u.solMass * u.Myr**2))

    T_ff = (new_G* (M_Gas / np.sum(A** (3/2))))**(-1/2)
    star_efficiency = (StarRate / (M_Gas / T_ff)).value
    e_ff.append(star_efficiency)
    
    fname = path1.split("/")[-1].replace(".YSOobjects.hdf5", ".e_ff" + str(CUTOFF) + "value.hdf5")
    if OUTPATH:
        imgpath = OUTPATH + fname
    else:
       outdir = str(pathlib.Path(path1).parent.resolve()) + "/E_FF" + str(CUTOFF) + "_values/"
    if not isdir(outdir):
        mkdir(outdir)
    imgpath = outdir + fname
    
    with h5py.File(imgpath, "w") as F:
        F.create_dataset("SFE", data=e_ff)

def main():
    options = docopt(__doc__)

    if options["--BoxSize"]:
        Box_Size = int(options["--BoxSize"])
    else:
        Box_Size = 100

    if options["--Cutoff"]:
        CUTOFF = int(options["--Cutoff"])
    else:
        CUTOFF = 500

    if options["--output_path"]:
        if not isdir(options["--output_path"]):
            mkdir(options["--output_path"])
            OUTPATH = options["--output_path"]
    else:
        OUTPATH = None

    file1= options["<file1>"]
    file2 = options["<file2>"]
    
    
    if isinstance(file1, list):
        file1 = file1[0]
    if isinstance(file2, list):
        file2 = file2[0]
    
    SFE_values(file1,file2, CUTOFF, Box_Size, OUTPATH)

if __name__ == "__main__":
    main()
