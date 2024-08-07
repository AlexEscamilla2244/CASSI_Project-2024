#!/usr/bin/env python      

"""

Generates a range of values of star formation efficiency 

Usage: e_ff.py <file1> <file2> ... [options]

Options:

-h --help                  Show this screen>
--BoxSize=<BX>             Box size of the simulation (defaults to 100)
--Starting_Range=<R_i>     Starting range of density cutoff (defaults to 65)
--Ending_Range=<R_f>       Ending range of density cutoff (defaults to 2500)
--output_path              output path for data (defaults to /E_FF_range in parent directory )         
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


def SFE_values(path1,path2, Box_Size, R_I, R_F, OUTPATH):
    print(path1)
    #Getting YSO Coordinates
    with h5py.File(path1) as A:
        YSO_X = (A["X_pc"][:])
        YSO_Y = (A["Y_pc"][:])

    #Getting Dust Coordinates
    with h5py.File(path2) as F:
        X = (F["X_pc"][:])
        Y = (F["Y_pc"][:])
        Surface_Den = F["SurfaceDensity_Msun_pc2"][:]

        #Flattening out Density and Dust Coordinates
        new_den = (Surface_Den.flatten()) * (u.solMass / u.pc**2)
        dust_cords = np.c_[X.flatten(), Y.flatten()]

    # Finding Surface Density of YSO From Nearest Pixel
    targets = np.c_[YSO_X.flatten(), YSO_Y.flatten()]
    T = KDTree(dust_cords)
    dist, idx = T.query(targets)
    df = pd.DataFrame({ "YSOs" : targets[:,0], "YSO Surface Den" : new_den[idx]})

   
    # Calculating Star Formation Rate
    NumYSOs = df["YSOs"]
    yso_surface_den = df["YSO Surface Den"]

    e_ff_range = []
    for i in np.arange(R_I,R_F,1):
        
        StarRate = (((yso_surface_den > i).sum() * 0.5) /  0.5) * (u.solMass / u.Myr)

    
        # Calculating Area 
        surface_unit = 1 * (u.solMass / u.pc**2)
        L = Box_Size / 5
        Mask = new_den > i * surface_unit
        A_i= (L / 1024)**2 * u.pc**2
        A = np.sum(Mask) *  A_i

        #Calculating Sigma_Gas
        M_Gas = np.sum(A_i * new_den[Mask] ) 
    
        #Calculating e_ff
        new_G = G.to(u.pc**3/(u.solMass * u.Myr**2))
    
        T_ff = (new_G* (M_Gas / np.sum(A** (3/2))))**(-1/2)
        star_efficiency = (StarRate / (M_Gas / T_ff)).value
        e_ff_range.append(star_efficiency)
    

    fname = path1.split("/")[-1].replace(".YSOobjects.hdf5", ".e_ff_range.hdf5")   
    if OUTPATH:
        imgpath = OUTPATH + fname
    else:
        outdir = "/work2/10071/alexescamilla2244/frontera/output" + "/E_FF/"
    if not isdir(outdir):
        mkdir(outdir)
    imgpath = outdir + fname
    with h5py.File(imgpath, "w") as F:
        F.create_dataset("SFE_Values", data=e_ff_range)


   


def main():
    options = docopt(__doc__)

    if options["--BoxSize"]:
        Box_Size = int(options["--BoxSize"])
    else:
        Box_Size = 100

    if options["--Starting_Range"]:
        R_I = int(options["--Starting_Range"])
    else:
        R_I = 65

    if options["--Ending_Range"]:
        R_F = int(options["--Ending_Range"])
    else:
        R_F = 2500

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
    
    SFE_values(file1, file2, Box_Size, R_I, R_F, OUTPATH)

if __name__ == "__main__":
    main()




    
    
    
    










        
        
    
