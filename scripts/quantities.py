#!/usr/bin/env python

"""

Generates Quantities For Data Analysis (NYSOs, M_Gas, Area, etc)

Usage: quantities.py <file1> <file2>  ... [options]

Options:

-h --help            Show this screen
--Cutoff=<Cutoff>    Cutoff density where quantities are calculated (defaults to 100)    
--output_path        output path for data (defaults to /E_FF<CUTOFF>_Quantities in parent directory)
--BoxSize=<Bx_Size>  Box size of simulation (defaults to 100)  
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

options = docopt(__doc__)
if options["--output_path"]:
    if not isdir(options["--output_path"]):
        mkdir(options["--output_path"])
        OUTPATH = options["--output_path"]
else:
    OUTPATH = None

CUTOFF = 100

def Quantities(path1, path2, CUTOFF, Box_Size, OUTPATH):
    #Getting YSO Coordinates                                                                                                                                                       
    with h5py.File(path1, "r") as A:
        YSO_X = A["X_pc"][:]
        YSO_Y = A["Y_pc"][:]

    with h5py.File(path2, "r") as F:
        X = F["X_pc"][:]
        Y = F["Y_pc"][:]
        Surface_Den = F["SurfaceDensity_Msun_pc2"][:]
        time =(F["Header"].attrs["Time"] * u.pc / (u.meter/u.second)).to(u.Myr)
    #lists for the quantities that will be displayed

    N_YSOS = []
    Mass_Gass = []
    Area = []
    tff = []
    t = []
    #Flattening out Density and Dust Coordinates                                                                                                                                    
    new_den = (Surface_Den.flatten()) * (u.solMass / u.pc**2)
    dust_cords = np.c_[X.flatten(), Y.flatten()]

    # Finding Surface Density of YSO From Nearest Pixel                                                                                                                           
                                                                                                                                                                                    
    
    targets = np.c_[YSO_X.flatten(), YSO_Y.flatten()]
    T = KDTree(dust_cords)
    dist, idx = T.query(targets)
    df = pd.DataFrame({ "YSOs" : targets[:,0], "YSO Surface Den" : new_den[idx]})

    #Calculating Star Formation Rate                                                                                                                                               
    NumYSOs = df["YSOs"]
    surface_unit = 1 * (u.solMass / u.pc**2)
    
    for i in np.arange(65,2500):

        yso_surface_den = df["YSO Surface Den"]
        num_yso_at_surface_den = (yso_surface_den > i).sum()
        StarRate = (((num_yso_at_surface_den * 0.5) /  0.5) * (u.solMass / u.Myr))
        L = Box_Size  / 5
        Mask = new_den > i * surface_unit
        A_i= (L / 1024)**2 * u.pc**2
        A = (np.sum(Mask) *  A_i).value

        #Calculating Mass_Gas                                                                                                                                                           
        M_Gas = (np.sum(A_i * new_den[Mask])).value

        #Calculating e_ff                                                                                                                                                               
        new_G = G.to(u.pc**3/(u.solMass * u.Myr**2)).value

        t_ff_num = np.sum(A**(3/2)) * np.sqrt(np.pi)
        t_ff_denom = (8 * new_G * M_Gas)
        T_ff = np.sqrt(t_ff_num/t_ff_denom)
    
    

    
        #Appending All Quantities
        N_YSOS.append(num_yso_at_surface_den)
        Mass_Gass.append(M_Gas)
        Area.append(A)
        tff.append(T_ff)
        t.append(time)
    #Creatiing File and Datasets
    fname =  os.path.basename(path1).replace(".YSOobjects.hdf5", ".quantities.hdf5")
    if OUTPATH:
        imgpath = OUTPATH + fname
    else:
       outdir = "/work2/10071/alexescamilla2244/frontera/CASSI_Project-2024/output" + "/E_FF" + str(CUTOFF) + "_Quantities/"
    if not isdir(outdir):
        mkdir(outdir)
    imgpath = outdir + fname

    with h5py.File(imgpath, "w") as F:
        F.create_dataset("NYSOs", data=N_YSOS)
        F.create_dataset("M_Gas", data=Mass_Gass)
        F.create_dataset("Area", data=Area)
        F.create_dataset("T_ff", data=tff)
        F.create_dataset("t", data=t)
    
        
def main():
    options = docopt(__doc__)
    if options["--output_path"]:
        if not isdir(options["--output_path"]):
            mkdir(options["--output_path"])
            OUTPATH = options["--output_path"]
    else:
        OUTPATH = None

    if options["--Cutoff"]:
        CUTOFF = options["--Cutoff"]
    else:
        CUTOFF = 100
    if options["--BoxSize"]:
        Box_Size = options["--BoxSize"]
    else:
        Box_Size = 100
        
    
    file1= options["<file1>"]
    file2 = options["<file2>"]


    if isinstance(file1, list):
        file1 = file1[0]
    if isinstance(file2, list):
        file2 = file2[0]

    
    Quantities(file1,file2, CUTOFF, Box_Size, OUTPATH)
    
if __name__ == "__main__":
    main()
