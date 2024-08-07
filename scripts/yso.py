#!/usr/bin/env python                                                                                                                                                                            

"""                                                                                                                                                                                                 
Generates a list of coordinations, luminosities of young stellar objects (YSO) from simulation snapshots.                                                                                            
                                                                                                                                                                                                     
Usage: yso.py <files> ... [options]                                                                                                                                                         
                                                                                                                                                                                                     
Options:                                                                                                                                                                                             
   -h --help                   Show this screen.    
   --age=<tau>                 Desired age of young stellar objects (YSO)
   --output_path               output path for data (defaults to /YSOobjects directory next to the snapshot)
   --num_jobs=<N>              Number of snapshots to process in parallel [default: 1]
"""                                                                                                                                                                                          
from os import mkdir
from os.path import isdir
import pathlib
import numpy as np
from docopt import docopt
import h5py     
import astropy.units as u  
from joblib import Parallel, delayed


np.random.seed(42)

options = docopt(__doc__)
if options["--age"]:
    TAU = float(options["--age"])
else:
    TAU = 0.5
if options["--output_path"]:
    if not isdir(options["--output_path"]):
        mkdir(options["--output_path"])
        OUTPATH = options["--output_path"]
else:
    OUTPATH = None

NUM_JOBS = int(options["--num_jobs"])


conversion_time = (1 * u.parsec / (u.meter/ u.second)).to(u.Myr)

def age(path):
    Xcords = []
    Ycords = []
    lum = []
    with h5py.File(path,"r") as F:
        Simulation_time = (F["Header"].attrs["Time"] * conversion_time).value
        X_Coordinates = (F["PartType5/Coordinates"][:,0])
        Y_Coordinates = (F["PartType5/Coordinates"][:,1])
        Luminosity = (F["PartType5/StarLuminosity_Solar"][:])
        Star_Formation  = (F["PartType5/ProtoStellarAge"][:] * conversion_time).value
        for i in range(len(Star_Formation)):
            t =  Simulation_time - Star_Formation[i]
            if t <= TAU:
                Xcords.append(X_Coordinates[i])
                Ycords.append(Y_Coordinates[i])
                lum.append(Luminosity[i])
            else:
                continue
    fname = path.split("/")[-1].replace(".hdf5", ".YSOobjects.hdf5")

    if OUTPATH:
        imgpath = OUTPATH + fname
    else:
        outdir ="/work2/10071/alexescamilla2244/frontera/output" + "/YSOobjects/"
    if not isdir(outdir):
        mkdir(outdir)
    imgpath = outdir + fname
    with h5py.File(imgpath, "w") as F:
        F.create_dataset("X_pc", data=Xcords)
        F.create_dataset("Y_pc", data=Ycords)
        F.create_dataset("Luminosity", data=lum)



def main():
    Parallel(n_jobs=NUM_JOBS)(
        delayed(age)(f) for f in options["<files>"]
    )


if __name__ == "__main__":
    main()
