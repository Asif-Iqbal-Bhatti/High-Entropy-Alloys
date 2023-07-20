#!/usr/bin/env python

import os, sys
import numpy as np
import py4vasp as p4
import matplotlib.pyplot as plt

print(p4.__version__)
print(dir(p4.data))

def ML_forces():
    os.system("cat ML_LOGFILE | grep BEEF > BEEF.dat")
    t1, beef = np.loadtxt("BEEF.dat",usecols=[1,3], unpack=True)
    t = p4.plot(t1, beef,
        xlabel="Time step",
        ylabel="Bayesian error",
        title="Bayesian error estimate of forces (max) (eV Angst^-1)")
    t.show()
    
def ML_RMSE():
    os.system("cat ML_LOGFILE | grep ERR > ERR.dat")
    t2, inerr = np.loadtxt("ERR.dat",usecols=[1,2], unpack=True)
    t = p4.plot(t2, inerr,
        xlabel="Time step",
        ylabel="RMSE",
        title="Root mean squared error of forces (eV Angst^-1)")
    t.show()


#calc = p4.Calculation.from_path(".")
#mycalc = p4.Calculation.from_path( "." )
#mycalc.structure[:].plot()