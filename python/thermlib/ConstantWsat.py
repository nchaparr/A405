#!/usr/bin/env python

import sys
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

"""This script finds the Temperature which returns zero for the
   following function:
   Wsat0Temp = Wsat - Wsat0
   where Wsat = 0.644*Esat/(press - Esat)
   and Esat = 6.112*exp(17.67*Temp/(Temp +243.3))
"""

def WsatShift(Temp, press, Wsat0):
    Esat = 6.11*np.exp(17.67*1.0*(Temp - 273.15)/((Temp - 273.15) + 243.5))
    Wsat = 0.644*(1.0*Esat/(press -Esat))
    if Wsat > 0.060:
        Wsat = 0.060
    elif Wsat < 0:
            Wsat = 0
    WsatShift = Wsat - Wsat0
    return WsatShift

def WsatTemp(press, Wsat0):
    WsatTemp = optimize.zeros.brenth(WsatShift, 200, 400, (press, Wsat0))
    return WsatTemp

if __name__=="__main__":
    
    Wsat0 = .01
    p = 80000 
   
    print WsatTemp(p, Wsat0)
