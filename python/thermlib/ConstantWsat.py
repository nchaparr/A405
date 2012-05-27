#!/usr/bin/env python

import sys
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

"""This script finds the Temperature which returns zero for the
   following function:
   Wsat0Temp = Wsat - Wsat0
   where Wsat = 0.622*Esat/(press - Esat)
   and Esat = 611.2*exp(17.67*Temp(C)/(Temp(C) +243.3))

   Inputs are press in Pa and Temp in K (converted to C)

"""

def WsatShift(Temp, press, Wsat0):
    Esat = 611.2*np.exp(17.67*1.0*(Temp - 273.15)/((Temp - 273.15) + 243.5))
    Wsat = 0.622*(1.0*Esat/(press -Esat))
    if Wsat > 0.060:
        Wsat = 0.060
    elif Wsat < 0:
            Wsat = 0
    WsatShift = Wsat - Wsat0
    return WsatShift

def WsatTemp(press, Wsat0):
    WsatTemp = optimize.zeros.brenth(WsatShift, 200, 330, (press, Wsat0))
    return WsatTemp

if __name__=="__main__":

    Wsat0 = .01
    p = 80000 
    print WsatShift(200, p, Wsat0), WsatShift(330, p, Wsat0)
    print WsatTemp(p, Wsat0)
