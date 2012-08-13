import numpy as np
from rootfinder import fzero
from constants import constants as c
import matplotlib.cbook as cbook
import numpy.testing as test

def esat(Temp):
    """
    esat(Temp)

    Calculates the saturation water vapor pressure over a flat
    surface of water at temperature 'T'.

    Parameters
    - - - - - -
    Temp : float or array_like
        Temperature of parcel (K).

    Returns
    - - - -
    esatOut : float or list
        Saturation water vapour pressure (Pa).

    Examples
    - - - - -
    >>> test.assert_almost_equal(esat(300.),3534.5196,decimal=3)
    >>> np.allclose(esat([300., 310.]),[3534.519, 6235.532])
    True

    References
    - - - - - -
    Emanuel 4.4.14 p. 117
      
    """
    # determine if Temp has been input as a vector
    is_scalar=True
    if cbook.iterable(Temp):
        is_scalar = False
    Temp=np.atleast_1d(Temp)
    Tc = Temp - c.Tc
    esatOut = 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))
    # if T is a vector
    if is_scalar:
        esatOut = esatOut[0]
    return esatOut

def LCLfind(Td, T, p):
    """
    LCLfind(Td, T, p)

    Finds the temperature and pressure at the lifting condensation
    level (LCL) of an air parcel.

    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa)

    Returns
    - - - -
    Tlcl : float
        Temperature at the LCL (K).
    plcl : float
        Pressure at the LCL (Pa).

    Raises
    - - - -
    NameError
        If the air is saturated at a given Td and T (ie. Td >= T)
    
    Examples
    - - - - -
    >>> [Tlcl, plcl] =  LCLfind(280., 300., 8.e4)
    >>> print [Tlcl, plcl]
    [275.76250387361404, 59518.928699453245]
    >>> LCLfind(300., 280., 8.e4)
    Traceback (most recent call last):
        ...
    NameError: parcel is saturated at this pressure

    References
    - - - - - -
    Emanuel 4.6.24 p. 130 and 4.6.22 p. 129
    
    """
    hit = Td >= T;
    if hit is True:
        raise NameError('parcel is saturated at this pressure');

    e = esat(Td);
    ehPa = e * 0.01; #Bolton's formula requires hPa.
    # This is is an empircal fit from for LCL temp from Bolton, 1980 MWR.
    Tlcl = (2840. / (3.5 * np.log(T) - np.log(ehPa) - 4.805)) + 55.;

    r = c.eps * e / (p - e);
    #disp(sprintf('r=%0.5g',r'))
    cp = c.cpd + r * c.cpv;
    logplcl = np.log(p) + cp / (c.Rd * (1 + r / c.eps)) * \
              np.log(Tlcl / T);
    plcl = np.exp(logplcl);
    #disp(sprintf('plcl=%0.5g',plcl))

    return Tlcl, plcl


def wsat(Temp, press):
    """
    wsat(Temp, press)

    Calculates the saturation vapor mixing ratio of an air parcel.

    Parameters
    - - - - - -
    Temp : float or array_like
        Temperature in Kelvin.
    press : float or array_like
        Pressure in Pa.

    Returns
    - - - -
    theWs : float or array_like 
        Saturation water vapor mixing ratio in (kg/kg).

    Raises
    - - - -
    IOError
        If both 'Temp' and 'press' are array_like.

    Examples
    - - - - -
    >>> test.assert_almost_equal(wsat(300, 8e4),0.02875,decimal=4)
    >>> test.assert_array_almost_equal(wsat([300,310], 8e4),[0.0287, 0.0525],decimal=4)
    >>> test.assert_array_almost_equal(wsat(300, [8e4, 7e4]),[0.0287, 0.0330],decimal=4)
    >>> wsat([300, 310], [8e4, 7e4])
    Traceback (most recent call last):
        ...
    IOError: Can't have two vector inputs.

    """
    is_scalar_temp=True
    if cbook.iterable(Temp):
        is_scalar_temp = False
    is_scalar_press=True
    if cbook.iterable(press):
        is_scalar_press = False
    Temp=np.atleast_1d(Temp)
    press=np.atleast_1d(press)
    if (np.size(Temp) !=1) and (np.size(press) != 1):
        raise IOError, "Can't have two vector inputs."
    es = esat(Temp);
    theWs=(c.eps * es/ (press - es))
    theWs[theWs > 0.060]=0.06
    theWs[theWs < 0.0] = 0.
    if is_scalar_temp and is_scalar_press:
        theWs=theWs[0]
    return theWs

def thetaep(Td, T, p):
    """
    thetaep(Td, T, p)

    Calculates the pseudo equivalent potential temperature of a
    parcel. 


    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    thetaepOut : float
        Pseudo equivalent potential temperature (K).

    """
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = LCLfind(Td, T, p);
        wv = wsat(Td, p);
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Tlcl = T;
        wv = wsat(T, p);
        
    Pd = 1.0*(c.eps/(wv+c.eps))*p
    thetaval = T*(c.p0/Pd)**(c.Rd/c.cpd)    
    ws = wsat(T, p)
    # $$$   disp('inside theate')
    # $$$   [Td,T,wv]
    thetaepOut = thetaval * np.exp((c.lv0*ws)/(c.cpd*T))
    #
    # peg this at 450 so rootfinder won't blow up
    #

    if (thetaepOut > 450.):
        thetaepOut = 450

    return thetaepOut

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

