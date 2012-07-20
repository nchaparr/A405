"""This is the docstring for the calcAdiabat.py module."""
import site
site.addsitedir('C:\Users\Den\mya405\python\\thermlib')
import numpy as np
from findTmoist import findTmoist

def calcAdiabat(press0, thetae0, topPress):
    """
    
    Calculates the temperature-pressure coordinates of a moist adiabat.
    
    Parameters
    - - - - - -
    press0: the initial pressure (Pa)
    thetae0: the equivalent potential temperature (K) of the adiabat
    topPress: the final pressure (Pa)
    
    Returns
    - - - - - -
    (pressVals, tempVals): pressVals (Pa) and tempVals (K) are type ndarray
                           and are the coordinates of the thetae0 adiabat
                
    
    Tests
    - - - - -
    >>> p,T = calcAdiabat(800*100, 300, 1000*100)
    >>> len(p)
    50
    >>> len(T)
    50
    

    """
    
    pressVals = np.linspace(press0, topPress, 50)
    
    tempVals = findTmoist(thetae0, pressVals)
    
    return pressVals, np.asarray(tempVals)

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

     
         
        
    
    
    
    