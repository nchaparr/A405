import scipy.optimize

from findWvWl import findWvWl
from my_thetaep import thetaep

def t_thetaep(thetaepVal, wT, p):
    
    """
    derived from tinvert
    
    """
    if p > 1.e5:
        raise NameError, \
              'expecting pressure level less than 100000 Pa'
    # The temperature has to be somewhere between 'thetal'
    # (T at surface) and -40 deg. C (no ice).
    
    T = scipy.optimize.zeros.brenth(Tchange, 233.15, thetaepVal, \
                                     (thetaepVal, wT, p));
    [wv, wl] = findWvWl(T, wT, p);
    return T, wv, wl


def Tchange(Tguess, thetaepVal, wT, p):
    [wv, wl] = findWvWl(Tguess, wT, p);
    # Iterate on Tguess until this function is zero to within the
    # default tolerance in brenth.
    return thetaepVal - thetaep(wv, Tguess, p);

def _test():
    #import doctest
    #doctest.testmod()
    thetaepVal = float(raw_input('thetaep?'))
    wT = float(raw_input('total water?'))
    p = float(raw_input('pressure?'))
    print t_thetaep(thetaepVal, wT, p)

if __name__ == "__main__":
    _test()

