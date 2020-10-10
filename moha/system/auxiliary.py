import numpy as np
import scipy.special as sc

def fact2(n):
    """Double factorial function.

    Parameters
    ----------
    n : int
        Input number of double factorial function

    Raises
    ------
    TypeError
        If input parameter is not an integer.
    """
    if not isinstance(n, int):
        raise TypeError("Input parameter must be an integer")
    if n <= 0:
        return 1
    else:
        return n * fact2(n-2)

def boys(n,t):
    """Boys function for the calculation of coulombic integrals.

    Parameters
    ----------
    n : int
        Order of boys function

    t : float
        Varible for boys function.

    Raises
    ------
    TypeError
        If boys function order is not an integer.

    ValueError
        If boys function order n is not a none negative number.
    """
    if not isinstance(n, int):
        raise TypeError("Boys function order n must be an integer")
    if n < 0:
        raise ValueError("Boys function order n must be a none negative number")    
    if not isinstance(t, float):
        raise TypeError("Boys function varible t must be integer or float")
    return sc.hyp1f1(n+0.5,n+1.5,-t)/(2.0*n+1.0)

def gaussian_product_center(a,A,b,B):
    """
    """
    A = np.array(A)
    B = np.array(B)
    return (a*A+b*B)/(a+b)

def eint(a,b,c,d):
    """
    """
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return int(abcd)
