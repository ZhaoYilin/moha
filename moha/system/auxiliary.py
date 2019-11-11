import numpy as np
import scipy.special as sc

def fact2(n):
    """
    """
    if n <= 0:
        return 1
    else:
        return n * fact2(n-2)

def boys(n,T):
    """
    """
    return sc.hyp1f1(n+0.5,n+1.5,-T)/(2.0*n+1.0)

def gaussian_product_center(a,A,b,B):
    """
    """
    A = np.array(A)
    B = np.array(B)
    return (a*A+b*B)/(a+b)

def eint(a,b,c,d):
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return int(abcd)