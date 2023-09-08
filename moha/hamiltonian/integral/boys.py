from math import exp, gamma

def boys(v,x):
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
    if x <= 25:
        i = 0
        ans = 0
        g_v = gamma(v + 0.5)
        while True:
            seq = (g_v / gamma(v + i + 1.5)) * x**i
            if seq < 1e-10:
                break
            ans += seq
            i += 1
        ans *= (1/2) * exp(-x)
        return ans

    # Approximation of the boys function for large x
    elif x > 25:
        i = 0
        ans = 0
        g_v = gamma(v + 0.5)
        while True:
            seq = (g_v / gamma(v - i + 1.5)) * x**(-i)
            if seq < 1e-10:
                break
            ans += seq
            i += 1
        ans *= (1/2) * exp(-x)
        ans = (g_v / (2*x**(v + 0.5))) - ans
        return ans


def boys_recursion(N, x, f_N):
    """Returns the answer to the boys function f_{v - 1}(x) using the
    answer for the boys function f_{v}(x).

    Parameters
    ----------
    v : {int, float}
    x : {int, float}
    f_v : float

    Returns
    -------
    result : float
        The boys function f_{N - 1}(x).
    """
    result =  (exp(-x) + 2 * x * f_N) / (2 * N - 1)
    return result
