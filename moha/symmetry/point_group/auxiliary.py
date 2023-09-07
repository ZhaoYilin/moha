import numpy as np


def divisors(n:int)->list:
    """Get the divisors of a given number, start from 2.

    Parameters
    ----------
    n : int
        The given number.

    Returns
    -------
    result : List
        A list of divisors.
    """
    result = set()
    for i in range(2, int(n**0.5)+1):
        if n % i == 0:
            result.add(i)
            result.add(n//i)
    result = list(result)
    return result

def angle(a:list, b:list)->float:
    """Evalute the angle between two vectors.

    Parameters
    ----------
    a : List[float,float,float]
        First vector.

    b : List[float,float,float]
        Second vector.

    Returns
    -------
    rad : float
        Angle with rad unit.
    """    
    inner = np.inner(a, b)
    norms = np.linalg.norm(a) * np.linalg.norm(b)

    cos = inner / norms
    rad = np.arccos(np.clip(cos, -1.0, 1.0))
    return rad
