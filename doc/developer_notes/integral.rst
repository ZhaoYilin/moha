Molecular Integrals
===================
Contains the information about how the integrals are calculated. In the equations in this section the following definitions is used.

.. math::
    p=a+b
    μ=aba+b
    Px=aAx+bBxp
    XAB=Ax−Bx

Here a and b are Gaussian exponent factors. Ax and Bx are the position of the Gaussians in one dimension. Further the basisset functions is of Gaussian kind and given as:

.. math::
    ϕA(r)=N(x−Ax)l(y−Ay)m(z−Az)nexp(−ζ(r⃗ −A⃗ )2)

with a normalization constant given as:

.. math::
    N=(2απ)3/4[(8α)l+m+nl!m!n!(2l)!(2m)!(2n)!]

Boys Function
-------------
The Boys function is given as:
.. math::
    Fn(x)=∫10exp(−xt2)t2ndt

FUNCTION:

MIcython.boys(m,T)
return value
Input:

m, subscript of the Boys function
T, argument of the Boys function
Output:

value, value corrosponding to given m and T

Expansion coefficients
----------------------
The expansion coefficient is found by the following recurrence relation:

Ei,jt=0,t<0ort>i+j
Ei+1,jt=12pEi,jt−1+XPAEi,jt+(t+1)Ei,jt+1
Ei,j+1t=12pEi,jt−1+XPBEi,jt+(t+1)Ei,jt+1
With the boundary condition that:

E0,00=exp(−pX2AB)
FUNCTION:

MolecularIntegrals.E(i,j,t,Qx,a,b,XPA,XPB,XAB)
return val
Input:

i, input values
j, input values
t, input values
Qx, input values
a, input values
b, input values
XPA, input values
XPB, input values
XAB, input values
Output:

val, value corrosponding to the given input

Overlap
-------
The overlap integrals are solved by the following recurrence relation:

.. math::
    Si+1,j=XPASij+12p(iSi−1,j+jSi,j−1)
    Si,j+1=XPBSij+12p(iSi−1,j+jSi,j−1)

With the boundary condition that:

.. math::

    S00=πp‾‾√exp(−μX2AB)

FUNCTION:

MolecularIntegrals.Overlap(a, b, la, lb, Ax, Bx)
return Sij
Input:

a, Gaussian exponent factor
b, Gaussian exponent factor
la, angular momentum quantum number
lb, angular momentum quantum number
Ax, position along one axis
Bx, position along one axis
Output:

Sij, non-normalized overlap element in one dimension


Kinetic energy
--------------
The kinetic energy integrals are solved by the following recurrence relation:

.. math::

    Ti+1,j=XPATi,j+12p(iTi−1,j+jTi,j−1)+bp(2aSi+1,j−iSi−1,j)
    Ti,j+1=XPBTi,j+12p(iTi−1,j+jTi,j−1)+ap(2bSi,j+1−iSi,j−1)

With the boundary condition that:

.. math::
    T00=[a−2a2(X2PA+12p)]S00

FUNCTION:

Kin(a, b, Ax, Ay, Az, Bx, By, Bz, la, lb, ma, mb, na, nb, N1, N2, c1, c2)
return Tij, Sij
Input:

a, Gaussian exponent factor
b, Gaussian exponent factor
Ax, position along the x-axis
Bx, position along the x-axis
Ay, position along the y-axis
By, position along the y-axis
Az, position along the z-axis
Bz, position along the z-axis
la, angular momentum quantum number
lb, angular momentum quantum number
ma, angular momentum quantum number
mb, angular momentum quantum number
na, angular momentum quantum number
nb, angular momentum quantum number
N1, normalization constant
N2, normalization constant
c1, Gaussian prefactor
c2, Gaussian prefactor
Output:

Tij, normalized kinetic energy matrix element
Sij, normalized overlap matrix element


Electron-nuclear attraction
---------------------------

The electron-nuclear interaction integral is given as:

V000ijklmn=2πp∑ti+jEijt∑uk+lEklu∑vm+nEmnvRtuv
FUNCTION:

MolecularIntegrals.elnuc(P, p, l1, l2, m1, m2, n1, n2, N1, N2, c1, c2, Zc, Ex, Ey, Ez, R1)
return Vij
Input:

P, Gaussian product
p, exponent from Guassian product
l1, angular momentum quantum number
l2, angular momentum quantum number
m1, angular momentum quantum number
n1, angular momentum quantum number
n2, angular momentum quantum number
N1, normalization constant
N2, normalization constant
c1, Gaussian prefactor
c2, Gaussian prefactor
Zc, Nuclear charge
Ex, expansion coefficients
Ey, expansion coefficients
Ez, expansion coefficients
R1, hermite coulomb integrals
Output:

Vij, normalized electron-nuclei attraction matrix element

Electron-electron repulsion
---------------------------
The electron-electron repulsion integral is calculated as:

gabcd=∑tl1+l2Eabt∑um1+m2Eabu∑vn1+n2Eabv∑τl3+l4Ecdτ∑νm3+m4Ecdν∑ϕn3+n4Ecdϕ(−1)τ+ν+ϕ2π5/2pqp+q‾‾‾‾‾√Rt+τ,u+ν,v+ϕ(α,RPQ)
FUNCTION:

MIcython.elelrep(p, q, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, N1, N2, N3, N4, c1, c2, c3, c4, E1, E2, E3, E4, E5, E6, Rpre)
return Veeijkl
Input:

p, Gaussian exponent factor from Gaussian product
q, Gaussian exponent factor from Gaussian product
l1, angular momentum quantum number
l2, angular momentum quantum number
l3, angular momentum quantum number
l4, angular momentum quantum number
m1, angular momentum quantum number
m2, angular momentum quantum number
m3, angular momentum quantum number
m4, angular momentum quantum number
n1, angular momentum quantum number
n2, angular momentum quantum number
n3, angular momentum quantum number
n4, angular momentum quantum number
N1, normalization constant
N2, normalization constant
N3, normalization constant
N4, normalization constant
c1, Gaussian prefactor
c2, Gaussian prefactor
c3, Gaussian prefactor
c4, Gaussian prefactor
E1, expansion coefficient
E2, expansion coefficient
E3, expansion coefficient
E4, expansion coefficient
E5, expansion coefficient
E6, expansion coefficient
Rpre, hermite coulomb integral
Output:

Veeijkl, normalized electron-electron repulsion matrix element

