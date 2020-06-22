#------------------------------------------------------------------------
#
# Jordan-Wigner transformation maps fermion problem into spin problem
# 
# |1> => |alpha> and |0> => |beta >:
#
#    a_j^+ => Prod_{l=1}^{j-1}(-sigma_z[l]) * sigma_+[j]
#    a_j   => Prod_{l=1}^{j-1}(-sigma_z[l]) * sigma_-[j] 
#
#------------------------------------------------------------------------
#
# Alternatively, we can use another convention [[[ Used ]]]
#
# |0> => |alpha> and |1> => |beta >: 
#
#    a_j^+ => Prod_{l=1}^{j-1}(sigma_z[l]) * sigma_-[j]
#    a_j   => Prod_{l=1}^{j-1}(sigma_z[l]) * sigma_+[j] 
#
#------------------------------------------------------------------------
import numpy
# Identity
sigma_0 = numpy.empty((2,2),dtype=complex)
sigma_0[0,0]=1
sigma_0[0,1]=0
sigma_0[1,0]=0
sigma_0[1,1]=1
# Sigma_x
sigma_x = numpy.empty((2,2),dtype=complex)
sigma_x[0,0]=0
sigma_x[0,1]=1
sigma_x[1,0]=1
sigma_x[1,1]=0
# Sigma_y
sigma_y = numpy.empty((2,2),dtype=complex)
sigma_y[0,0]=0
sigma_y[0,1]=-1.j
sigma_y[1,0]=1.j
sigma_y[1,1]=0
# Sigma_z
sigma_z = numpy.empty((2,2),dtype=complex)
sigma_z[0,0]=1
sigma_z[0,1]=0
sigma_z[1,0]=0
sigma_z[1,1]=-1
# Sigma_+/-
sigma_p = 0.5*(sigma_x + 1.j*sigma_y)
sigma_m = 0.5*(sigma_x - 1.j*sigma_y)
# Spin_a/b
spin_a = numpy.array([1.0,0.0])
spin_b = numpy.array([0.0,1.0])

#==========================
# Second convention (used)
#==========================
#cre = sigma_m.real 
#ann = sigma_p.real 
#sgn = sigma_z.real
#idn = sigma_0.real
#nii = cre.dot(ann)
sgn  = numpy.array([[ 1., 0.],[ 0.,-1.]]) 
idn  = numpy.array([[ 1., 0.],[ 0., 1.]])
idnt = numpy.array([[ 1., 0.],[ 0.,-1.]])
cre  = numpy.array([[ 0., 0.],[ 1., 0.]])
cret = numpy.array([[ 0., 0.],[ 1., 0.]])
ann  = numpy.array([[ 0., 1.],[ 0., 0.]])
annt = numpy.array([[ 0.,-1.],[ 0., 0.]])
nii  = numpy.array([[ 0., 0.],[ 0., 1.]])
niit = numpy.array([[ 0., 0.],[ 0.,-1.]])

if __name__ == '__main__':
    # test
    print(sigma_p)
    print(sigma_m)
    print(spin_a)
    print(sigma_p.dot(spin_a))
    print(sigma_p.dot(spin_b))
    print(numpy.einsum('i,j',spin_a,spin_a).reshape(len(spin_a)*2))
    print(numpy.einsum('i,j',spin_a,spin_b).reshape(len(spin_a)*2))
    
    # Matrix representation of operators in a direct product basis
    print('[a1]_{12}')
    a1=numpy.einsum('ij,kl->ikjl',sigma_p,sigma_0).reshape(4,4)
    print(a1.real) 
    print('[a2]_{12}') 
    # This is because we choose |alpha>=|1>, so exchange 1 time
    a2=numpy.einsum('ij,kl->ikjl',-sigma_z,sigma_p).reshape(4,4)
    print(a2.real)
    print("Check {a1^+,a2^+}=0")
    print(a1.dot(a2).real)
    print(a2.dot(a1).real)
    print(a1.dot(a2)+a2.dot(a1))
    print("(-sz)(-sz)")
    print(numpy.kron(-sigma_z,-sigma_z).real)
