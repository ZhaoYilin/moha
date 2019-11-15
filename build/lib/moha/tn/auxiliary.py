import numpy as np
import scipy
import scipy.sparse.linalg
import scipy.sparse as sparse
from scipy.linalg import expm
from copy import deepcopy

##################################################
#  auxiliary function for time evolution method  #
##################################################
def TEO_two_sites(MPO,t_interval):
    """
    #build the two sites operator
    """
    d = MPO[0].shape[2]
    D = MPO[0].shape[1]
    O = np.zeros((d**2,d**2))
    for i in range(D):
        O += np.kron(MPO[0][0,i,:,:],MPO[-1][i,0,:,:])
    #Build the two sites time evolution operator
    TEO = np.reshape(expm(-t_interval*O),(d,d,d,d))
    return TEO

def TE_two_sites(bonds,vertices,TEO,i,N,d,chi):
    """
      -+THETA+-
       |     |
       +-TEO-+
       |     |

    """
    #coarse grain
    theta = np.einsum("ij,jsk->isk",np.diag(bonds[(i-1)%N][:]),vertices[i][:,:,:])
    theta = np.einsum("isj,jk->isk",theta,np.diag(bonds[i][:]))
    theta = np.einsum("isj,jtk->istk",theta,vertices[(i+1)%N][:,:,:])
    theta = np.einsum("istj,jk->istk",theta,np.diag(bonds[(i+1)%N][:]))

    #multiple operator with wavefunction
    theta = np.einsum("istk,stuv->iuvk",theta,TEO)
    theta = np.reshape(theta,(chi*d,d*chi))

    #svd
    X,Y,Z = np.linalg.svd(theta)

    bonds[i][0:chi] = Y[0:chi]/np.sqrt(np.sum(Y[0:chi]**2))

    X = np.reshape(X[0:chi*d,0:chi],(chi,d,chi))
    vertices[i][:,:,:] = np.tensordot(np.diag(bonds[(i-1)%N][:]**(-1)),X,axes=(1,0))

    Z= np.reshape(Z[0:chi,0:d*chi],(chi,d,chi))
    vertices[(i+1)%N][:,:,:] = np.tensordot(Z,np.diag(bonds[(i+1)%N][:]**(-1)),axes=(2,0))

    return theta
###############################################
#  auxiliary function for variational method  #
###############################################
def OL_update(A, L, B):
    """
    tensor contraction from the left hand side
    +-    +--A-
    L' =  L  |
    +-    +--B-
    """
    Temp = np.einsum("sij,ik->sjk", A, L)
    L_prime = np.einsum("sjk,skl->jl", Temp, B)
    return L_prime

def OL(MPS1, MPS2, index):
    """

    """
    ## initial the left vacuum states
    L_dummy = np.zeros((1,1))
    L_dummy[0] = 1
    L = [L_dummy]
    #build L up to the left of given index 
    for i in range(0,index,1):
        L.append(OL_update(MPS1[i], L[-1], MPS2[i]))
    return L

def OR_update(A, R, B):
    """
    tensor contraction from the right hand side
     -+     -A--+
      R' =   |  R
     -+     -B--+
    """
    Temp = np.einsum("sij,jl->sil", A, R)
    R_prime = np.einsum("sil,skl->ik", Temp, B)
    return R_prime

def OR(MPS1, MPS2, index):
    """
    """
    ## initial the right vacuum states
    R_dummy = np.zeros((1,1))
    R_dummy[-1] = 1
    R = [R_dummy]
    #build R up to the right of given index 
    for i in range(len(MPS1)-1, index, -1):
        R.append(OR_update(MPS1[i], R[-1], MPS2[i]))
    return R

def overlap(MPS1,MPS2):
    """
    Function the evaluate the expectation value on tow given MPS
    <MPS1|MPS2>
    """
    return OL(MPS1,MPS2,len(MPS1))[-1]



def EL_update(W, A, L, B):
    """
    tensor contraction from the left hand side
    +-    +--A-
    |     |  |
    L' =  L--W-
    |     |  |
    +-    +--B-
    """
    Temp = np.einsum("sij,aik->sajk", A, L)
    Temp = np.einsum("sajk,abst->tbjk", Temp, W)
    L_prime = np.einsum("tbjk,tkl->bjl", Temp, B)
    return L_prime

def EL(MPS1, MPO, MPS2, index):
    """
    """
    ## initial the left vacuum states
    L_dummy = np.zeros((MPO[0].shape[0],1,1))
    L_dummy[0] = 1
    L = [L_dummy]
    #build L up to the left of given index 
    for i in range(0,index,1):
        L.append(EL_update(MPO[i], MPS1[i], L[-1], MPS2[i]))
    return L

def ER_update(W, A, R, B):
    """
    tensor contraction from the right hand side
     -+     -A--+
      |      |  |
     -R' =  -W--R
      |      |  |
     -+     -B--+
    """
    Temp = np.einsum("sij,bjl->sbil", A, R)
    Temp = np.einsum("sbil,abst->tail", Temp, W)
    R_prime = np.einsum("tail,tkl->aik", Temp, B)
    return R_prime

def ER(MPS1, MPO, MPS2, index):
    """
    """
    ## initial the right vacuum states
    R_dummy = np.zeros((MPO[-1].shape[1],1,1))
    R_dummy[-1] = 1
    R = [R_dummy]
    #build R up to the right of given index 
    for i in range(len(MPO)-1, index, -1):
        R.append(ER_update(MPO[i], MPS1[i], R[-1], MPS2[i]))
    return R
 
def expectation(MPS1, MPO, MPS2):
    """
    Function the evaluate the expectation value of an MPO on a given MPS
    <MPS1|MPO|MPS2>
    """
    return EL(MPS1,MPO,MPS2,len(MPO))[-1]

def Energy(MPS,MPO):
    """
    Function the evaluate the energy
             <MPS|MPO|MPS>
    Energy = ---------------
               <MPS|MPS>
    """
    E = expectation(MPS,MPO,MPS)
    O = overlap(MPS,MPS)
    return np.asscalar(E/O)


class HamiltonianMultiply(sparse.linalg.LinearOperator):
    """
    Functor to evaluate the Hamiltonian matrix-vector multiply
           +--A--+
           |  |  |
    -M- =  L--W--R
     |     |  |  |
           +-   -+
    """
    def __init__(self, L, W, R):
        self.L = L
        self.W = W
        self.R = R
        self.dtype = np.dtype('d')
        self.req_shape = [W.shape[2], L.shape[1], R.shape[2]]
        self.size = self.req_shape[0]*self.req_shape[1]*self.req_shape[2]
        self.shape = [self.size, self.size]
 
    def _matvec(self, A):
        M = np.einsum("aij,sik->ajsk", self.L, np.reshape(A, self.req_shape))
        M = np.einsum("ajsk,abst->bjtk", M, self.W)
        M = np.einsum("bjtk,bkl->tjl", M, self.R)
        return M

def coarse_grain_MPO(W, X):
    """
    2-1 coarse-graining of two site MPO into one site
     |     |  |
    -R- = -W--X-
     |     |  |
    """
    return np.reshape(np.einsum("abst,bcuv->acsutv",W,X),
                      [W.shape[0], X.shape[1],
                       W.shape[2]*X.shape[2],
                       W.shape[3]*X.shape[3]])
 
def product_W(W, X):
    """
    'vertical' product of MPO W-matrices
           |
     |    -W-
    -R- =  |
     |    -X-
           |
    """
    return np.reshape(np.einsum("abst,cdtu->acbdsu", W, X), [W.shape[0]*X.shape[0],
                                                             W.shape[1]*X.shape[1],
                                                             W.shape[2],X.shape[3]])
 
 
def product_MPO(M1, M2):
    assert len(M1) == len(M2)
    Result = []
    for i in range(0, len(M1)):
        Result.append(product_W(M1[i], M2[i]))
    return Result
 
 
def coarse_grain_MPS(A,B):
    """
    2-1 coarse-graining of two-site MPS into one site
      |     |  |
     -R- = -A--B-
    """
    return np.reshape(np.einsum("sij,tjk->stik",A,B),
                      [A.shape[0]*B.shape[0], A.shape[1], B.shape[2]])
 
def fine_grain_MPS(A, dims):
    assert A.shape[0] == dims[0] * dims[1]
    Theta = np.transpose(np.reshape(A, dims + [A.shape[1], A.shape[2]]),
                         (0,2,1,3))
    M = np.reshape(Theta, (dims[0]*A.shape[1], dims[1]*A.shape[2]))
    U, S, V = np.linalg.svd(M, full_matrices=0)
    U = np.reshape(U, (dims[0], A.shape[1], -1))
    V = np.transpose(np.reshape(V, (-1, dims[1], A.shape[2])), (1,0,2))
    return U, S, V
 
def truncate_SVD(U, S, V, m):
    """
    # truncate the matrices from an SVD to at most m states
    """
    m = min(len(S), m)
    trunc = np.sum(S[m:])
    S = S[0:m]
    U = U[:,:,0:m]
    V = V[:,0:m,:]
    return U,S,V,trunc,m

def optimize_one_site(A, B, W, E, F, dir):
    """
    optimize a single site given the MPO matrix W, and tensors E,F
    """
    H = HamiltonianMultiply(E,W,F)
    E,V = sparse.linalg.eigsh(H,1,v0=A,which='SA', tol=1E-8)
    A = np.reshape(V[:,0], H.req_shape)
    if (dir == 'right'):
        M = np.reshape(A,(H.req_shape[1],H.req_shape[0]*H.req_shape[2]))
        U,S,V = np.linalg.svd(M, full_matrices=0)
        A = np.reshape(V, [H.req_shape[1],H.req_shape[0],H.req_shape[2]])
        A = np.transpose(A,(1,0,2))
        US = np.einsum("ij,jk->ik", U, np.diag(S))
        B = np.einsum("sij,jk->sik", B, US)
    elif (dir == 'left'):
        M = np.reshape(A,(H.req_shape[0]*H.req_shape[1],H.req_shape[2]))
        U,S,V = np.linalg.svd(M, full_matrices=0)
        A = np.reshape(U, H.req_shape)
        SV = np.einsum("ij,jk->ik", np.diag(S),V)
        B = np.einsum("ij,sjk->sik", SV, B)
    return E[0], A, B
 
def optimize_two_sites(A, B, W1, W2, E, F, m, dir):
    """
    two-site optimization of MPS A,B with respect to MPO W1,W2 and
    environment tensors E,F
    dir = 'left' or 'right' for a left-moving or right-moving sweep
    """
    W = coarse_grain_MPO(W1,W2)
    AA = coarse_grain_MPS(A,B)
    H = HamiltonianMultiply(E,W,F)
    E,V = sparse.linalg.eigsh(H,1,v0=AA,which='SA')
    AA = np.reshape(V[:,0], H.req_shape)
    A,S,B = fine_grain_MPS(AA, [A.shape[0], B.shape[0]])
    A,S,B,trunc,m = truncate_SVD(A,S,B,m)
    if (dir == 'left'):
        B = np.einsum("ij,sjk->sik", np.diag(S), B)
    else:
        assert dir == 'right'
        A = np.einsum("sij,jk->sik", A, np.diag(S))
    return E[0], A, B, trunc, m

#############################################
#  auxiliary function for projected method  #
#############################################
def Mu(configuration):
    d =2
    D = 3
    N = len(configuration)
    vertices = []
    vertices.append(np.zeros((d,1,D)))
    for i in range(N-2):
        vertices.append(np.zeros((d,D,D)))
    vertices.append(np.zeros((d,D,1)))
    for index,content in enumerate(configuration):
        vertices[index][content][0][0] = 1
    return vertices

def d_overlap(MPS1,MPS2,index):
    """
    Functor to evaluate the Hamiltonian matrix-vector multiply

    -M-     +--A--+
     |   =  L  |  R
            +-   -+
    """
    L = OL(MPS1, MPS2, index)[-1]
    R = OR(MPS1, MPS2, index)[-1]
    A = MPS1[index]
    M = np.einsum("ij,sik->sjk", L, A)
    M = np.einsum("sjk,kl->sjl", M, R)
    return M

def d_expectation(MPS1,MPO,MPS2,index):
    """
    Functor to evaluate the Hamiltonian matrix-vector multiply
           +--A--+
           |  |  |
    -M- =  L--W--R
     |     |  |  |
           +-   -+
    """
    L = EL(MPS1, MPO, MPS2, index)[-1]
    R = ER(MPS1, MPO, MPS2, index)[-1]
    A = MPS1[index]
    W = MPO[index]
    M = np.einsum("aij,sik->ajsk", L, A)
    M = np.einsum("ajsk,abst->bjtk", M, W)
    M = np.einsum("bjtk,bkl->tjl", M, R)
    return M



def f_mu(mu,MPO,MPS):
    """
    +-+-i-+-+    
    | |   | |    +-+-i-+-+
    O-O---O-O - E| |   | |
    | |   | |    +-+MPS+-+
    +-+MPS+-+    

    """
    Exp = expectation(mu,MPO,MPS)
    E = Energy(MPS,MPO)
    Over = overlap(mu,MPS)
    B = Exp-E*Over
    return np.asscalar(B)

def D_f_mu(mu,MPO,MPS,i):
    """
    i: the index of configuration
    j: the index of A_i 
    """
    D_exp = d_expectation(mu,MPO,MPS,i)
    E = Energy(MPS,MPO)
    D_over = d_overlap(mu,MPS,i)
    C = D_exp-E*D_over
    return C

def Linear_Equation(configurations,MPO,MPS):
    M = len(configurations)
    N = len(MPS)
    Jacobian = np.zeros([M,N],dtype=object)
    g = np.zeros([M,1])
    for i,configuration in enumerate(configurations):
        mu = Mu(configuration)
        g[i] = f_mu(mu,MPO,MPS)
        for j,A in enumerate(MPS):
            Jacobian[i][j] = D_f_mu(mu,MPO,MPS,j)
    return Jacobian,g


def Jacobian_ravel(Jacobian):
    Matrix = []
    shape = []
    for i in range(Jacobian.shape[1]):
        shape.append(Jacobian[0][i].shape)

    for i in range(Jacobian.shape[0]):
        tmp = []
        for j in range(Jacobian.shape[1]):
            tmp.extend(Jacobian[i][j].ravel())
        Matrix.append(tmp)
    return np.array(Matrix),shape

def Jacobian_fold(Jacobian_ravel,shape):
    Jacobian = []
    for i in range(Jacobian_ravel.shape[0]):
        index = 0
        for j in shape:
            A = np.reshape(Jacobian_ravel[i][index:index+np.prod(j)],j)
            Jacobian.append(np.array(A))
            index +=np.prod(j)
    return Jacobian

def MPS_ravel(MPS):
    vector = []
    shape = []
    for i in range(len(MPS)):
        shape.append(MPS[i].shape)
    for i in range(len(MPS)):
        vector.extend(MPS[i].ravel())
    MPS = np.array(vector)[:,None]
    return MPS,shape

def MPS_fold(MPS_ravel,shape):
    MPS = []
    index = 0
    for i,S in enumerate(shape):
        A = np.reshape(MPS_ravel[index:index+np.prod(S)],S)
        MPS.append(np.array(A))
        index +=np.prod(S)
    return MPS

def pmps_step(configurations,MPO,MPS):
    Jacobian,g = Linear_Equation(configurations,MPO,MPS)
    Jacobian = Jacobian_ravel(Jacobian)[0]
    Jacobian_pinv = np.linalg.pinv(Jacobian)
    step = np.dot(Jacobian_pinv,g)
    MPS,shape = MPS_ravel(MPS)
    MPS = MPS-step
    MPS = MPS_fold(MPS,shape)
    return MPS










