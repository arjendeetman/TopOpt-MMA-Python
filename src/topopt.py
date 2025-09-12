# A 165 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
# MMA OPTIMIZER ADDED BY ARJEN DEETMAN, NOVEMBER 2019
from __future__ import division
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
from mmapy import mmasub,subsolv
# MAIN DRIVER
def main(nelx,nely,volfrac,penal,rmin,ft,xsolv): 
    print("Minimum compliance problem with OC")
    print("ndes: " + str(nelx) + " x " + str(nely))
    print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", penal: " + str(penal))
    print("Filter method: " + ["Sensitivity based","Density based"][ft])
    print("Optimizer: " + ["OC method","MMA"][xsolv])
    # Max and min stiffness
    Emin = 1e-9
    Emax = 1.0
    # dofs:
    ndof = 2*(nelx+1)*(nely+1)
    # Allocate design variables (as array), initialize and allocate sens.
    n = nely*nelx
    x = volfrac*np.ones(nely*nelx, dtype=float)
    xPhys = x.copy()
    dc = np.zeros((nely,nelx), dtype=float)
    # Initialize OC
    if xsolv == 0:
        xold1 = x.copy()
        g = 0 # must be initialized to use the NGuyen/Paulino OC approach
    # Initialize MMA
    elif xsolv == 1:
        m = 1 
        xmin = np.zeros((n,1))
        xmax = np.ones((n,1)) 
        xval = x[np.newaxis].T 
        xold1 = xval.copy() 
        xold2 = xval.copy() 
        low = np.ones((n,1))
        upp = np.ones((n,1))
        a0 = 1.0 
        a = np.zeros((m,1)) 
        c = 10000*np.ones((m,1))
        d = np.zeros((m,1))
        move = 0.2 
    # FE: Build the index vectors for the for coo matrix format.
    KE = lk()
    edofMat = np.zeros((nelx*nely,8),dtype=int)
    for elx in range(nelx):
        for ely in range(nely):
            el = ely+elx*nely
            n1 = (nely+1)*elx+ely
            n2 = (nely+1)*(elx+1)+ely
            edofMat[el,:] = np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
    # Construct the index pointers for the coo format
    iK = np.kron(edofMat,np.ones((8,1))).flatten()
    jK = np.kron(edofMat,np.ones((1,8))).flatten()    
    # Filter: Build (and assemble) the index+data vectors for the coo matrix format
    nfilter = int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc = 0
    for i in range(nelx):
        for j in range(nely):
            row = i*nely+j
            kk1 = int(np.maximum(i-(np.ceil(rmin)-1),0))
            kk2 = int(np.minimum(i+np.ceil(rmin),nelx))
            ll1 = int(np.maximum(j-(np.ceil(rmin)-1),0))
            ll2 = int(np.minimum(j+np.ceil(rmin),nely))
            for k in range(kk1,kk2):
                for l in range(ll1,ll2):
                    col = k*nely+l
                    fac = rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
                    iH[cc] = row
                    jH[cc] = col
                    sH[cc] = np.maximum(0.0,fac)
                    cc = cc+1
    # Finalize assembly and convert to csc format
    H = coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()    
    Hs = H.sum(1)
    # BC's and support
    dofs = np.arange(2*(nelx+1)*(nely+1))
    fixed = np.union1d(dofs[0:2*(nely+1):2], np.array([2*(nelx+1)*(nely+1)-1]))
    free = np.setdiff1d(dofs,fixed)
    # Solution and RHS vectors
    f = np.zeros((ndof,1))
    u = np.zeros((ndof,1))
    # Set load
    f[1,0] = -1
    # Set loop counter and gradient vectors 
    loop = 0
    change = 1
    dv = np.ones(nely*nelx)
    dc = np.ones(nely*nelx)
    ce = np.ones(nely*nelx)
    while (change>0.001) and (loop<2000):
        loop = loop+1
        # Setup and solve FE problem
        sK = ((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
        K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
        # Remove constrained dofs from matrix
        K = K[free,:][:,free]
        # Solve system 
        u[free,0] = spsolve(K,f[free,0])    
        # Objective and sensitivity
        ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE) * u[edofMat].reshape(nelx*nely,8) ).sum(1)
        obj = ((Emin+xPhys**penal*(Emax-Emin))*ce ).sum()
        dc[:] = (-penal*xPhys**(penal-1)*(Emax-Emin))*ce
        dv[:] = np.ones(nely*nelx)
        # Sensitivity filtering:
        if ft == 0:
            dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
        elif ft == 1:
            dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
            dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
        # Optimality criteria
        if xsolv == 0:
            xold1[:] = x
            (x[:],g) = oc(nelx,nely,x,volfrac,dc,dv,g)
        # Method of moving asymptotes
        elif xsolv == 1:
            mu0 = 1.0 # Scale factor for objective function
            mu1 = 1.0 # Scale factor for volume constraint function
            f0val = mu0*obj 
            df0dx = mu0*dc[np.newaxis].T
            fval = mu1*np.array([[xPhys.sum()/(n*volfrac)-1]])
            dfdx = mu1*(dv/(n*volfrac))[np.newaxis]
            xval = x.copy()[np.newaxis].T 
            xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp = \
                mmasub(m,n,k,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,move)
            xold2 = xold1.copy()
            xold1 = xval.copy()
            x = xmma.copy().flatten()
        # Filter design variables
        if ft == 0: xPhys[:] = x
        elif ft == 1: xPhys[:] = np.asarray(H*x[np.newaxis].T/Hs)[:,0]
        # Compute the change by the inf. norm
        change = np.linalg.norm(x.reshape(nelx*nely,1)-xold1.reshape(nelx*nely,1),np.inf)
        # Write iteration history to screen (req. Python 2.6 or newer)
        print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(loop,obj,x.sum()/n,change))
    # Plot result
    fig,ax = plt.subplots()
    im = ax.imshow(-xPhys.reshape((nelx,nely)).T, cmap='gray',\
        interpolation='none', norm=colors.Normalize(vmin=-1,vmax=0))
    plt.show()    
#element stiffness matrix
def lk():
    E = 1
    nu = 0.3
    k = np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
    KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
    [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
    [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
    [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
    [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
    [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
    [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
    [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
    return (KE)
# Optimality criterion
def oc(nelx,nely,x,volfrac,dc,dv,g):
    l1 = 0
    l2 = 1e9
    move = 0.2
    # reshape to perform vector operations
    xnew = np.zeros(nelx*nely)
    while (l2-l1)/(l1+l2)>1e-3:
        lmid = 0.5*(l2+l1)
        xnew[:] = np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
        gt = g+np.sum((dv*(xnew-x)))
        if gt > 0:
            l1 = lmid
        else:
            l2 = lmid
    return (xnew,gt)
# The real main driver    
if __name__ == "__main__":
    # Default input parameters
    nelx = 150
    nely = 50
    volfrac = 0.5
    rmin = 1.5
    penal = 3.0
    ft = 1 # ft==0 -> sens, ft==1 -> dens
    xsolv = 1 # xsolv==0 -> OC, xsolv==1 -> MMA
    import sys
    if len(sys.argv)>1: nelx    = int(sys.argv[1])
    if len(sys.argv)>2: nely    = int(sys.argv[2])
    if len(sys.argv)>3: volfrac = float(sys.argv[3])
    if len(sys.argv)>4: rmin    = float(sys.argv[4])
    if len(sys.argv)>5: penal   = float(sys.argv[5])
    if len(sys.argv)>6: ft      = int(sys.argv[6])
    if len(sys.argv)>7: xsolv   = int(sys.argv[7])
    main(nelx,nely,volfrac,penal,rmin,ft,xsolv)