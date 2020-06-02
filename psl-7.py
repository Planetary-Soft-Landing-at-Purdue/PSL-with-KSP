import planetarySoftLanding as psl 
import ecos, time, math
import numpy as np  
import scipy as sp
import scipy.sparse
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csc_matrix
from operator import add

def writeData(sol):
    dataFile, dataText = open("dataFile.csv", 'w'), ""
    for r in range(tSteps):
        for c in range(r*11, r*11 + 11):
            dataText += str(sol[c]) + ','
        dataText += '\n'
    dataFile.write(dataText)
    dataFile.close()

def matExp(A,x):        # finds matrix exponential using taylor polynomials, A -> matrix, x -> scalar that A is multiplied by
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

def equalityConstraints(tSteps):
    startTime = time.time()
   
    psi   = matExp( psl.A, psl.dt )                                                                         # matrix that accounts for rotating surface
    phi   = np.trapz([np.dot( matExp(psl.A, tau*.1), psl.B ) for tau in range(10)], axis=0, dx=.1)          # matrix that accounts for effects of control: thrust and gravity
    omega = np.concatenate((psi, phi), axis=1)                                                              # combining these two matrices
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)                                          

    rowList, colList, valList = [], [], []
    
    for r in range(len(E)):                     # finds row and column index of each non-zero value in E
        for c in range(len(E[0])):
            if E[r, c] != 0:
                rowList.append(r)
                colList.append(c)
                valList.append(E[r, c])
                
    n_k, nVals               = len(rowList), len(rowList)*tSteps + 16               # number of non-zero values per time step, number of non-zero values in A
    A_row, A_col, A_val      = np.zeros(nVals), np.zeros(nVals), np.zeros(nVals)    # initializing row, column, and value lists for A                                            
    A_row[0:11], A_col[0:11] = np.linspace(0, 10, 11), np.linspace(0, 10, 11)       # 11x11 identity matrix for initial conditions
    A_val[0:11]              = np.ones(11)

    for k in range(tSteps-1):                                                   # fills in lists for sparse A matrix
        A_row[11 + n_k*k : 11 + n_k*(k+1)] = [11 + a + k*7 for a in rowList]
        A_col[11 + n_k*k : 11 + n_k*(k+1)] = [a + k*11     for a in colList]
        A_val[11 + n_k*k : 11 + n_k*(k+1)] = valList
            
    start, stop, tN = 11 + n_k*(tSteps-1), 14 + n_k*(tSteps-1), tSteps-1
    A_row[start : stop] = [11+7*tN, 12+7*tN, 13+7*tN]                       # identity matrix for setting final velocity vector to zero
    A_col[start : stop] = [3+11*tN, 4+11*tN, 5+11*tN]
    A_val[start : stop] = [1, 1, 1]

    A_row[stop] = 14+7*tN                                                   # setting final vertical position to zero
    A_col[stop] = 11*tN
    A_val[stop] = 1
    
    b = np.concatenate((psl.x_0, np.zeros(7*tSteps)))
    bVect = np.dot(omega, psl.g)                                    # bVect is repeating component of b

    for t in range(tSteps-1):                                       # fills in A and eVec
        b[11+t*7 : 18+t*7] = bVect
    b[11:18] = np.zeros(7)                                          # temporary fix to system's z-velo changing by g after one time step

    print("Setting up equality constraints took: ", time.time()-startTime)

    return csc_matrix((A_val, (A_row, A_col)), shape=(11+7*tSteps, 11*tSteps)), b.astype('float64') # stores A matrix as scipy csc_matrix

def socConstraints(tSteps):
    startTime = time.time()

    q_k   = [3, 4, 2, 3]
    dim_k = sum(q_k)
    l     = 2 * tSteps
    
    e7, e11, e8 = np.zeros((1, 11)), np.zeros((1, 11)), np.zeros((1, 11))
    e7[0, 6]   = 1
    e8[0, 7]   = 1
    e11[0, 10] = 1
    nCol = 11*(tSteps-1) # starting column for final tStep
    
    G = np.zeros((l + dim_k*tSteps+3, 11*tSteps))
    h = np.zeros(l + dim_k*tSteps+3)
        
    for k in range(tSteps):
        kRow, kCol = l + dim_k*k, 11*k          # starting row/column for timestep
        
        z_0 = log(psl.m_0 - psl.alpha*psl.rho_2 * k*psl.dt)       # values used for thrust bounds
        A  =  e7 * (.5 * psl.rho_1*exp(-z_0))**.5
        bT = -psl.rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  psl.rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)
        
        G[k, kCol:kCol+11] = -(psl.rho_2*exp(-z_0)*e7 + e11)    # thrust upper bound
        h[k]               = psl.rho_2*exp(-z_0)*(1 + z_0)
        G[tSteps+k, kCol:kCol+11] = e8 - e11 * cos(psl.theta)   # pointing constraint, no constant values so h stays as zero
        
        G[kRow:kRow+3, kCol:kCol+11] = np.concatenate((-bT/2, bT/2, A), axis=0)                       # thrust lower bound
        h[kRow:kRow+3]               = np.array([.5*(1-c), .5*(1+c), 0])
        
        G[kRow+3, kCol:kCol+11]          = e11              # thrust magnitude
        G[kRow+4:kRow+7, kCol+7:kCol+10] = np.identity(3)
        
        G[kRow+9:kRow+12, kCol:kCol+3]  = np.identity(3)    # glideslope constraint
        G[kRow+9:kRow+12, nCol:nCol+3] -= np.identity(3)    # no glideslope constraint for final tStep (z=0), if nCol = kCol, no SOC constraint
        G[kRow+9, kCol] /= tan(psl.gamma)
        G[kRow+9, nCol] /= tan(psl.gamma)
    
    G[l + dim_k*tSteps+1:l + dim_k*tSteps+3, nCol+1:nCol+3] = np.identity(2)
    h[l + dim_k*tSteps] = 5
    
    q = []
    for k in range(tSteps):
        for a in q_k:
            q.append(a)
    q.append(3)
    
    print("Setting up SOC constraints took: ", time.time()-startTime)
    return csc_matrix(-1*G), h, q, l

def setMinFunc(tSteps):
    c = np.zeros(11*tSteps)
    c[11*(tSteps-1) + 6] = -1
        
    return c
            
tSteps = 80
G, h, q, l = socConstraints(tSteps)
A, b       = equalityConstraints(tSteps)    
c          = setMinFunc(tSteps)

solution = ecos.solve(c, G, h, {'l':l, 'q':q}, A=A, b=b, feastol=.5, reltol=.5, abstol=.5)

writeData(solution['x'])
