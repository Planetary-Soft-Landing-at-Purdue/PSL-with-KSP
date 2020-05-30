import planetarySoftLanding as psl 
import ecos, time, math
import numpy as np  
import scipy as sp
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csc_matrix

startTime = time.time()

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
    psi   = matExp( psl.A, psl.dt )                                                                         # matrix that accounts for rotating surface
    phi   = np.trapz([np.dot( matExp(psl.A, tau*.01), psl.B ) for tau in range(100)], axis=0, dx=.01)       # matrix that accounts for effects of control: thrust and gravity
    omega = np.concatenate((psi, phi), axis=1)                                                              # combining these two matrices
    D = np.concatenate((np.identity(7), np.zeros((7, 4))), axis=1)                                          
    E = np.concatenate((-omega, D), axis=1)                                          

    A = np.zeros((11 + 7*tSteps, 11*tSteps))                 # main linear equation Ax = b, A := A, x := [stateVector eta], b := b
    A[0:11, 0:11] = np.identity(11)                          # identity matrix for initial conditions
    b = psl.x_0                                             
    bVect = np.dot(omega, psl.g)                                # bVect is repeating component of b

    for t in range(tSteps-1):                                   # fills in A and eVec
        A[11+t*7 : 18+t*7, 11*t : 11*(t+2)] = E
        b = np.concatenate((b, bVect), axis=0) 

    A[5+7*tSteps : 8+7*tSteps, 11*(tSteps-1)+3 : 11*(tSteps-1)+6] = np.identity(3)   # identity matrix for final velocity (set to zero)
    A[4+7*tSteps, 11*(tSteps-1)] = 1                                                 # sets final z-position to zero

    b = np.concatenate((b, np.zeros(7)), axis=0)     
    
    return csc_matrix(A), b.astype('float64')

def socConstraints(tSteps):
    q_k   = [3, 2, 4, 2, 3]
    dim_k = sum(q_k)
    e7, e11, e8 = np.zeros((1, 11)), np.zeros((1, 11)), np.zeros((1, 11))
    e7[0, 6]   = 1
    e8[0, 7]   = 1
    e11[0, 10] = 1
    
    G = np.zeros((dim_k*tSteps+3, 11*tSteps))
    h = np.zeros(dim_k*tSteps+3)
    
    for k in range(tSteps-1):
        kRow, kCol = dim_k*k, 11*k
        
        z_0 = log(psl.m_0 - psl.alpha*psl.rho_2 * k*psl.dt)       # values used for thrust bounds
        A  =  e7 * (.5 * psl.rho_1*exp(-z_0))**.5
        bT = -psl.rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  psl.rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)
        
        G[kRow:kRow+3, kCol:kCol+11] = np.concatenate((-bT/2, bT, A), axis=0)                       # thrust lower bound
        h[kRow:kRow+3]               = np.array([.5*(1-c), .5*(1+c), 0])
        G[kRow+4, kCol:kCol+11] = e11 + psl.rho_2*exp(-z_0)*e7                                      # thrust upper bound
        h[kRow+3]               = psl.rho_2*exp(-z_0) * (1 + z_0)
        
        G[kRow+5, kCol:kCol+11]          = e11              # thrust magnitude
        G[kRow+6:kRow+9, kCol+7:kCol+10] = np.identity(3)
        G[kRow+10, kCol:kCol+11] = e8                       # pointing constraint
        G[kRow+11, kCol:kCol+11] = e11 * cos(psl.theta)
    
    G[dim_k*tSteps+1:dim_k*tSteps+3, 11*(tSteps-1)+1:11*(tSteps-1)+3] = np.identity(2)
    h[dim_k*tSteps] = 2
    
    q = []
    for k in range(tSteps):
        for a in q_k:
            q.append(a)
    q.append(3)
    
    return csc_matrix(-1*G), h, q

def setMinFunc(tSteps):
    c = np.zeros(11*tSteps)
    c[11*(tSteps-1) + 6] = 1
        
    return c
            
tSteps = 80
G, h, q = socConstraints(tSteps)
A, b    = equalityConstraints(tSteps)
c       = setMinFunc(tSteps)

print("Setting up the constraints took: ", str(time.time() - startTime))
solution = ecos.solve(c, G, h, {'l':0, 'q':q}, A=A, b=b)

writeData(solution['x'])

