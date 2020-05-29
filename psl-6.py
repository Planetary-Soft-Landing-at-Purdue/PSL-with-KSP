import planetarySoftLanding as psl 
import ecos, time, math
import numpy as np  
import scipy as sp
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csr_matrix

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
    
    return csr_matrix(A), b.astype('float64')

def socConstraints(tSteps):
    G = np.zeros((0, 11*tSteps))
    h = np.zeros(0)
    q = []
    q_k = [4, 2, 3, 3, 1]

    thrustMag = np.concatenate(( np.zeros((4, 7)), np.flip( np.identity(4), axis=1 ) ), axis=1)     # isolates thrust components from y_k
    pointCons = np.zeros((2, 11))                                                                   # isolates thrust_z and sigma, used for pointing constraint
    pointCons[0, 10] = cos(psl.theta)
    pointCons[1, 7]  = 1
    
    e7, e11 = np.zeros((1, 11)), np.zeros((1, 11))      # basis vectors e_7 and e_11
    e7[0, 6], e11[0, 10] = 1, 1
    
    for k in range(tSteps):
        magThr, conPoi = np.zeros((4, 11*tSteps)), np.zeros((2, 11*tSteps))                     # magnitude/thrust constraints      
        magThr[:, 11*k:11*(k+1)] = thrustMag        
        conPoi[:, 11*k:11*(k+1)] = pointCons

        z_0 = log(psl.m_0 - psl.alpha*psl.rho_2 * k*psl.dt)                                     # calculated lower bound for mass loss
        A  =  e7 * (.5 * psl.rho_1*exp(-z_0))**.5
        bT = -psl.rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  psl.rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)

        low_k, high_k = np.zeros((3, 11*tSteps)), np.zeros((1, 11*tSteps))

        low   = np.concatenate((-bT/2, bT, A), axis=0)                                          # lower thrust bound
        low_h = np.array([.5*(1-c), .5*(1+c), 0])
        low_k[:, 11*k:11*(k+1)] = low
        
        high   = e11 + psl.rho_2*exp(-z_0)*e7                                                   # upper thrust bound
        high_h = np.array([ psl.rho_2*exp(-z_0) * (1+z_0) ])
        high_k[:, 11*k:11*(k+1)] = high

        gSlope = np.zeros((3, 11*tSteps))                                                       # glideslope constraints
        gSlope[:, 11*k:11*k + 3]                    = np.identity(3)
        gSlope[:, 11*(tSteps-1):11*(tSteps-1) + 3] += -1*np.identity(3)
        gSlope[0][11*(tSteps-1)]                   /= tan(psl.gamma)
        gSlope[0][11*k]                            /= tan(psl.gamma)

        G = np.concatenate((G, magThr, conPoi, gSlope, low_k, high_k), axis=0)
        h = np.concatenate((h, np.zeros(9), low_h, high_h))
        
        for a in q_k:
            q.append(a)
    
    finDist   = np.zeros((2, 11*tSteps))
    finDist_h = np.array([2, 0])
    finDist[:, 11*(tSteps-1) + 1:11*(tSteps-1) + 3] = np.identity(2)
    q.append(2)
    
    G = np.concatenate((G, finDist), axis=0)
    h = np.concatenate((h, finDist_h))
    
    return csr_matrix(-1*G), h, q

def setMinFunc(tSteps):
    c = np.zeros(11*tSteps)
    c[11*tSteps-4]
        
    return c
            
tSteps = 80
G, h, q = socConstraints(tSteps)
A, b    = equalityConstraints(tSteps)
c       = setMinFunc(tSteps)

print("Setting up the constraints took: ", str(time.time() - startTime))
solution = ecos.solve(c, G, h, {'l':0, 'q':q}, A=A, b=b)

writeData(solution['x'])