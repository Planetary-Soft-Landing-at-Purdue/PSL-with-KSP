import planetarySoftLanding as psl 
import ecos, time, math
import numpy as np  
from math import tan, pi, exp, cos, log, factorial

startTime = time.time()

def matExp(A,x):        # finds matrix exponential using taylor polynomials, A -> matrix, x -> scalar that A is multiplied by
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

def equalityConstraints(tSteps):
    psi   = matExp( psl.A, psl.dt )                                                                         # matrix that accounts for rotating surface
    phi   = np.trapz([np.dot( matExp(psl.A, tau*.01), psl.B ) for tau in range(100)], axis=0, dx=.01)       # matrix that accounts for effects of control: thrust and gravity
    omega = np.concatenate((psi, phi), axis=1)                                                              # combining these two matrices
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)                                          

    eMat = np.zeros((11 + 7*tSteps, 11*tSteps))                 # main linear equation Ax = b, A := eMat, x := [stateVector eta], b := eVect
    eMat[0:11, 0:11] = np.identity(11)                          # identity matrix for initial conditions
    eVect = psl.x_0                                             
    bVect = np.dot(omega, psl.g)                                # bVect is repeating component of eVect

    for t in range(tSteps-1):                                   # fills in eMat and eVec
        eMat[11+t*7 : 18+t*7, 11*t : 11*(t+2)] = E
        eVect = np.concatenate((eVect, bVect), axis=0) 

    eMat[5+7*tSteps : 8+7*tSteps, 11*(tSteps-1)+3 : 11*(tSteps-1)+6] = np.identity(3)   # identity matrix for final velocity (set to zero)
    eMat[4+7*tSteps, 11*(tSteps-1)] = 1                                                 # sets final z-position to zero

    eVect = np.concatenate((eVect, np.zeros(7)), axis=0)     

    return eMat, eVect

def socConstraints(tSteps):
    G_main = np.zeros((1, 11*tSteps))
    h_main = np.zeros(1)
    
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

        z_0 = log(psl.m_0 - psl.alpha*psl.rho_2 * (k//11)*psl.dt)                               # calculated lower bound for mass loss
        A  =  e7 * (.5 * psl.rho_1*exp(-z_0))**.5
        bT = -psl.rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  psl.rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)
        
        low_k, high_k = np.zeros((3, 11*tSteps)), np.zeros((1, 11*tSteps))

        low   = np.concatenate((bT, -bT, -A), axis=0)                                           # lower thrust bound
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

        G_main = np.concatenate((G_main, magThr, conPoi, gSlope, low_k, high_k), axis=0)
        h_main = np.concatenate((h_main, np.zeros(9), low_h, high_h))
    
    finDist   = np.zeros((2, 11*tSteps))
    finDist_h = np.array([0, 2, 2])
    finDist[:, 11*(tSteps-1) + 1:11*(tSteps-1) + 3] = np.identity(2)
    
    G_main = np.concatenate((G_main, finDist), axis=0)
    h_main = np.concatenate((h_main, finDist_h))
      
    return -1*G_main, h_main 
                
tSteps = 2
G_main, h_main = socConstraints(tSteps)
eMat, eVect    = equalityConstraints(tSteps)
    
print("This took: ", str(time.time() - startTime))

G_txt = ''
for r in range(len(G_main)):
    for c in G_main[r]:
        G_txt += ' '*(8-len(str(round(c, 3)))) + str(round(c, 3))
    G_txt += '  |' + ' '*(8-len(str(round(h_main[r], 3)))) + str(round(h_main[r], 3)) +  '\n'
txtFile = open("txtFile.txt", 'w')
txtFile.write(G_txt)
