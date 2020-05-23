import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
import time
from math import tan, pi, exp, cos, log, factorial

start_time = time.time()

def matExp(A,x):        # finds matrix exponential using taylor polynomials, A -> matrix, x -> scalar that A is multiplied by
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

def runProb(constraints, objective, start_time):
    print("\nSet up the constraints in %s seconds." % (time.time() - start_time))
    start_time = time.time()
    prob = cp.Problem(objective, constraints)

    try:
        prob.solve(verbose=True, solver=cp.ECOS, feastol=.005, reltol=100)
        #prob.solve(verbose=True, solver=cp.ECOS, feastol=.005, reltol=4000, abstol=50)
        print("Solved problem in %s seconds.\n" % (time.time() - start_time))

        if str(y.value) != "None":
            dataFile, dataText = open("dataFile.csv", 'w'), ""
            for r in range(tSteps):
                for c in range(r*11, r*11 + 11):
                    dataText += str(y.value[c]) + ','
                dataText += '\n'
            dataFile.write(dataText)
            dataFile.close()

        else:
            print("No solution was found.\n")
        
    except cp.error.SolverError:
        print("The solver failed.\n")

def equalityConstraints(tSteps):
    psi   = matExp( psl.A, psl.dt )                                                                         # matrix that accounts for rotating surface
    phi   = np.trapz([np.dot( matExp( psl.A, tau*.001 ), psl.B ) for tau in range(1000)], axis=0, dx=.001)  # matrix that accounts for effects of control: thrust and gravity
    omega = np.concatenate((psi, phi), axis=1)                                                              # combining these two matrices
    D = np.concatenate((np.identity(7), np.zeros((7, 4))), axis=1)                                          
    E = np.concatenate((-omega, D), axis=1)
    
    eMat = np.zeros((11 + 7*tSteps, 11*tSteps))                 # main linear equation Ax = b, A := eMat, x := [stateVector eta], b := eVect
    eMat[0:11, 0:11] = np.identity(11)
    eVect = psl.x_0
    bVect = np.dot(omega, psl.g)                                # bVect is repeating component of eVect

    for t in range(tSteps-1):                                   # fills in eMat and eVec
        eMat[11+t*7 : 18+t*7, 11*t : 11*(t+2)] = E
        eVect = np.concatenate((eVect, bVect), axis=0) 
    eMat[5+7*tSteps : 8+7*tSteps, 11*(tSteps-1)+3 : 11*(tSteps-1)+6] = np.identity(3)
    eMat[4+7*tSteps, 11*(tSteps-1)] = 1

    eVect = np.concatenate((eVect, np.zeros(7)), axis=0)     

    return eMat, eVect    

def setConstraints(y, eMat, eVect):
    constraints = [cp.matmul(eMat, y) - eVect <= 2, cp.matmul(eMat, y) - eVect >= -2]
    y_N = y[11*(tSteps-1) : 11*(tSteps-1) + 7]

    for k in range(11, 11 * tSteps, 11):
        z_0   = log(psl.m_0 - psl.alpha*psl.rho_2 * (k//11)*psl.dt)

        sLow  = psl.rho_1*exp(-z_0) * (1 - (y[k+6] - z_0) + .5*(y[k+6] - z_0)**2)
        sHigh = psl.rho_2*exp(-z_0) * (1 - (y[k+6] - z_0)) 

        constraints += [y[k+10] >= cp.norm(y[k+7:k+10]),                                         # thrust magnitude constraint   
                        y[k+7]  >= y[k+10] * cos(psl.theta),                                    # thrust pointing constraint
                        y[k+10] >= sLow,                                                        # lower bound on thrust magnitude
                        y[k+10] <= sHigh,                                                      # upper bound on thrust magnitude
                        cp.SOC( psl.cE @ (y[k:k+7]-y_N), psl.SE @ (y[k:k+7]-y_N) )                   # second order cone constraint (glideslope constraint)
                        ]
    return constraints
    
timeTotal = 120                     # total time of operation (seconds)
tSteps = int(timeTotal//psl.dt)     # number of time steps, total time divided by size of time step

y = cp.Variable(tSteps * 11)

parameters = cp.norm( y[tSteps*11+1 : tSteps*11+3] )
objective  = cp.Minimize(parameters)

eMat, eVect = equalityConstraints(tSteps)
constraints = setConstraints(y, eMat, eVect)

runProb(constraints, objective, start_time)