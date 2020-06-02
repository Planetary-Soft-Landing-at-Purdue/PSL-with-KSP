import planetarySoftLanding as psl 
import ecos, time, math
import numpy as np  
import scipy.sparse
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csc_matrix

def writeData(sol):
    # This function writes the solution to an .csv file. It is formatted such that each line represents an individual tine step, and each
    # column is associated with an individual variable:
    #     x, y, z, dx, dy, dz, z(t), xThrust, yThrust, zThrust, thrustMagn
            
    dataFile, dataText = open("dataFile.csv", 'w'), ""
    #dataText = "x, y, z, dx, dy, dz, z(t), xThrust, yThrust, zThrust, thrustMagn, \n"
    for r in range(tSteps):
        for c in range(r*11, r*11 + 11):
            dataText += str(sol[c]) + ','
        dataText += '\n'         
    dataFile.write(dataText)
    dataFile.close()

def matExp(A,x):        
    # finds the exponential of a matrix, A, multiplied by a scaler, x
    # Uses a taylor polynomial:
    #       matExp(A,x) = x + xA + (xA)^2/2! + (xA)^3/3! + ... + (xA)^N/N!  
    
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

def equalityConstraints(tSteps):
    # Creates the equality constraints, Ax = b that will be used in the ecos solver. 
    # Defines state vector at every time step, along with initial and final conditions.
    
    startTime = time.time()
   
    # Here we define the ODE system in discrete time
    #   -psi accounts for a rotating reference frame
    #   -phi accounts for control affects: thrust and gravity
    #   -omega combines these two matrices
    #   -E shows that the state vector at each time step is a function of only the state 
    #    vector at the previous time step 
    
    psi   = matExp( psl.A, psl.dt )                                                                         
    phi   = np.trapz([np.dot( matExp(psl.A, tau*.1), psl.B ) for tau in range(10)], axis=0, dx=.1)          
    omega = np.concatenate((psi, phi), axis=1)                                                              
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)                                          
    
    # The matrix E is used repeatedly in matrix A. Three seperate lists are created that contain 
    # data on non-zero values and their locations in E
    
    rowList, colList, valList = [], [], []
    
    for r in range(len(E)):                     
        for c in range(len(E[0])):
            if E[r, c] != 0:
                rowList.append(r)
                colList.append(c)
                valList.append(E[r, c])
      
    # number of non-zero values per time step, number of non-zero values in A     
    # initializing row, column, and value lists for A      
    # 11x11 identity matrix for defining initial conditions
    
    n_k, nVals               = len(rowList), len(rowList)*tSteps + 16               
    A_row, A_col, A_val      = np.zeros(nVals), np.zeros(nVals), np.zeros(nVals)                                              
    A_row[0:11], A_col[0:11] = np.linspace(0, 10, 11), np.linspace(0, 10, 11)       
    A_val[0:11]              = np.ones(11)

    # Iterates through every time step (except for final time step), fills in 
    # data for row, column, and value for every non-zero element in A
    
    for k in range(tSteps-1):                                                   
        A_row[11 + n_k*k : 11 + n_k*(k+1)] = [11 + a + k*7 for a in rowList]
        A_col[11 + n_k*k : 11 + n_k*(k+1)] = [a + k*11     for a in colList]
        A_val[11 + n_k*k : 11 + n_k*(k+1)] = valList
    
    # Defines final conditions. Final velocity components are all set to zero,
    # final vertical position also set to zero. 
           
    start, stop, tN = 11 + n_k*(tSteps-1), 14 + n_k*(tSteps-1), tSteps-1
    A_row[start : stop] = [11+7*tN, 12+7*tN, 13+7*tN]                       
    A_col[start : stop] = [3+11*tN, 4+11*tN, 5+11*tN]
    A_val[start : stop] = [1, 1, 1]

    A_row[stop] = 14+7*tN                                                   
    A_col[stop] = 11*tN
    A_val[stop] = 1
    
    # b represents the right hand side of the equality constraint, Ax=b. As with
    # the matrix E, bVect is repeated for every time step (except for the last time
    # step). The first 11 elements in b represent the initial conditions of the 
    # system.
    
    b = np.concatenate((psl.x_0, np.zeros(7*tSteps)))
    bVect = np.dot(omega, psl.g)                                    

    for t in range(tSteps-1):                                       
        b[11+t*7 : 18+t*7] = bVect
    b[11:18] = np.zeros(7)                                          

    print("Setting up equality constraints took: ", time.time()-startTime)

    # returns A matrix as scipy csc_matrix, returns b as 'float64' so that ecos
    # could use it
    
    return csc_matrix((A_val, (A_row, A_col)), shape=(11+7*tSteps, 11*tSteps)), b.astype('float64') 

def socConstraints(tSteps):
    # Creates the second order cone constraints used in the ecos solver.
    # These constraints are of the type, Gx <=_k h. First l rows in G are
    # linear inequality constraints, every row > l are SOC constraints.
    
    startTime = time.time()
    
    # Linear constraints are ordered:
    #   thrust upper bound (0 ... N)
    #   pointing constraints (0 ... N)
    #   [thrust lower bound, thrust magnitude, glideslope constraint].T (0 ... N)
    # q_k is a list of the three SOC dimensions, in order for each time step
    
    q_k   = [3, 4, 3]
    dim_k = sum(q_k)
    l     = 2 * tSteps
    
    # e7, e9, e11 are each transposed basis vectors used to isolate individual
    # components of y_k
    # nCol is the starting column index for the final time step
    
    e7, e8, e11 = np.zeros((1, 11)), np.zeros((1, 11)), np.zeros((1, 11))
    e7[0, 6]   = 1
    e8[0, 7]   = 1
    e11[0, 10] = 1
    nCol = 11*(tSteps-1) 
    
    G = np.zeros((l + dim_k*tSteps+3, 11*tSteps))
    h = np.zeros(l + dim_k*tSteps+3)
        
    for k in range(tSteps):
        # kRow, kCol are starting row/column for timestep
        
        kRow, kCol = l + dim_k*k, 11*k          
        
        # initializing thrust upper bounds and pointing constraints.
        # For pointing constraint, there are no constant values so h 
        # stays as zero
        
        G[k, kCol:kCol+11] = -(psl.rho_2*exp(-z_0)*e7 + e11)    
        h[k]               = psl.rho_2*exp(-z_0)*(1 + z_0)
        G[tSteps+k, kCol:kCol+11] = e8 - e11 * cos(psl.theta)   
        
        # values and vectors used for thrust bounds
        
        z_0 = log(psl.m_0 - psl.alpha*psl.rho_2 * k*psl.dt)       
        A  =  e7 * (.5 * psl.rho_1*exp(-z_0))**.5
        bT = -psl.rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  psl.rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)
        
        G[kRow:kRow+3, kCol:kCol+11] = np.concatenate((-bT/2, bT/2, A), axis=0)                       
        h[kRow:kRow+3]               = np.array([.5*(1-c), .5*(1+c), 0])
        
        # thrust magnitude constraints. Magnitude of thrust components must
        # be less than or equal to sigma
        
        G[kRow+3, kCol:kCol+11]          = e11              
        G[kRow+4:kRow+7, kCol+7:kCol+10] = np.identity(3)
        
        # glideslope constraint, no glideslope constraint for final tStep (z=0), 
        # if nCol = kCol, no SOC constraint
        
        G[kRow+7:kRow+10, kCol:kCol+3]  = np.identity(3)    
        G[kRow+7:kRow+10, nCol:nCol+3] -= np.identity(3)    
        G[kRow+7, kCol] /= tan(psl.gamma)
        G[kRow+7, nCol] /= tan(psl.gamma)
    
    # final landing location constraint. Distance (magnitude of final x and y
    # positions) must be less than or equal to 5
    
    G[l + dim_k*tSteps+1:l + dim_k*tSteps+3, nCol+1:nCol+3] = np.identity(2)
    h[l + dim_k*tSteps] = 5
    
    # q stores a list of the dimension of each SOC constraint, in order. Needed
    # for ecos solver.
    
    q = []
    for k in range(tSteps):
        for a in q_k:
            q.append(a)
    q.append(3)
    
    print("Setting up SOC constraints took: ", time.time()-startTime)
    return csc_matrix(-1*G), h, q, l

def setMinFunc(tSteps):
    # creates function for ecos to minimize. Minimizes negative of final mass,
    # results in maximization of final mass.
    
    c = np.zeros(11*tSteps)
    c[11*(tSteps-1) + 6] = -1
        
    return c
            
tSteps = 80
G, h, q, l = socConstraints(tSteps)
A, b       = equalityConstraints(tSteps)    
c          = setMinFunc(tSteps)

solution = ecos.solve(c, G, h, {'l':l, 'q':q}, A=A, b=b, feastol=.1, abstol=.1, reltol=.1)

writeData(solution['x'])
