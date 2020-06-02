import ecos, time
import numpy as np  
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csc_matrix

# here we define some constants that are used throughout the program

tSteps = 80                 # number of time steps
dt = .5                     # time step length (seconds)
w = [2.53e-5, 0, 6.62e-5]   # planet angular velocity    
m_f = 1505                  # dry mass
m_0 = 1905                  # wet mass
alpha = 4.53e-4             # fuel consumption rate kg/N
rho_1 = 4972                # lower bound on thrust
rho_2 = 13260               # upper bound on thrust
theta = pi/4                # pointing constraint
gamma = pi/3                # glideslope constraint

# initial state of spacecraft [position, velocity, ln(initial mass), thrust, thrust magnitude]
# position, velocity, and thrust all have three components in the z, x, and y directions, in 
# that order
# g is the gravity vector, gravitation acceleration only acts in the vertical direction

x_0 = np.array([1000, 150, -40, 0, 0, 0, np.log(m_0), 0, 0, 0, 0])      
g   = np.array([0, 0, 0, -3.7114, 0, 0, 0, 0, 0, 0, 0])

# s_W represents the rotating reference frame, A and B are used in y_k definitions

s_W = np.array([[(w[2]**2+w[1]**2),-1*w[1]*w[0]     ,-1*w[2]*w[0]     , 0     , 2*w[2],-2*w[1]], 
                [-1*w[1]*w[0]     ,(w[2]**2+w[0]**2),-1*w[2]*w[1]     ,-2*w[2], 0     , 2*w[0]], 
                [-1*w[2]*w[0]     ,-1*w[2]*w[1]     ,(w[1]**2+w[0]**2), 2*w[1],-2*w[0], 0     ],
               ])

A = np.concatenate(( np.zeros((3, 3)), np.identity(3) ), axis=1)
A = np.concatenate(( A, s_W, np.zeros((1, 6)) ))
A = np.concatenate(( A, np.zeros((7, 1)) ), axis=1)

B = np.concatenate(( np.zeros((3, 3)), np.identity(3), np.zeros((1, 3)) ), axis=0)
B = np.concatenate(( B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T),         axis=1)

def writeData(sol):
    # This function writes the solution to an .csv file. It is formatted such that each line 
    # represents an individual tine step, and each column is associated with an individual variable.
            
    dataFile, dataText = open("dataFile.csv", 'w'), ""
    dataText = "x, y, z, dx, dy, dz, z(t), xThrust, yThrust, zThrust, thrustMagn, \n"
    for r in range(tSteps):
        for c in range(r*11, r*11 + 11):
            dataText += str(sol[c]) + ','
        dataText += '\n'         
    dataFile.write(dataText)
    dataFile.close()

def matExp(A,x):        
    # approximates the exponential of a matrix, A, multiplied by a scaler, x
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
   
    # Here we define the ODE system in discrete time.
    # D*y_k = omega*(y_(k-1) + g) is reformulated as:
    # D*y_k - omega*y_(k-1) = omega*g
    
    psi   = matExp( A, dt )                                                                         
    phi   = np.trapz([np.dot( matExp(A, tau*.1), B ) for tau in range(10)], axis=0, dx=.1)          
    omega = np.concatenate((psi, phi), axis=1)                                                              
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)                                          
    
    # The matrix E is used repeatedly in matrix A. Three seperate lists are created that contain 
    # data on non-zero values and their respective row and column indices in E
    
    rowList, colList, valList = [], [], []
    
    for r in range(len(E)):                     
        for c in range(len(E[0])):
            if E[r, c] != 0:
                rowList.append(r)
                colList.append(c)
                valList.append(E[r, c])
      
    # n_k represents the number of non-zero values per time step, 
    # nVals represents number of non-zero values in A     
    # initializing row, column, and value lists for A      
    # An 11x11 identity matrix is created for defining initial conditions
    
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
    # final vertical position also set to zero 
           
    start, stop, tN = 11 + n_k*(tSteps-1), 14 + n_k*(tSteps-1), tSteps-1
    A_row[start : stop] = [11+7*tN, 12+7*tN, 13+7*tN]                       
    A_col[start : stop] = [3+11*tN, 4+11*tN, 5+11*tN]
    A_val[start : stop] = [1, 1, 1]

    A_row[stop] = 14+7*tN                                                   
    A_col[stop] = 11*tN
    A_val[stop] = 1
    
    # b represents the right hand side of the equality constraint, Ax=b. As with
    # the matrix E, bVect is repeated for every time step (except for the last time
    # step). The first 11 elements in b are the initial conditions of the system.
    
    b = np.concatenate((x_0, np.zeros(7*tSteps)))
    bVect = np.dot(omega, g)                                    

    for t in range(tSteps-1):                                       
        b[11+t*7 : 18+t*7] = bVect
    b[11:18] = np.zeros(7)                                          

    print("Setting up equality constraints took: ", time.time()-startTime)

    # returns A matrix as scipy csc_matrix, returns b as 'float64' so that ecos
    # could use it
    
    return csc_matrix((A_val, (A_row, A_col)), shape=(11+7*tSteps, 11*tSteps)), b.astype('float64') 

def socConstraints(tSteps):
    # Creates linear and second order cone constraints used in the ecos solver.
    # These constraints are of the type, Gx <=_k h. First 2*tSteps rows in G are
    # linear inequality constraints, every row after the first 2*tSteps are SOC 
    # constraints.
    
    startTime = time.time()
    
    # constraints are ordered as such:
    #   rows(0 -> tSteps-1): all thrust upper bound linear constraints
    #   rows(tSteps -> 2*tSteps-1): all thrust pointing constraints
    #   rows(2*tSteps -> 12*tSteps-1): SOC constraints:
    #           - thrust lower bound
    #           - thrust magnitude constraint
    #           - glideslope constraint
    #   rows(12*tSteps-1 -> 12*tSteps+2): final landing location constraints
    
    q_k   = [3, 4, 3]
    dim_k = sum(q_k)
    l     = 2 * tSteps
    
    # e7, e9, e11 are each vectors used to isolate individual components of y_k
    # nCol is the starting column index for the final time step
    # G and h are initialized to zero matrices of the respective correct sizes,
    # and then filled in with the correct values.

    e7, e8, e11 = np.zeros((1, 11)), np.zeros((1, 11)), np.zeros((1, 11))
    e7[0, 6]   = 1
    e8[0, 7]   = 1
    e11[0, 10] = 1
    nCol = 11*(tSteps-1) 
    
    G = np.zeros((l + dim_k*tSteps+3, 11*tSteps))
    h = np.zeros(l + dim_k*tSteps+3)
        
    for k in range(tSteps):
        # kRow, kCol are starting row/column for timestep
        # z_0, A, bT, and c are all values and vectors used for 
        # the upper and lower thrust bounds
        
        kRow, kCol = l + dim_k*k, 11*k         
        z_0 = log(m_0 - alpha*rho_2 * k*dt)       
        A  =  e7 * (.5 * rho_1*exp(-z_0))**.5
        bT = -rho_1*exp(-z_0) * (e7 + z_0*e7) - e11
        c  =  rho_1*exp(-z_0) * (1 + z_0 + .5*z_0**2)
        
        # initializing thrust upper bounds and pointing constraints.
        # For pointing constraint, there are no constant values so h 
        # stays as zero
        
        G[k, kCol:kCol+11] = -(rho_2*exp(-z_0)*e7 + e11)    
        h[k]               = rho_2*exp(-z_0)*(1 + z_0)
        G[tSteps+k, kCol:kCol+11] = e8 - e11 * cos(theta)   
        
        # filling in G and h with lower thrust bound, thrust magnitude, and 
        # glideslope SOC constraints
        
        G[kRow:kRow+3, kCol:kCol+11] = np.concatenate((-bT/2, bT/2, A), axis=0)                       
        h[kRow:kRow+3]               = np.array([.5*(1-c), .5*(1+c), 0])
        
        G[kRow+3, kCol:kCol+11]          = e11              
        G[kRow+4:kRow+7, kCol+7:kCol+10] = np.identity(3)
        
        G[kRow+7:kRow+10, kCol:kCol+3]  = np.identity(3)    
        G[kRow+7:kRow+10, nCol:nCol+3] -= np.identity(3)    
        G[kRow+7, kCol] /= tan(gamma)
        G[kRow+7, nCol] /= tan(gamma)
    
    # final landing location constraint. Distance (magnitude of final x and y
    # positions) must be less than or equal to 5
    # q stores a list of the dimension of each SOC constraint, in order. Needed
    # for ecos solver. 
    
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
    # creates function for ecos to minimize. Minimizes 
    # negative of final mass, results in maximization 
    # of final mass.
    
    c = np.zeros(11*tSteps)
    c[11*(tSteps-1) + 6] = -1
        
    return c
            
if __name__ == "__main__":
    G, h, q, l = socConstraints(tSteps)
    A_mat, b   = equalityConstraints(tSteps)    
    c          = setMinFunc(tSteps)

    solution = ecos.solve(c, G, h, {'l':l, 'q':q}, A=A_mat, b=b, feastol=.05, abstol=.05, reltol=.05)

    writeData(solution['x'])
