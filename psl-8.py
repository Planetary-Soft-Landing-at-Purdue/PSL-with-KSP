import ecos
import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt
from scipy.sparse import csc_matrix
from numpy import linalg

m_0 = 2000                      # wet mass
w = [2.53e-5, 0, 6.62e-5]       # planet angular velocity
theta = pi / 4                  # pointing constraint
gamma = pi / 6                  # glideslope constraint
m_f = 300                       # dry mass
alpha = 2.25e-4                 # fuel consumption rate kg/N
rho_1 = 4600.0                  # lower bound on thrust
rho_2 = 14200.0                 # upper bound on thrust
tFinal = 0                      # final time found after golden searcg
tElapsed = 0                    # amount of time already elapsed

g = np.array([0, 0, 0, 0, 0, 0, 0, -3.7114, 0, 0, 0])

s_W = np.array([[w[2]**2 + w[1]**2, -w[1] * w[0]     , -w[2] * w[0]     , 0        , 2 * w[2] , -2 * w[1]],
                [-w[1] * w[0]     , w[2]**2 + w[0]**2, -w[2] * w[1]     , -2 * w[2], 0        , 2 * w[0] ],
                [-w[2] * w[0]     , -w[2] * w[1]     , w[1]**2 + w[0]**2, 2 * w[1] , -2 * w[0],         0],
                ])

A = np.concatenate((np.zeros((3, 3)), np.identity(3)), axis=1)
A = np.concatenate((A, s_W, np.zeros((1, 6))))
A = np.concatenate((A, np.zeros((7, 1))), axis=1)

B = np.concatenate((np.zeros((3, 3)), np.identity(3), np.zeros((1, 3))), axis=0)
B = np.concatenate((B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T), axis=1)

def matExp(A, x):
    # approximates the exponential of a matrix, A, multiplied by a scaler, x
    # Uses a taylor polynomial:
    #       matExp(A,x) = x + xA + (xA)^2/2! + (xA)^3/3! + ... + (xA)^N/N!
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A, i) * x**i / factorial(i)
    return expMat

def equalityConstraints(tSteps):  
    # Creates the equality constraints, Ax = b that will be used in 
    # the ecos solver. Defines state vector at every time step, along 
    # with initial and final conditions.
    
    n_k, nVals = len(rowList), len(rowList) * tSteps + 16
    A_row, A_col, A_val = np.zeros(nVals), np.zeros(nVals), np.zeros(nVals)
    A_row[0:11], A_col[0:11] = np.linspace(0, 10, 11), np.linspace(0, 10, 11)
    A_val[0:11] = np.ones(11)
    
    # Iterates through every time step (except for final time step), fills in
    # data for row, column, and value for every non-zero element in A
    
    for k in range(tSteps - 1):
        start, stop = 11 + n_k * k, 11 + n_k * (k + 1)
        
        A_row[start:stop] = [11 + a + k * 7 for a in rowList]
        A_col[start:stop] = [a + k * 11 for a in colList]
        A_val[start:stop] = valList
    
    # Defines final conditions. Final velocity components are all set to zero,
    # final vertical position also set to zero
    
    start, stop = 11 + n_k * (tSteps - 1), 14 + n_k * (tSteps - 1)
    nRow, nCol = 7 * (tSteps - 1) + 11, 11 * (tSteps - 1) + 3
    
    A_row[start: stop] = [nRow, nRow + 1, nRow + 2]
    A_col[start: stop] = [nCol, nCol + 1, nCol + 2]
    A_val[start: stop] = [1, 1, 1]

    A_row[stop] = nRow + 3
    A_col[stop] = nCol - 3
    A_val[stop] = 1
    
    # b represents the right hand side of the equality constraint, Ax=b. As with
    # the matrix E, bVect is repeated for every time step (except for the last time
    # step). The first 11 elements in b are the initial conditions of the system.
    
    b = np.concatenate((x_0, np.zeros(7 * tSteps)))

    for t in range(tSteps - 1):
        b[11 + t * 7: 18 + t * 7] = bVect
    b[11:18] = np.zeros(7)
        
    return csc_matrix((A_val, (A_row, A_col)), shape=(11 + 7 * tSteps, 11 * tSteps + 1)), b.astype('float64')

def socConstraints(tSteps, goldSearch):
    # Creates linear, second order cone, and exponential cone constraints used 
    # in the ecos solver. These constraints are of the type, Gx <=_k h. First 
    # 2*tSteps rows in G are linear inequality constraints, the last 3*tSteps rows
    # are the exponential cones constraints, and everyting in between are SOC
    # constraints
    
    # q stores a list of the dimension of each SOC constraint, 
    # in order. Needed for ecos solver.
    
    q = [4, 3] * tSteps
    q.append(3)
    
    l = 2 * tSteps + 1          # number of linear inequalities
    e = tSteps                  # number of exponential cones
    dim_k = 7                   # number of rows in SOC constraints for one time step
    nCol  = 11 * (tSteps - 1)   # starting column for last time step
    nRow  = l + dim_k * tSteps  # starting row for final SOC constraints
    eRow  = nRow + 3            # starting row for exponetial cone constraints
    
    # the matrix G is initialized as a csc matrix, h is a vector. 
    # G_row, G_col, and G_val represent the row index, column index, and scalar
    # value, respectively, of each non-zero element in G. G and h are initially
    # set to the final mass constraint.
    
    h = np.zeros(eRow + 3*tSteps)
    
    G_row, G_col, G_val = [0], [nCol + 7], [1]
    h[0] = log(m_f)
    
    for t in range(tSteps):        
        k, kRow, kCol = t + 1, l + dim_k * t, 11 * t     # kRow, kCol are starting row/column for timestep
        z_0 = log(m_0 - alpha * rho_2 * t * dt)          # z_0 is used form lower/upper thrust bounds
        
        G_row.extend([k, k])                             # upper thrust constraint
        G_col.extend([kCol + 6, kCol + 10])
        G_val.extend([-rho_2 * exp(-z_0), -1])
        h[k] = rho_2 * exp(-z_0) * (1 + z_0)
        
        G_row.extend([tSteps+k, tSteps+k])               # thrust pointing constraint
        G_col.extend([kCol + 7, kCol + 10])
        G_val.extend([1, -cos(theta)])
        
        G_row.extend([kRow, kRow+1, kRow+2, kRow+3])     # thrust magnitude constraint
        G_col.extend([kCol+10, kCol+7, kCol+8, kCol+9])
        G_val.extend([1, 1, 1, 1])
        
        if t != tSteps - 1:                              # glideslope constraint for every time step except
            G_row.extend(range(kRow+4, kRow+7))          # for the last time step
            G_col.extend(range(kCol, kCol+3))
            G_val.extend([1/tan(gamma), 1, 1])
            
            G_row.extend(range(kRow+4, kRow+7))
            G_col.extend(range(nCol, nCol+3))
            G_val.extend([-1/tan(gamma), -1, -1])
                
        G_row.extend([eRow+3*t, eRow+3*t + 1])          # lower thrust bound, represented as an exponential cone
        G_col.extend([kCol + 6, kCol + 10])
        G_val.extend([-1, 1.0 / rho_1])
        h[eRow+3*t + 2] = 1
    
    # if golden search is being conducted, then final landing distance
    # must be less than slack variable. If no golden search, then final
    # distance must be less than distance found in golden search

    if goldSearch:    
        G_row.extend([nRow, nRow + 1, nRow + 2])
        G_col.extend([11 * tSteps, nCol + 1, nCol + 2])
        G_val.extend([1, 1, 1])
    else:
        G_row.extend([nRow + 1, nRow + 2])
        G_col.extend([nCol + 1, nCol + 2])
        G_val.extend([1, 1]) 

        h[nRow] = finDist

    G_val = [-g for g in G_val]     # g matrix has to be negated for ecos solver
    
    return csc_matrix((G_val, (G_row, G_col)), shape=(eRow + 3*tSteps, 11 * tSteps + 1)), h, q, l, e

def findPath(delta_t, stateVect0, initialSearch=False):
    global rowList, colList, valList, bVect, solution, x_0, dt    
    
    x_0 = stateVect0        # x_0 is the state of the spacecraft prior to finding a path 
    dt  = delta_t           # dt is the length of time for each time step
    
    # Here we define the ODE that is used to determine the vessel's state vector 
    # at every time step. This matrix is used once per time step in the equality 
    # constraints. Three seperate lists are created that contain data on non-zero 
    # values and their respective row and column indices in E. The bVect represents 
    # the right side of Ax = b for every time step.
    
    phi = matExp(A, dt)
    psi = np.trapz([np.dot(matExp(A, tau * .002 * dt), B)
                    for tau in range(500)], axis=0, dx=.002 * dt)
    omega = np.concatenate((phi, psi), axis=1)
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)

    rowList, colList, valList = [], [], []

    for r in range(len(E)):
        for c in range(len(E[0])):
            if E[r, c] != 0:
                rowList.append(r)
                colList.append(c)
                valList.append(E[r, c])
    
    bVect = np.dot(omega, g)
    
    if initialSearch:      
        goldenSearch()    
    else:
        solution = runEcos(tFinal)

def goldenSearch():
    # Conducts a golden search to find the optimal descent time. The search
    # tries to find the descent time that will result in the landing spot
    # being as close to the origin as possible. Returns distance from origin
    # and number of time steps.
    
    global finDist, tFinal
    
    zeta = .618     # golden ratio
    tLow = 0        # starting lower bound on thrust
    tHigh = 2000    # starting upper bound on thrust

    t1 = int(tHigh - zeta * (tHigh - tLow))     # initial lower search time
    t2 = int(tLow + zeta * (tHigh - tLow))      # initial upper search time
    inf_1, f_t1 = runEcos(t1, goldSearch=True)
    inf_2, f_t2 = runEcos(t2, goldSearch=True)
    
    while (inf_1 != 0 and inf_1 != 10) or (inf_2 != 0 and inf_2 != 10) or abs(t1 - t2) > int(1/dt):
        moveRight = False

        # If t1 results in an INFEASABLE solution, or both t1 and t2 result in
        # OPTIMAL solutions, then the lower boundary, tLow, is moved to the right
        # Otherwise, the higher boundary, tHigh, is moved to the left.
        
        if (inf_1 != 0 and inf_1 != 10) or ((inf_2 == 0 or inf_2 == 10) and f_t1 > f_t2):
            moveRight = True

        if moveRight:
            tLow = t1
            t1 = t2
            t2 = int(tLow + zeta * (tHigh - tLow))

            inf_1, f_t1 = inf_2, f_t2
            inf_2, f_t2 = runEcos(t2, goldSearch=True)

        else:
            tHigh = t2
            t2 = t1
            t1 = int(tHigh - zeta * (tHigh - tLow))

            inf_2, f_t2 = inf_1, f_t1
            inf_1, f_t1 = runEcos(t1, goldSearch=True)
        
        print(t1, inf_1, t2, inf_2)
              
    tFinal, finDist = t1, f_t1 
        
def runEcos(tSteps, goldSearch=False):
    G, h, q, l, e = socConstraints(tSteps, goldSearch)      # find linear, SOC, and exponential inequality constraints  
    A_mat, b = equalityConstraints(tSteps)                  # find equality constraints
    c = np.zeros(11 * tSteps + 1)                           # initializes function
    
    if goldSearch:      
        c[-1] = 1   
        
        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                              verbose=False, abstol=.00005, feastol=.00005, reltol=.00005)
        
        finDist = linalg.norm(solution['x'][11 * (tSteps - 1) + 1:11 * (tSteps - 1) + 3])
        return solution['info']['exitFlag'], finDist
    
    else:
        c[11*(tSteps-1) + 6] = -1

        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                              verbose=False, abstol=.00005, feastol=.00005, reltol=.00005)
        
        return solution['x']  

def writeData():
    dataFile, dataText = open("dataFile.csv", 'w'), ""
    for r in range(len(solution) // 11):
        for c in range(r * 11, r * 11 + 11):
            dataText += str(solution[c]) + ','
        dataText += '\n'
    dataFile.write(dataText)
    dataFile.close()

stateVect0 = np.array([4500, 600, -120, 123, -5, 11, np.log(m_0), 0, 0, 0, 0])

startTime = time.time()

findPath(.1, stateVect0, initialSearch=True)
findPath(.1, stateVect0)
writeData()

print("This took: ", time.time() - startTime)

thrustList = []
for t in range(tFinal):
    thrustVect = solution[11*t + 7:11*t + 10]
    thrustMagn = linalg.norm(thrustVect)
    thrustDire = thrustVect / thrustMagn
    
    thrustList.extend(thrustDire)
    thrustList.append(thrustMagn)
    
for t in range(0, len(thrustList), 4):
    print(thrustList[t:t+4])    
    