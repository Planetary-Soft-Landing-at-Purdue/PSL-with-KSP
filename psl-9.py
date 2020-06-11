import ecos
import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt
from scipy.sparse import csc_matrix
from numpy import linalg, concatenate

dt = 0
zeta = .618
w = [2.53e-5, 0, 6.62e-5]       # planet angular velocity
g = np.array([0, 0, 0, 0, 0, 0, 0, -3.7114, 0, 0, 0])

s_W = np.array([[w[2]**2 + w[1]**2, -w[1] * w[0]     , -w[2] * w[0]     , 0        , 2 * w[2] , -2 * w[1]],
                [-w[1] * w[0]     , w[2]**2 + w[0]**2, -w[2] * w[1]     , -2 * w[2], 0        , 2 * w[0] ],
                [-w[2] * w[0]     , -w[2] * w[1]     , w[1]**2 + w[0]**2, 2 * w[1] , -2 * w[0],         0],
                ])

A = concatenate((np.zeros((3, 3)), np.identity(3)), axis=1)
A = concatenate((A, s_W, np.zeros((1, 6))))
A = concatenate((A, np.zeros((7, 1))), axis=1)

def newVessel(m_f0, rho_10, rho_20, alpha0, gamma0, theta0):
    global m_f, rho_1, rho_2, alpha, gamma, theta, B
    m_f, rho_1, rho_2, alpha, gamma, theta = m_f0, rho_10, rho_20, alpha0, gamma0, theta0
    
    B = concatenate((np.zeros((3, 3)), np.identity(3), np.zeros((1, 3))), axis=0)
    B = concatenate((B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T), axis=1)

def matExp(A, x):
    # approximates the exponential of a matrix, A, multiplied by a scaler, x
    # Uses a taylor polynomial:
    #       matExp(A,x) = x + xA + (xA)^2/2! + (xA)^3/3! + ... + (xA)^N/N!
    
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A, i) * x**i / factorial(i)
    return expMat

def findPath(delta_t, x_0, initialSearch=False, tWait=None, tDesc=None, tDist=None):   
    print("Finding path")
    global rowList, colList, valList, bVect, dt, omega, m_0
    
    m_0 = exp(x_0[6])           # initial wet mass   

    if delta_t != dt:           # only recalculate E matrix if dt changes
        print("Finding phi/psi")
        dt  = delta_t           # dt is the length of time for each time step
        
        phi = matExp(A, dt)
        psi = np.trapz([np.dot(matExp(A, tau * .002 * dt), B)
                        for tau in range(500)], axis=0, dx=.002 * dt)
        omega = concatenate((phi, psi), axis=1)
        E = concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)

        rowList, colList, valList = [], [], []

        for r in range(len(E)):
            for c in range(len(E[0])):
                if E[r, c] != 0:
                    rowList.append(r)
                    colList.append(c)
                    valList.append(E[r, c])
        
        bVect = np.dot(omega, g)
        
    if initialSearch:      
        return findInitialPath(x_0)                                             # returns solve time from initial path
    else: 
        x_s, eta, m_s = x_0, x_0[7:11], m_0                                     # finds the state vector if the vessel waits tWait before    
        for _ in range(int(tWait / dt)):                                        # starting its descent
            x_s = concatenate(( np.dot(omega, (x_s + g)), eta ))
        m_s -= alpha * rho_2 * (tWait * dt)                                     # wet mass of vessel after waiting tWait
        
        return runEcos(int(tDesc/dt), int(tWait/dt), x_s, m_s, tDist=tDist)     # returns thrust vector with T time steps

def findInitialPath(x_0):
    print("Finding initial path")
    tHigh, tLow = 0, 10
    
    tWait_1 = int(tHigh - zeta * (tHigh - tLow))     # initial lower search time
    tWait_2 = int(tLow + zeta * (tHigh - tLow))      # initial upper search time
    tDesc_1, f_t1 = goldenSearch(tWait_1, x_0)
    tDesc_2, f_t2 = goldenSearch(tWait_2, x_0)
    
    while abs(tWait_1 - tWait_2) > int(1/dt):
        moveRight = f_t1 > f_t2

        if moveRight:
            tLow = tWait_1
            tWait_1, tDesc_1, f_t1 = tWait_2, tDesc_2, f_t2
            
            tWait_2 = int(tLow + zeta * (tHigh - tLow))
            tDesc_2, f_t2 = goldenSearch(tWait_2, x_0)

        else:
            tHigh = tWait_2
            tWait_2, tDesc_2, f_t2 = tWait_1, tDesc_1, f_t1
            
            tWait_1 = int(tHigh - zeta * (tHigh - tLow))
            tDesc_1, f_t1 = goldenSearch(tWait_1, x_0)
            
    return tWait_1 * dt, tDesc_1 * dt, f_t1      # returns wait time, descent time, and max final landing distance
 
def goldenSearch(tWait, x_0):
    print("Running golden search")
    
    x_s, eta, m_s = x_0, x_0[7:11], m_0                         # finds the state vector if the vessel waits tWait before    
    for _ in range(tWait):                                      # starting its descent
        x_s = concatenate(( np.dot(omega, (x_s + g)), eta ))
    m_s -= alpha * rho_2 * (tWait * dt)                         # wet mass of vessel after waiting tWait
    
    tLow = 0                                                    # starting lower bound on time
    tHigh = (m_s - m_f) / (alpha * rho_2 * dt)                  # starting upper bound on time

    t1 = int(tHigh - zeta * (tHigh - tLow))                     # initial lower search time
    t2 = int(tLow + zeta * (tHigh - tLow))                      # initial upper search time
    err_1, f_t1 = runEcos(t1, tWait, x_s, m_s, goldSearch=True)
    err_2, f_t2 = runEcos(t2, tWait, x_s, m_s, goldSearch=True)
    
    while (err_1 != 0 and err_1 != 10) or (err_2 != 0 and err_2 != 10) or abs(t1 - t2) > int(1/dt):
        moveRight = (err_1 != 0 and err_1 != 10) or ((err_2 == 0 or err_2 == 10) and f_t1 > f_t2)

        if moveRight:
            tLow = t1
            t1, err_1, f_t1 = t2, err_2, f_t2 
            
            t2 = int(tLow + zeta * (tHigh - tLow))
            err_2, f_t2 = runEcos(t2, tWait, x_s, m_s, goldSearch=True)

        else:
            tHigh = t2
            t2, err_2, f_t2 = t1, err_1, f_t1
            
            t1 = int(tHigh - zeta * (tHigh - tLow))
            err_1, f_t1 = runEcos(t1, tWait, x_s, m_s, goldSearch=True)
        print(tWait, '  ', t1, t2, '  ', err_1, err_2, '  ', f_t1, f_t2, '  ')                
    return t1, f_t1 
           
def runEcos(tDesc, tWait, x_s, m_s, goldSearch=False, tDist=None):
    #print("Running ecos")   
    
    G, h, q, l, e = socConstraints(tDesc, goldSearch, m_s, tDist=tDist)     # find linear, SOC, and exponential inequality constraints  
    A_mat, b = equalityConstraints(tDesc, x_s)                              # find equality constraints
    c = np.zeros(11 * tDesc + 1)                                            # initializes function
    
    if goldSearch:      
        c[-1] = 1   
        
        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                              verbose=False, abstol=1e-4, feastol=1e-4, reltol=1e-4)
        
        finDist = linalg.norm(solution['x'][11 * (tDesc - 1) + 1:11 * (tDesc - 1) + 3])
        return solution['info']['exitFlag'], finDist
    
    else:
        c[11*(tDesc-1) + 6] = -1

        solution = ecos.solve(c, G, h, {'l': l, 'q': q, 'e': e}, A=A_mat, b=b, 
                              verbose=False, abstol=1e-4, feastol=1e-4, reltol=1e-4)
        print("The solution was: ", solution['info']['exitFlag'])
        return solution['x'] 
    
def equalityConstraints(tSteps, x_s):      
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
    
    b = np.concatenate((x_s, np.zeros(7 * tSteps)))

    for t in range(tSteps - 1):
        b[11 + t * 7: 18 + t * 7] = bVect
    b[11:18] = np.zeros(7)
        
    return csc_matrix((A_val, (A_row, A_col)), shape=(11 + 7 * tSteps, 11 * tSteps + 1)), b.astype('float64')

def socConstraints(tSteps, goldSearch, m_s, tDist=None):
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
        z_0 = log(m_s - alpha * rho_2 * t * dt)          # z_0 is used form lower/upper thrust bounds
        
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

        h[nRow] = tDist * 1.1

    G_val = [-g for g in G_val]     # g matrix has to be negated for ecos solver
    
    return csc_matrix((G_val, (G_row, G_col)), shape=(eRow + 3*tSteps, 11 * tSteps + 1)), h, q, l, e
 
newVessel(300, 4600, 14200, 2.25e-4, pi/6, pi/4)

delta_t = .5
m_0 = 2000

stateVect0 = np.array([4000, 5300, -4900, 10, -50, 55, np.log(m_0), 0, 0, 0, 0])

tWait, tDesc, tDist = findPath(delta_t, stateVect0, initialSearch=True)
print(tWait, tDesc, tDist)
 
sol = findPath(delta_t, stateVect0, tWait=tWait, tDesc=tDesc, tDist=tDist)

dataFile, dataText = open("dataFile.csv", 'w'), ""
for r in range(len(sol) // 11):
    for c in range(r * 11, r * 11 + 11):
        dataText += str(sol[c]) + ','
    dataText += '\n'
dataFile.write(dataText)
dataFile.close()
