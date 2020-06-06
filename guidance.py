import ecos
import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt
from scipy.sparse import csc_matrix

# here we define some constants that are used throughout the program

dt = .1                         # time step length (seconds)
w = [2.53e-5, 0, 6.62e-5]       # planet angular velocity
m_f = 1505                      # dry mass
m_0 = 1905                      # wet mass
alpha = 4.53e-4                 # fuel consumption rate kg/N
rho_1 = 4972                    # lower bound on thrust
rho_2 = 13260                   # upper bound on thrust
theta = 2 * pi / 16             # pointing constraint
gamma = pi / 3                  # glideslope constraint
  
x_0 = np.array([1000, 1000, -2000, 0, 0, 0, np.log(m_0), 0, 0, 0, 0])
g = np.array([0, 0, 0, 0, 0, 0, 0, -3.7114, 0, 0, 0])

s_W = np.array([[(w[2]**2 + w[1]**2), -1 * w[1] * w[0], -1 * w[2] * w[0], 0, 2 * w[2], -2 * w[1]],
                [-1 * w[1] * w[0], (w[2]**2 + w[0]**2), -
                 1 * w[2] * w[1], -2 * w[2], 0, 2 * w[0]],
                [-1 * w[2] * w[0], -1 * w[2] * w[1],
                    (w[1]**2 + w[0]**2), 2 * w[1], -2 * w[0], 0],
                ])

A = np.concatenate((np.zeros((3, 3)), np.identity(3)), axis=1)
A = np.concatenate((A, s_W, np.zeros((1, 6))))
A = np.concatenate((A, np.zeros((7, 1))), axis=1)

B = np.concatenate(
    (np.zeros((3, 3)), np.identity(3), np.zeros((1, 3))), axis=0)
B = np.concatenate(
    (B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T), axis=1)


def writeData(sol, tSteps):
    dataFile, dataText = open("dataFile.csv", 'w'), ""
    #dataText = "x, y, z, dx, dy, dz, z(t), xThrust, yThrust, zThrust, thrustMagn, \n"
    for r in range(tSteps):
        for c in range(r * 11, r * 11 + 11):
            dataText += str(sol[c]) + ','
        dataText += '\n'
    dataFile.write(dataText)
    dataFile.close()


def matExp(A, x):
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A, i) * x**i / factorial(i)
    return expMat


phi = matExp(A, dt)
psi = np.trapz([np.dot(matExp(A, tau * .001 * dt), B)
                for tau in range(1000)], axis=0, dx=.001 * dt)
omega = np.concatenate((phi, psi), axis=1)
E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)

rowList, colList, valList = [], [], []

for r in range(len(E)):
    for c in range(len(E[0])):
        if E[r, c] != 0:
            rowList.append(r)
            colList.append(c)
            valList.append(E[r, c])


def equalityConstraints(tSteps):  
    n_k, nVals = len(rowList), len(rowList) * tSteps + 16
    A_row, A_col, A_val = np.zeros(nVals), np.zeros(nVals), np.zeros(nVals)
    A_row[0:11], A_col[0:11] = np.linspace(0, 10, 11), np.linspace(0, 10, 11)
    A_val[0:11] = np.ones(11)
    
    for k in range(tSteps - 1):
        A_row[11 + n_k * k: 11 + n_k *
              (k + 1)] = [11 + a + k * 7 for a in rowList]
        A_col[11 + n_k * k: 11 + n_k * (k + 1)] = [a + k * 11 for a in colList]
        A_val[11 + n_k * k: 11 + n_k * (k + 1)] = valList
    
    start, stop, tN = 11 + n_k * \
        (tSteps - 1), 14 + n_k * (tSteps - 1), tSteps - 1
    A_row[start: stop] = [11 + 7 * tN, 12 + 7 * tN, 13 + 7 * tN]
    A_col[start: stop] = [3 + 11 * tN, 4 + 11 * tN, 5 + 11 * tN]
    A_val[start: stop] = [1, 1, 1]

    A_row[stop] = 14 + 7 * tN
    A_col[stop] = 11 * tN
    A_val[stop] = 1
    
    b = np.concatenate((x_0, np.zeros(7 * tSteps)))
    bVect = np.dot(omega, g)

    for t in range(tSteps - 1):
        b[11 + t * 7: 18 + t * 7] = bVect
    b[11:18] = np.zeros(7)
        
    return csc_matrix((A_val, (A_row, A_col)), shape=(11 + 7 * tSteps, 11 * tSteps + 1)), b.astype('float64')


def socConstraints(tSteps):
    q_k = [3, 4, 3]
    dim_k = sum(q_k)
    l = 2 * tSteps + 1
    nCol, nRow = 11 * (tSteps - 1), l + dim_k * tSteps

    h = np.zeros(l + dim_k * tSteps + 3)

    G_row, G_col, G_val = [0], [nCol + 7], [1]
    h[0] = log(m_f)
    
    for t in range(tSteps):
        k, kRow, kCol = t + 1, l + dim_k * t, 11 * t
        
        z_0 = log(m_0 - alpha * rho_2 * t * dt)
        c   = 1 + z_0 + .5 * z_0**2
        
        G_row.extend([k, k])
        G_col.extend([kCol + 6, kCol + 10])
        G_val.extend([-rho_2 * exp(-z_0), -1])
        h[k] = rho_2 * exp(-z_0) * (1 + z_0)
        
        G_row.extend([tSteps+k, tSteps+k])
        G_col.extend([kCol + 7, kCol + 10])
        G_val.extend([1, -cos(theta)])
        
        #G_row.extend([kRow, kRow, kRow+1, kRow+1, kRow+2])    
        #G_col.extend([kCol+6, kCol+10, kCol+6, kCol+10, kCol+6])    
        #G_val.extend([ .5 * (1 + z_0),  .5 / (rho_1 * exp(-z_0)),
        #              -.5 * (1 + z_0), -.5 / (rho_1 * exp(-z_0)),
        #              1 / (sqrt(2))
        #              ])
        #h[kRow:kRow + 3] = [.5 * (1 - c), .5 * (1 + c), 0]

        G_row.extend([kRow+3, kRow+4, kRow+5, kRow+6])
        G_col.extend([kCol+10, kCol+7, kCol+8, kCol+9])
        G_val.extend([1, 1, 1, 1])
        
        if t != tSteps - 1:  
            G_row.extend(range(kRow+7, kRow+10))
            G_col.extend(range(kCol, kCol+3))
            G_val.extend([1/tan(gamma), 1, 1])
            
            G_row.extend(range(kRow+7, kRow+10))
            G_col.extend(range(nCol, nCol+3))
            G_val.extend([-1/tan(gamma), -1, -1])
        
    G_row.extend([nRow, nRow + 1, nRow + 2])
    G_col.extend([11 * tSteps, nCol + 1, nCol + 2])
    G_val.extend([1, 1, 1])

    q = []
    for k in range(tSteps):
        for a in q_k:
            q.append(a)
    q.append(3)
    
    G_val = [-g for g in G_val]
    
    return csc_matrix((G_val, (G_row, G_col)), shape=(l + dim_k * tSteps + 3, 11 * tSteps + 1)), h, q, l


def runEcos(tSteps, getAns=False):
    G, h, q, l = socConstraints(tSteps)
    A_mat, b = equalityConstraints(tSteps)
    c = np.zeros(11 * tSteps + 1)
    c[-1] = 1
    
    solution = ecos.solve(
        c, G, h, {'l': l, 'q': q}, A=A_mat, b=b, verbose=False)

    if getAns == False:
        return solution['info']['exitFlag'], -solution['x'][11 * (tSteps - 1) + 6]
    else:
        return solution


def runGoldenSearch():
    zeta = .618
    tLow = 0
    tHigh = 1500

    t1 = int(tHigh - zeta * (tHigh - tLow))
    t2 = int(tLow + zeta * (tHigh - tLow))
    inf_1, f_t1 = runEcos(t1)
    inf_2, f_t2 = runEcos(t2)
    
    while (inf_1 != 0 and inf_1 != 10) or (inf_2 != 0 and inf_2 != 10) or abs(t1 - t2) > int(1/dt):
        moveRight = False

        if (inf_1 != 0 and inf_1 != 10) or ((inf_2 == 0 or inf_2 == 10) and f_t1 > f_t2):
            moveRight = True

        if moveRight == True:
            tLow = t1
            t1 = t2
            t2 = int(tLow + zeta * (tHigh - tLow))

            inf_1, f_t1 = inf_2, f_t2
            inf_2, f_t2 = runEcos(t2)

        else:
            tHigh = t2
            t2 = t1
            t1 = int(tHigh - zeta * (tHigh - tLow))

            inf_2, f_t2 = inf_1, f_t1
            inf_1, f_t1 = runEcos(t1)
        print(t1, t2)
        
    tFinal = int(.5 * float(t1 + t2))

    solution = runEcos(tFinal, getAns=True)
    writeData(solution['x'], tFinal)
    
    print("Total time steps: ", tFinal, "\n")
    print(solution['info']['exitFlag'])

