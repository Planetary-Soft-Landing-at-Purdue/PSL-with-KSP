import ecos
import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial
from scipy.sparse import csc_matrix

# here we define some constants that are used throughout the program

dt = .25                     # time step length (seconds)
w = [2.53e-5, 0, 6.62e-5]   # planet angular velocity
m_f = 1505                  # dry mass
m_0 = 1905                  # wet mass
alpha = 4.53e-4             # fuel consumption rate kg/N
rho_1 = 4972                # lower bound on thrust
rho_2 = 13260               # upper bound on thrust
theta = pi / 2                # pointing constraint
gamma = pi / 3                # glideslope constraint
  
x_0 = np.array([1500, 400, -220, 0, 0, 0, np.log(m_0), 0, 0, 0, 0])
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
    (B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T),         axis=1)


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


def equalityConstraints(tSteps):
    startTime = time.time()
    
    psi = matExp(A, dt)
    phi = np.trapz([np.dot(matExp(A, tau * .002), B)
                    for tau in range(500)], axis=0, dx=.002)
    omega = np.concatenate((psi, phi), axis=1)
    E = np.concatenate((-omega, np.identity(7), np.zeros((7, 4))), axis=1)
    
    rowList, colList, valList = [], [], []

    for r in range(len(E)):
        for c in range(len(E[0])):
            if E[r, c] != 0:
                rowList.append(r)
                colList.append(c)
                valList.append(E[r, c])
    
    n_k, nVals = len(rowList), len(rowList) * tSteps + 16
    A_row, A_col, A_val = np.zeros(nVals+1), np.zeros(nVals+1), np.zeros(nVals+1)
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
    
    A_row[-1], A_col[-1] = 14 + 7 * tN, 11 + 11 * tN
    
    print("Setting up equality constraints took: ", time.time() - startTime)

    return csc_matrix((A_val, (A_row, A_col)), shape=(11 + 7 * tSteps, 1 + 11 * tSteps)), b.astype('float64')


def socConstraints(tSteps):
    startTime = time.time()

    q_k = [3, 4, 3]
    dim_k = sum(q_k)
    l = 2 * tSteps + 1
    
    e7, e8, e11 = np.zeros((1, 11)), np.zeros((1, 11)), np.zeros((1, 11))
    e7[0, 6] = 1
    e8[0, 7] = 1
    e11[0, 10] = 1
    nCol = 11 * (tSteps - 1)

    G = np.zeros((l + dim_k * tSteps + 3, 11 * tSteps + 1))
    h = np.zeros(l + dim_k * tSteps + 3)

    G[0, 11*(tSteps-1):11*tSteps] = e7
    h[0] = log(m_f)
    
    for i in range(tSteps):
        k, kRow, kCol = i + 1, l + dim_k * i, 11 * i
        
        z_0 = log(m_0 - alpha * rho_2 * k * dt)
        A = e7 * (.5 * rho_1 * exp(-z_0))**.5
        bT = -rho_1 * exp(-z_0) * (e7 + z_0 * e7) - e11
        c = rho_1 * exp(-z_0) * (1 + z_0 + .5 * z_0**2)

        G[k, kCol:kCol + 11] = -(rho_2 * exp(-z_0) * e7 + e11)
        h[k] = rho_2 * exp(-z_0) * (1 + z_0)
        G[tSteps + k, kCol:kCol + 11] = e8 - e11 * cos(theta)
        
        G[kRow:kRow + 3, kCol:kCol +
            11] = np.concatenate((-bT / 2, bT / 2, A), axis=0)
        h[kRow:kRow + 3] = np.array([.5 * (1 - c), .5 * (1 + c), 0])

        G[kRow + 3, kCol:kCol + 11] = e11
        G[kRow + 4:kRow + 7, kCol + 7:kCol + 10] = np.identity(3)

        G[kRow + 7:kRow + 10, kCol:kCol + 3] = np.identity(3)
        G[kRow + 7:kRow + 10, nCol:nCol + 3] -= np.identity(3)
        G[kRow + 7, kCol] /= tan(gamma)
        G[kRow + 7, nCol] /= tan(gamma)
        
    G[l + dim_k * tSteps + 1:l + dim_k * tSteps + 3, nCol + 1:nCol + 3] = np.identity(2)
    h[l + dim_k * tSteps] = 5

    q = []
    for k in range(tSteps):
        for a in q_k:
            q.append(a)
    q.append(3)

    print("Setting up SOC constraints took: ", time.time() - startTime)
    return csc_matrix(-1 * G), h, q, l


def setMinFunc(tSteps):
    c = np.zeros(11 * tSteps + 1)
    c[-1] = 0

    return c


def runEcos(tSteps, getAns=False):
    G, h, q, l = socConstraints(tSteps)
    A_mat, b = equalityConstraints(tSteps)
    c = setMinFunc(tSteps)

    solution = ecos.solve(
        c, G, h, {'l': l, 'q': q}, A=A_mat, b=b, feastol=.05, abstol=.05, reltol=.05)

    if getAns == False:
        return solution['info']['exitFlag'], -solution['x'][11 * (tSteps - 1) + 6]
    else:
        return solution['x']


if __name__ == "__main__":
    zeta = .618
    tLow = 0
    tHigh = 500
    exitCond = .001
    count = 0

    t1 = int(tHigh - zeta * (tHigh - tLow))
    t2 = int(tLow + zeta * (tHigh - tLow))
    inf_1, f_t1 = runEcos(t1)
    inf_2, f_t2 = runEcos(t2)
    
    while (inf_1 != 0 or inf_2 != 0 or abs(f_t1 - f_t2) >= exitCond) and count < 15:
        print(t1, t2)
        moveRight = False

        if inf_1 != 0 or (inf_2 == 0 and f_t1 > f_t2):
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
        count += 1
    
    tFinal = int(.5 * float(t1 + t2))

    print(tFinal)
    solution = runEcos(tFinal, getAns=True)
    writeData(solution, tFinal)

