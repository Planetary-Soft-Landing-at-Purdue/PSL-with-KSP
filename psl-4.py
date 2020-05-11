import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
import matplotlib.pyplot as plt 
import time
import matplotlib.pyplot as plt
from math import tan, pi, exp, cos, log, factorial

start_time = time.time()

def matExp(A,x):
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

timeTotal = 10                      # total time of operation (seconds)
timeSteps = int(timeTotal//psl.dt)  # number of time steps, total time divided by size of time step
dt        = psl.dt                  # time step (seconds)

x_k = cp.Variable((timeSteps, 7))                    # list of state vectors
x_0, x_N      = x_k[0, :], x_k[-1, :]                # initial stateVector, final stateVector
z_k           = x_k[:, 6]                            # fuel consumption??
p_x, p_y, p_z = x_k[:, 1], x_k[:, 2], x_k[:, 0]      # x, y, z components of the position vector

eta = cp.Variable((timeSteps * 4))
T_x   = eta[1::4]
T_y   = eta[2::4]
T_z   = eta[0::4]
sigma = eta[3::4]

parameters  = cp.sum_squares(x_N[1:3])  # distance from final landing spot (xN, yN) to desired landing spot (0, 0) 
objective   = cp.Minimize(parameters)   # objective of convex solver is to minimize this distance
constraints = [x_0      == psl.x_0,                                       # initial state vector
               #x_N[0]   == 0,                                            # final z coordinate has to be 0
               x_N[3:6] == [0, 0, 0],                                     # final vessel velocity has to be 0
               x_N[6]   >= log(psl.m_0 - psl.alpha*psl.rho_2*timeTotal),  # final z_k constraint
               p_z      >= 0
               ]

for k in range(1, timeSteps):
    x_r = matExp( psl.A, k*dt )
    
    x_c = np.zeros((7, 4*timeSteps))
    for tao in range(k):
        x_c[:, 4*tao:4*tao+4] = np.dot( matExp( -psl.A, k*dt ), psl.B )
    x_c = np.dot( matExp(psl.A, k*dt ), x_c )

    x_g = np.sum([np.dot( matExp(-psl.A, (dt*tao)), psl.g ) for tao in range(k)], axis=0)
    x_g = cp.matmul( matExp(psl.A, k*dt), x_g )
    
    z_0   = log(psl.m_0 - psl.alpha*psl.rho_2 * k*dt)
    sLow  = psl.rho_1*exp(-z_0) * (1 - (z_k[k] - z_0) + .5*(z_k[k] - z_0)**2)
    sHigh = psl.rho_2*exp(-z_0) * (1 - (z_k[k] - z_0))
    
    constraints += [sigma[k] >= cp.norm(T_x[k]**2 + T_y[k]**2 + T_z[k]**2),
                    T_z[k]   >= sigma[k] * cos(psl.theta),
                    sigma[k] >= sLow,
                    sigma[k] <= sHigh,
                    x_k[k]   == cp.matmul( x_r, x_0 ) + cp.matmul( x_c, eta ) + x_g,
                    cp.SOC( psl.cE @ (x_k[k]-x_N), psl.SE @ (x_k[k]-x_N) )
                    ]

print("\nSet up the constraints in %s seconds." % (time.time() - start_time))
start_time = time.time()
prob = cp.Problem(objective, constraints)

try:
    prob.solve()
    print("\nSolved problem in %s seconds.\n" % (time.time() - start_time))

    if str(x_k.value) != "None":
        print(x_k.value[:, 0:3])
        
        dataFile, dataText = open("dataFile.txt", 'w'), ""
        for r in range(timeSteps):
            for c in range(7):
                dataText += str(x_k.value[r][c]) + ','
            dataText += '\n'
        dataFile.write(dataText)
        dataFile.close()
        
    else:
        print("\nNo solution was found.")
     
except cp.error.SolverError:
    print("The ECOS solver failed.")