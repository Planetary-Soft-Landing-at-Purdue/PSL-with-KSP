import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
import matplotlib.pyplot as plt 
import time
from mpl_toolkits.mplot3d import Axes3D
from math import tan, pi, exp, cos, log, factorial

start_time = time.time()

def matExp(A,x):        # finds matrix exponential using taylor polynomials, A -> matrix, x -> scalar that A is multiplied by
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A,i) * x**i/factorial(i)
    return expMat

timeTotal = 1050                      # total time of operation (seconds)
timeSteps = int(timeTotal/psl.dt)  # number of time steps, total time divided by size of time step
dt        = psl.dt                  # time step (seconds)

x_k = cp.Variable((timeSteps, 7))                    # list of state vectors
x_0, x_N      = x_k[0, :], x_k[-1, :]                # initial stateVector, final stateVector
z_k           = x_k[:, 6]                            # fuel consumption??
p_x, p_y, p_z = x_k[:, 1], x_k[:, 2], x_k[:, 0]      # x, y, z components of the position vector

eta = cp.Variable((timeSteps, 4))  # vector with thrust components and sigma, [T_x0, T_y0, T_z0, sigma0, ... T_xN, T_yN, T_zN, sigmaN]
T_x   = eta[:, 1]
T_y   = eta[:, 2]
T_z   = eta[:, 0]
sigma = eta[:, 3]

parameters  = cp.sum_squares(x_N[0:3]-psl.q)  # distance from final landing spot (xN, yN) to desired landing spot (0, 0) 
objective   = cp.Minimize(parameters)   # objective of convex solver is to minimize this distance
constraints = [x_0      == psl.x_0,                                       # initial state vector
               x_N[3:6] == [0, 0, 0],                                     # final vessel velocity has to be 0
               x_N[6]   >= log(psl.m_0 - psl.alpha*psl.rho_2*timeTotal),  # final z_k constraint
               p_z      >= 0, # the vessel has to stay above ground
               ]

x_r = matExp( psl.A, dt )
#phi = np.trapz([np.dot( matExp( psl.A, tau*.001 ), psl.B ) for tau in range(1000)], axis=0, dx=.001)

for k in range(timeSteps):   
    z_0   = log(psl.m_0 - psl.alpha*psl.rho_2 * k*dt)
    sLow  = psl.rho_1*exp(-z_0) * (1 - (z_k[k] - z_0) + .5*(z_k[k] - z_0)**2)
    sHigh = psl.rho_2*exp(-z_0) * (1 - (z_k[k] - z_0))
    
    constraints += [sigma[k] >= cp.norm(eta[k, 0:3]),                 # thrust magnitude constraint   
                    T_z[k]   >= sigma[k] * cos(psl.theta),                                  # thrust pointing constraint
                    sigma[k] >= sLow,                                                       # lower bound on thrust magnitude
                    sigma[k] <= sHigh,                                                      # upper bound on thrust magnitude
                    cp.SOC( psl.cE @ (x_k[k]-x_N), psl.SE @ (x_k[k]-x_N) )                  # second order cone constraint (glideslope constraint)
                    ]
    if(k!=(0)):
        constraints+=[x_k[k]   == cp.matmul(x_r, x_k[k-1]) + cp.matmul(psl.x_c, eta[k-1] + psl.g),# state vector
        z_k[k]<=z_k[k-1]
        ]

print("\nSet up the constraints in %s seconds." % (time.time() - start_time))
start_time = time.time()
prob = cp.Problem(objective, constraints)

try:
    prob.solve(verbose=True, solver=cp.ECOS)
    print("\nSolved problem in %s seconds.\n" % (time.time() - start_time))

    if str(x_k.value) != "None":
        print("x_k")
        print(x_k.value[:,0:6])#[:, 0:3]
        fig = plt.figure()
        ax = Axes3D(fig, azim = -60)
        #plt.plot(np.array(range(timeSteps))*dt , x_k.value[:, 0])
        ax.plot(x_k.value[:, 1], x_k.value[:, 2], x_k.value[:, 0])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        #print(x_k.value[:, 3:6])
        plt.figure(2)
        #plt.plot(np.array(range(timeSteps))*dt , np.exp(x_k.value[:, 6]))
        plt.plot(np.array(range(timeSteps))*dt, np.multiply(np.linalg.norm(eta.value[:, 0:3], axis=1),np.exp(x_k.value[:, 6])))
        #print(np.multiply(np.linalg.norm(eta.value[-1, 0:3]),np.exp(x_k.value[-1, 6])))
        x = np.linspace(0, timeTotal, timeSteps)
        plt.plot(x,psl.rho_1*np.ones((timeSteps,1)))
        plt.plot(x,psl.rho_2*np.ones((timeSteps,1)))
        plt.ylabel("thrust magnitude")
        
        plt.figure(3)
        plt.plot(np.array(range(timeSteps))*dt , np.exp(x_k.value[:, 6]))
        plt.ylabel("fuel mass")
        
        plt.figure(4)
        plt.plot(np.array(range(timeSteps))*dt, np.linalg.norm(x_k.value[:, 3:6]-psl.q, axis=1))
        plt.ylabel("speed magnitude")
            
        plt.figure(5)
        plt.plot(np.array(range(timeSteps))*dt, np.linalg.norm(np.divide(eta.value[:, 0:3]*(1e-3*np.exp(x_k.value[:, 6])).reshape((timeSteps,1)), np.exp(x_k.value[:, 6]).reshape((timeSteps,1))*abs(9.81)), axis=1))
        plt.ylabel("gforce")
        
        
        plt.show()
        
        for i in range(timeSteps):
            if(np.linalg.norm(x_k.value[i, 1:3]-psl.q[1:3])<0.01):
                print(int(i*dt))
                break
        
        
        
        dataFile, dataText = open("dataFile.csv", 'w'), ""
        for r in range(timeSteps):
            for c in range(7):
                dataText += str(x_k.value[r][c]) + ','
            for c in range(4):
                dataText += str(eta.value[r][c]) + ','
            dataText += '\n'
        dataFile.write(dataText)
        dataFile.close()
        
    else:
        print("\nNo solution was found.")
     
except cp.error.SolverError:
    print("The ECOS solver failed.")
    