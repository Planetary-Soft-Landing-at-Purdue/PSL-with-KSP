import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
import matplotlib.pyplot as plt 
from math import tan, pi, exp, cos, log

timeTotal = 25      # number of time steps

mainVector = cp.Variable((timeTotal, 11))            # 50 steps, [stateVector(dim=7), controlVector(dim=4)]
x_k, eta   = mainVector[:, 0:7], mainVector[:, 7:11] # stateVector(dim=7), controlVector(dim=4)
x_0, x_N      = x_k[0, :], x_k[-1, :]                # initial stateVector, final stateVector
z_t, sigma    = x_k[:, 6],  eta[:, 3]                # fuel consumption??, slack variable??
p_x, p_y, p_z = x_k[:, 1], x_k[:, 2], x_k[:, 0]      # x, y, z components of the position vector
T_x, T_y, T_z = eta[:, 1], eta[:, 2], eta[:, 0]      # x, y, z components of the thrust vector

parameters  = cp.sum_squares(x_N[1:3])
objective   = cp.Minimize(parameters)
constraints = [x_0      == psl.x_0, 
               p_z[-1]  == 0, 
               x_N[3:6] == [0, 0, 0], 
               z_t[-1]  >= log(psl.m_0 - psl.alpha*psl.rho_2*10)
               ]

for t in range(timeTotal):
    x_r = cp.matmul( np.exp(psl.A * t), x_0 )
    x_c = np.sum([.01 * np.dot( np.exp(psl.A * (1-.01*tao)), psl.B ) for tao in range(100)], axis=0)
    x_g = np.sum([.01 * np.dot( np.exp(psl.A * (t-.01*tao)), psl.g ) for tao in range(100)], axis=0)
    
    z_0   = log(psl.m_0 - psl.alpha*psl.rho_2*t)
    sLow  = psl.rho_1*np.exp(-z_0) * (1 - (z_0-z_t[t]) + .5*(z_0-z_t[t])**2)
    sHigh = psl.rho_2*np.exp(-z_0) * (1 - (z_0-z_t[t]))
    
    constraints += [sLow <= sigma[t], sigma[t] <= sHigh, 
                    sigma[t]         >= cp.norm(T_x[t]**2 + T_y[t]**2 + T_z[t]**2),
                    x_k[t]           == x_r + cp.matmul(x_c, eta[t]) + x_g,
                    T_z[t]           >= sigma[t] * cos(psl.theta),
                    p_z[t]/tan(pi/4) >= cp.norm( (p_x[t] - p_x[-1])**2 + (p_y[t] - p_y[-1])**2 )
                    ]

prob = cp.Problem(objective, constraints)

try:
    prob.solve()
    if str(mainVector.value) != "None":
        print(mainVector.value[:, 0:3])
    else:
        print("No solution was found.")
     
except cp.error.SolverError:
    print("The ECOS solver failed.")
    
