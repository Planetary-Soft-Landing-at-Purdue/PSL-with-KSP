import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
from math import tan, pi, exp, cos, log
from numpy import dot
from numpy import multiply as mult

mainVector = cp.Variable((50, 11)) #50 steps, [stateVector(dim=7), controlVector(dim=4)]
x_k, eta   = mainVector[:, 0:7], mainVector[:, 7:11] 
x_0, x_N, xN, yN  = x_k[0, :], x_k[-1, :], x_k[-1, 1], x_k[-1, 2]
z, x, y, z_t      = x_k[:, 0], x_k[:, 1], x_k[:, 2], x_k[:, 6]
zT, yT, xT, sigma = eta[:, 0], eta[:, 1], eta[:, 2], eta[:, 3]
thrust, velo      = eta[:, 0:3], x_k[:, 3:6]

parameters  = cp.sum_squares(x_N[0:3])
objective   = cp.Minimize(parameters)
constraints = [x_0 == psl.x_0, x_N[0] == 0, x_N[3:6] == [0, 0, 0], x_N >= log(psl.m_0 - psl.alpha*psl.rho_2*50)]

for t in range(1):
    x_g = np.sum([.001 * np.dot( np.exp(psl.A * (t-.01*t1)), psl.g ) for t1 in range(100)], axis=0)
    x_c = np.sum([.001 * np.dot( np.exp(psl.A * (1-.01*t1)), psl.B ) for t1 in range(100)], axis=0)
    x_r = np.dot( np.exp(psl.A * t), psl.x_0 )
 
    z_0   = log(psl.m_0 - psl.alpha*psl.rho_2*t)
    sHigh = psl.rho_2*np.exp(-z_0) * (1 - (z_0-z_t[t]))
    sLow  = psl.rho_1*np.exp(-z_0) * (1 - (z_0-z_t[t]) + .5*(z_0-z_t[t])**2)
                                          
    constraints += [sLow <= sigma[t], sigma[t] <= sHigh, 
                    x_k[t]      == x_r + cp.matmul(x_c, eta[t]) + x_g,
                    zT[t]**2    <= sigma * cos(45),
                    z/tan(pi/4) >= cp.norm(x**2 + y**2)
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
    
