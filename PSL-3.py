import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
from math import tan, pi, exp, cos, log
from numpy import dot
from numpy import multiply as mult

alpha, m_0 = 1, 5000
rho1, rho2 = 10, 20      #for thrust upper/lower bounds, 10 and 20 are just placeholders
xInitial = np.array([100, 19, -10, -21, -3, 2, log(m_0)])
g = np.array([0, 0, 0, -3.7114, 0, 0, 0])

mainVector = cp.Variable((50, 11))
x_k, eta = mainVector[:, 0:7], mainVector[:, 7:11] 
s_i, s_f, xf, yf  = x_k[0, :], x_k[-1, :], x_k[-1, 1], x_k[-1, 2]
z, x, y, z_t      = x_k[:, 0], x_k[:, 1], x_k[:, 2], x_k[:, 6]
zT, yT, xT, gamma = eta[:, 0], eta[:, 1], eta[:, 2], eta[:, 3]
thrust, uV        = eta[:, 0:3], x_k[:, 3:6]

parameters  = cp.sum_squares(s_f[0:3])
objective   = cp.Minimize(parameters)
constraints = [s_i == xInitial, s_f[0] <= .01, s_f[3:6] <= [.01, .01, .01]]

for t in range(50):
    x_g = np.sum([.001 * np.dot( np.exp(psl.A * (t-.01*t1)), g     ) for t1 in range(100)], axis=0)
    x_c = np.sum([.001 * np.dot( np.exp(psl.A * (1-.01*t1)), psl.B ) for t1 in range(100)], axis=0)
    x_r = np.dot( np.exp(psl.A * t), xInitial )
    
    z_0   = log(m_0 - alpha*rho2*t)
    tHigh = rho2*np.exp(-z_0) * (1 - (z_0-z_t[t]))
    tLow  = rho1*np.exp(-z_0) * (1 - (z_0-z_t[t]) + .5*(z_0-z_t[t])**2)

    constraints += [tLow <= gamma[t], gamma[t] <= tHigh, 
                    x_k[t] == x_g + x_r + cp.matmul(x_c, eta[t, :])]

prob = cp.Problem(objective, constraints)

try:
    prob.solve()
    if str(mainVector.value) != "None":
        print(mainVector.value[:, 0:3])
    else:
        print("No solution was found.")
     
except cp.error.SolverError:
    print("The solver failed.")
    