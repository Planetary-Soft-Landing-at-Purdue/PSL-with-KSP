import planetarySoftLanding as psl 
import cvxpy as cp  
import numpy as np 
from math import tan, pi, exp, cos
from numpy.linalg import norm
from numpy import dot

z_0 = np.zeros((100, 1)) #min fuel use contraint or something, should be a function of #t_k
rho1, rho2 = 10, 20      #for thrust upper/lower bounds, 10 and 20 are just placeholders

xInitial = np.array([100, .5, .1, -21, -3, 2, 5500])
xFinal   = np.array([0, 0, 0, 0, 0, 0, 0])

mainVector = cp.Variable((100, 11))
sVect, eta = mainVector[:, 0:7], mainVector[:, 7:11] 
s_i, s_f, xf, yf = sVect[0, :], sVect[-1, :], sVect[-1, 1], sVect[-1, 2]
z, x, y          = sVect[:, 0], sVect[:, 1], sVect[:, 2]
dz, dx, dy       = sVect[:, 3], sVect[:, 4], sVect[:, 5]
z_t              = sVect[:, 6]
zT, yT, xT = eta[:, 0], eta[:, 1], eta[:, 2]
thrust = eta[:, 0:3]
gamma  = eta[:, 3]
#direction vector -- unit vector of velocity
uV = sVect[:, 3:6] / cp.norm(dz**2 + dx**2 + dy**2)

parameters  = cp.sum_squares(s_f[0:3])
objective   = cp.Minimize(parameters)

constraints = [s_i == xInitial, s_f[0] <= .01, s_f[3:6] <= [.01, .01, .01],          #<--- initial/final conditions
               #cp.norm(x**2 + y**2) <= z*100,                                       #<--- glideslope contraint
               cp.norm(zT**2 + yT**2 + xT**2) <= gamma,                              #<--- Thrust magnitude constraint
               cp.sum( cp.multiply(uV, thrust), axis=1 ) <= .1 * gamma,                                      #<--- Pointing contraint
               #rho2*exp(-z_0) * (1 - (z_t - z_0))                                                               
               ]

prob = cp.Problem(objective, constraints)
prob.solve()

print(mainVector.value[:, 0:3])
