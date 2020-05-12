import numpy as np
import math

dt = 1                     # time step length (seconds)

w = [0,0,0]                 # planet angular velocity    
m_f = 1505                  # dry mass
m_0 = 2905                  # wet mass
#alpha = 4.53e-4             # fuel consumption rate kg/N
alpha = 4.53e-20            # fuel consumption rate kg/N
rho_1 = 4972                # lower bound on thrust
rho_2 = 13260               # upper bound on thrust
r_N = [0,0,0]               # final position, ECEF     
rdot_N = [0,0,0]            # final velocity, ECEF
V_max = 0
theta = 180                 # pointing constraint
gamma = 4*(math.pi)/180

x_0 = np.array([600, 0, 0, 0, 0, 0, np.log(m_0)])
g  = np.array([-3.7114, 0, 0, 0])

s_W = np.array([[(w[2]**2+w[1]**2),-1*w[1]*w[0]     ,-1*w[2]*w[0]     , 0     , 2*w[2],-2*w[1]], 
                [-1*w[1]*w[0]     ,(w[2]**2+w[0]**2),-1*w[2]*w[1]     ,-2*w[2], 0     , 2*w[0]], 
                [-1*w[2]*w[0]     ,-1*w[2]*w[1]     ,(w[1]**2+w[0]**2), 2*w[1],-2*w[0], 0     ],
               ])

A = np.concatenate(( np.zeros((3, 3)), np.identity(3) ), axis=1)
A = np.concatenate(( A, s_W, np.zeros((1, 6)) ))
A = np.concatenate(( A, np.zeros((7, 1)) ), axis=1)

B = np.concatenate(( np.zeros((3, 3)), np.identity(3), np.zeros((1, 3)) ), axis=0)
B = np.concatenate(( B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T),         axis=1)

c = np.array([[1,0,0]])/math.tan(gamma)

S = np.array([[0,1,0],
              [0,0,1]])

E = np.array([[1,0,0,0,0,0,0],
              [0,1,0,0,0,0,0],
              [0,0,1,0,0,0,0]])

F = np.array([[0,0,0,0,0,0,1]])

E_u = np.array([[1,0,0,0],
                [0,1,0,0],
                [0,0,1,0]])

E_v = np.array([[0,0,0,1,0,0,0],
                [0,0,0,0,1,0,0],
                [0,0,0,0,0,1,0]])

cE = np.matmul(c,E)
SE = np.matmul(S,E)
