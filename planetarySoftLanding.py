import numpy as np

w = [0,0,0]                 # planet angular velocity    
m_f = 1505                  # dry mass
m_0 = 1905                  # wet mass
alpha = 4.53e-4             # fuel consumption rate kg/N
rho_1 = 32000               # lower bound on thrust
rho_2 = 80000               # upper bound on thrust
r_N = [0,0,0]               # final position, ECEF     
rdot_N = [0,0,0]            # final velocity, ECEF
V_max = 0
theta = 45                  # pointing constraint

x_0 = np.array([1500, 200, 500, -75, 0, 100, np.log(m_0)])
g  = np.array([0, 0, 0, -1.62, 0, 0, 0])

s_W = np.array([[(w[2]**2+w[1]**2),-1*w[1]*w[0]     ,-1*w[2]*w[0]     , 0     , 2*w[2],-2*w[1]], 
                [-1*w[1]*w[0]     ,(w[2]**2+w[0]**2),-1*w[2]*w[1]     ,-2*w[2], 0     , 2*w[0]], 
                [-1*w[2]*w[0]     ,-1*w[2]*w[1]     ,(w[1]**2+w[0]**2), 2*w[1],-2*w[0], 0     ],
               ])

A = np.concatenate(( np.zeros((3, 3)), np.identity(3) ), axis=1)
A = np.concatenate(( A, s_W, np.zeros((1, 6)) ))
A = np.concatenate(( A, np.zeros((7, 1)) ), axis=1)

B = np.concatenate(( np.zeros((3, 3)), np.identity(3), np.zeros((1, 3)) ), axis=0)
B = np.concatenate(( B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T),         axis=1)


