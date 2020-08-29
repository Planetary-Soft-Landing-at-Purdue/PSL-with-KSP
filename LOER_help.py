import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt, sin
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt


Omega = 7.2921159e-5  # Earth's rotating rate
m     = 29e3
A     = 50

#mu = 1
G  = 6.67430e-11
M  = 5.9724e24
mu = G * (M+m)
R_0 = 6378135
g_0 = 9.80665

# for air density calculations
R = 287
h_list   = np.array([0, 11000, 25000, 47000, 53000, 79000, 90000])
T_list   = np.array([288.16, 216.65, 216.65, 282.66, 282.66, 165.66, 165.66])
p_list   = np.array([1.01e5, 2.26e4, 1.8834e1, 8.3186e-1, 3.8903e-1, 7.9019e-3])
rho_list = np.array([1.225, 3.648e-1, 4.0639e-2, 1.4757e-3, 7.1478e-4, 2.5029e-05, 3.351623e-6])
a_list   = np.array([-6.5e-3, 0, 3e-3, 0,-4.5e-3, 0, 4e-3])

def rho(h_g):
    h   = (R_0 / (R_0 + h_g)) * h_g

    if(h_g < h_list[1]):
        T_1 = T_list[0]
        T   = T_1 + a_list[0] * (h - h_list[0])
        return rho_list[0] * (T / T_1) ** -(g_0 / (a_list[0] * R) + 1)
    
    elif(h<h_list[2]):
        T_1 = T_list[1]
        return rho_list[1] * exp(-(g_0 / (R * T_1)) * (h - h_list[1]))
    
    elif(h<h_list[3]):
        T_1 = T_list[2]
        T   = T_1 + a_list[2] * (h - h_list[2])
        return rho_list[2] * (T / T_1) ** -(g_0 / (a_list[2] * R) + 1)
    
    elif(h<h_list[4]):
        T_1 = T_list[3]
        return rho_list[3] * exp(-(g_0 / (R * T_1)) * (h - h_list[3]))
    
    elif(h<h_list[5]):
        T_1 = T_list[4]
        T   = T_1 + a_list[4] * (h - h_list[4])
        return rho_list[4] * (T / T_1) ** -(g_0 / (a_list[4] * R) + 1)
    
    elif(h<=h_list[6]):
        T_1 = T_list[5]
        return rho_list[5] * exp(-(g_0 / (R * T_1)) * (h - h_list[5]))
    
    else:
        T_1 = T_list[6]
        T   = T_1 + a_list[6] * (h - h_list[6])
        return rho_list[6] * (T / T_1) ** -(g_0 / (a_list[6] * R) + 1)

def beta_r(h_g): return (rho(h_g) - rho(h_g - 1)) / (1 * rho(h_g))     
 
# energy-like variable used for path-finding
def E(r, v): return (mu / r) - (0.5 * v ** 2)

def floor(x, n): return n if x > n else x

def CL(alpha, M): return 0.25

def CD(alpha, M): return 0.5

def find_dyde(y, e, Omega, sigma, m, A):
    r, theta, phi, gamma, psi, s = y

    V = sqrt(2 * (mu / r - e))

    D = 0.5 * rho(r - R_0) * V**2 * A * 2
    L = 0.5 * rho(r - R_0) * V**2 * A * 0.25
    #print(L, g_0 * m, V)
    #print('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f' % (r/R_0, theta, V, phi, gamma, psi, L))

    # scale = D * V / sqrt(R_0 / g_0)
    # dyde=[(V * sin(gamma)),
          
    #       (V * cos(gamma) * sin(psi) / (r * cos(phi))),
          
    #       (V * cos(gamma) * cos(psi) / r),
          
    #       (1 / V * (L * cos(sigma) / m + (V**2 - mu / r) * (cos(gamma) / r) + 2 * Omega * V * cos(phi) * sin(psi)   
    #         + Omega**2 * r * cos(phi) * (cos(gamma) * cos(phi) - sin(gamma) * cos(psi) * sin(phi)))), #sign of sin(gamma) is suspect
          
    #       (1 / V * ( (L * sin(sigma) / (m * cos(gamma)) ) + V**2/ r * cos(gamma) * sin(psi) * tan(phi)      
    #         - 2 * Omega * V * (tan(gamma) * cos(psi) * cos(phi) - sin(phi)) 
    #         + Omega**2 * r / cos(gamma) * sin(psi) * sin(phi) * cos(phi))), #r term is suspect
          
    #       (-V * cos(gamma) / r)
    #     ]
    ''' dyde = [r-dot, theta-dot, phi-dot, gamma-dot, psi-dot, s-dot] '''
    #scale = D * V / sqrt(R_0 / g_0)
    scale = D * V / m
    dyde = [(V * sin(gamma)) / scale,
            
            ((V * cos(gamma) * sin(psi)) / (r * cos(phi))) / scale,
            
            ((V * cos(gamma) * cos(psi)) / r) / scale,
                
            #(-D / m - (mu * sin(gamma) / (r ** 2)) + Omega**2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * cos(psi))) / scale,
            
            ((1 / V) * (L * cos(sigma) / m + (V ** 2 - mu / r) * (cos(gamma) / r) + 2 * Omega * V * cos(phi) * sin(psi)
                + Omega ** 2 * r * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * cos(psi) * sin(phi)))) / scale,
            
            ((1 / V) * (L * sin(sigma) / (cos(gamma) * m) + ((V ** 2 / r) * cos(gamma) * sin(psi) * tan(phi))
                - 2 * Omega * V * (tan(gamma) * cos(psi) * cos(phi) - sin(phi))
                + (Omega ** 2 * r / cos(gamma) * sin(psi) * sin(phi) * cos(phi)))) / scale,
            
            (-V * cos(gamma) / r) / scale
        
            ]
    
    return dyde

def find_flight_path(sigma, y_0):
    # initial and final conditions for the vessel
    r_0 = y_0[0]
    v_0 = 7500
    e_0 = E(r_0, v_0)
    r_f = R_0 + 10e3
    v_f = 250
    e_f = E(r_f, v_f)

    t_0 = 0
    t_f = 150
    print(e_0, r_0, v_0)
    print(e_f, r_f, v_f, '\n')
    
    eList, de = np.linspace(e_0, e_f, 200), (e_f - e_0) / 200
    flightPath = np.array(odeint(find_dyde, y_0, eList, args=(Omega, sigma, m, A)))
    '''flightPath = []
    
    for e in eList:
        dyde = find_dyde(y_0, e, Omega, sigma, m, A)
        y_0  = [y_0[i] + dyde[i] * de for i in range(6)]
        
        flightPath.append(y_0)
        #print(np.array(y_0).round(3))'''

    return flightPath

def optimize_sigma(y_0):
    epsilon = 1e-12

    sigma_0, sigma_1 = pi / 4, pi / 4 - pi / 64

    z_0 = find_flight_path(sigma_0, y_0)[-1][5]
    z_1 = find_flight_path(sigma_1, y_0)[-1][5]

    # continues to search until the df/d_sigma is close enough to zero
    while abs(z_1 * (z_1 - z_0)) > epsilon:
        sigma_2 = sigma_1 - (z_1 / (z_1 - z_0)) * (sigma_1 - sigma_0) 
        sigma_0 = sigma_1
        sigma_1 = sigma_2

        z_0 = z_1
        z_1 = find_flight_path(sigma_1, y_0)[-1][5]

    return sigma_1

def graph_flight_path(flightPath, figureNum=0):
    xList, yList = [], []
    for stateVect in flightPath:
        r, s = stateVect[0], stateVect[5]
        x, y = r * sin(s), r * cos(s)

        xList.append(x)
        yList.append(y)

    earthX, earthY = [], []
    sList = np.linspace(0, pi/8, len(xList))
    for d in range(len(sList)):
        x, y = R_0 * sin(sList[d]), R_0 * cos(sList[d])
        earthX.append(x)
        earthY.append(y)

    plt.figure(figureNum)
    plt.plot(xList, yList, earthX, earthY)
    plt.savefig('flighPath' + str(figureNum) + '.png')
 
if __name__ == "__main__":
    r     = ((R_0 + 121000)) # radius
    theta = 0                     # longitude
    phi   = 0                     # latitude
    V     = 7500
    gamma = -pi / 32              # flight-path angle
    psi   = 0                     # heading angle (clockwise from north)
    s     = 2300                 # great circle angle
    
    y_0 = [r, theta, phi, gamma, psi, s]    
    #sigma      = optimize_sigma(y_0)
    #flightPath = find_flight_path(sigma, y_0)
    
    #graph_flight_path(flightPath)  
    #print("Sigma: ", sigma)
    graph_flight_path(find_flight_path(pi / 16, y_0), figureNum=1)
    #graph_flight_path(find_flight_path(pi/2, y_0), figureNum=2)
    

