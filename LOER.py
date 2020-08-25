import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt, sin
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# constants used for the algorithm

M  = 5.9724e24
m  = 100e3
G  = 6.67430e-11
mu = G * (m + M)
g_0   = 9.80665
R_0   = 6378135
p_0   = 1.01325e5
rho_0 = 1.2250
T_0   = 288.16
R_m   = 8.3144598
M     = 0.0289644
torr_c= 101325/760
h_list   = np.array([0, 11000, 25000, 47000, 53000, 79000, 90000])
T_list   = np.array([288.16, 216.65, 216.65, 282.66, 282.66, 165.66, 165.66])
p_list   = np.array([1.01e5, 2.26e4, 1.8834e1, 8.3186e-1, 3.8903e-1, 7.9019e-3])
rho_list = np.array([1.225, 3.648e-1, 4.0639e-2, 1.4757e-3, 7.1478e-4, 2.5029e-05, 3.351623e-6])
a_list   = np.array([-6.5e-3, 0, 3e-3, 0,-4.5e-3, 0, 4e-3])
R     = 287
v_0   = 2e2
dt=0.1


# energy-like variable used for path-finding
def E(r,v): return mu/r-0.5*v**2

def floor(x, n):
    if(x>n): 
        return n
    else: 
        return x
    
def rho(h_g):
    h=(R_0/(R_0+h_g))*h_g
    if(h_g<h_list[1]):
        T_1=T_list[0]
        T=T_1+a_list[0]*(h-h_list[0])
        return rho_list[0]*(T/T_1)**-(g_0/(a_list[0]*R)+1)
    
    elif(h<h_list[2]):
        T_1=T_list[1]
        return rho_list[1]*exp(-(g_0/(R*T_1))*(h-h_list[1]))
    
    elif(h<h_list[3]):
        T_1=T_list[2]
        T=T_1+a_list[2]*(h-h_list[2])
        return rho_list[2]*(T/T_1)**-(g_0/(a_list[2]*R)+1)
    
    elif(h<h_list[4]):
        T_1=T_list[3]
        return rho_list[3]*exp(-(g_0/(R*T_1))*(h-h_list[3]))
    
    elif(h<h_list[5]):
        T_1=T_list[4]
        T=T_1+a_list[4]*(h-h_list[4])
        return rho_list[4]*(T/T_1)**-(g_0/(a_list[4]*R)+1)
    
    elif(h<=h_list[6]):
        T_1=T_list[5]
        return rho_list[5]*exp(-(g_0/(R*T_1))*(h-h_list[5]))
    
    else:
        T_1=T_list[6]
        T=T_1+a_list[6]*(h-h_list[6])
        return rho_list[6]*(T/T_1)**-(g_0/(a_list[6]*R)+1)
        
def beta_r(h_g):
    h=(R_0/(R_0+h_g))*h_g
    return 0
    
def CL(alpha,M): return 0.25
def CD(alpha,M): return 0.5

def ode_e(y, e, Omega, sigma, m, A):
  '''
    Parameters:
      y (vector):     state vector -> [r, theta, phi, gamma, psi, s]_T
      e (double):     energylike variable
      Omega (double): Earth's rotation rate
      sigma (double): bank angle
      m (double):     vessel's mass
      A (double):     surface area of vessel
    Returns:
      dyde (vector): changes in y for an incremental change in e

  '''
  # an ode function, finds the change in the state vector, y, for each
  # incremental change in e. Used to find the path that the vessel will take
  # the state vector is broken up into its components, V, D, and L are
  # the velocity, Drag acceleration, and Lift acceleration magnitudes
  r, theta, phi, gamma, psi, s = y

  V = (2 * (mu / r - e))**0.5
  D = 0.5 * rho(r - R_0) * V**2 * A * 1 / m
  L = 0.5 * rho(r - R_0) * V**2 * A * 0.5 / m
 
  dyde=[#r-dot
        (V * sin(gamma)) / (D * V),
        #theta-dot
        (V * cos(gamma) * sin(psi) / (r * cos(phi))) / (D * V),
        #phi-dot
        (V * cos(gamma) * cos(psi) / r) / (D * V),
        #V-dot
        #-D-(mu*sin(gamma)/r**2)+Omega**2*r*cos(phi)*(sin(gamma)*cos(phi)- cos(gamma)*sin(phi)*cos(psi))
        #gamma-dot
        (V**-1 * (L * cos(sigma) + (V**2 - mu / r) * (cos(gamma) / r) + 2 * Omega * cos(phi) * sin(psi)
          + Omega**2 * r * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * cos(psi) * sin(phi)))) / (D * V),
        #psi-dot
        (V**-1 * ((L * sin(sigma) / cos(gamma)) + V**2/ r * cos(gamma) * sin(psi) * tan(phi)
          -2 * Omega * V * (tan(gamma) * cos(psi) * cos(phi) - sin(phi)) 
          + Omega**2 * r / cos(gamma) * sin(psi) * sin(phi) * cos(phi))) / (D * V),
        #s-dot
        (-V * cos(gamma) / r) / (D * V)
        ]
  return dyde

def ode_t(y, Omega, sigma, m, A):
  '''
    Parameters:
      y (vector):     state vector -> [r, theta, phi, gamma, psi, s]_T
      e (double):     energylike variable
      Omega (double): Earth's rotation rate
      sigma (double): bank angle
      m (double):     vessel's mass
      A (double):     surface area of vessel
    Returns:
      dydt (vector): the new state vector after one time step

  '''
  # calculates the change in the state vector after an incremental change in time

  r, theta, phi, V, gamma, psi, s = y

  D = 0.5 * rho(r - R_0) * V**2 * A * 0.50 / m
  L = 0.5 * rho(r - R_0) * V**2 * A * 0.25 / m

  y = [#r-dot
       r     + dt*(V * sin(gamma)),
       #theta-dot
       theta + dt*(V * cos(gamma) * sin(psi) / (r * cos(phi))),
       #phi-dot
       phi   + dt*(V * cos(gamma) * cos(psi) / r),
       #V-dot
       V     + dt*(-D - (mu * sin(gamma) / r**2)
        + Omega**2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * cos(psi))),
       #gamma-dot
       gamma + dt*(V**-1 * (L * cos(sigma) + (V**2 - mu / r) * (cos(gamma) / r) + 2 * Omega * cos(phi) * sin(psi)
        + Omega**2 * r * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * cos(psi) * sin(phi)))),
       #psi-dot
       psi   + dt*(V**-1 * ((L * sin(sigma) / cos(gamma)) + V**2/ r * cos(gamma) * sin(psi) * tan(phi)
        -2 * Omega * V * (tan(gamma) * cos(psi) * cos(phi) - sin(phi)) 
        + Omega**2 * r / cos(gamma) * sin(psi) * sin(phi) * cos(phi))),
       #s-dot
       s     + dt*(-V * cos(gamma) / r)
      ]
  return y

def solve_for_sigma(sigma, y_0, greatCircle):
  '''
    Parameters:
      sigma (double):       bank angle
      y_0 (vector):         initial state vector
      greatCircle (double): range to go (in radians)
    Returns:
      sol (2D array):      list of statevectors at each steo
      sol[-1][5] (double): final range to desired location

  '''
  A     = 49         # surface area of vessel
  Omega = 7.2921159e-5  # Earth's rotating rate

  # initial and final conditions for the vessel
  r_0 = y_0[0]
  e_0 = E(r_0, v_0)
  r_f = (2e2+R_0)
  v_f = 0
  e_f = E(r_f,v_f)

  # finds the path that the vessel takes, returns the state vector
  # at each step and the final great circle range
  e   = np.linspace(e_0, e_f, 1000)
  sol = np.array(odeint(ode_e, y_0, e, args=(Omega, sigma, m, A)))
  #sol[:,0]=sol[:,0]/R_0
  #sol[:,3]=sol[:,3]/((R_0*g_0)**0.5)

  # finds the time elapsed for this descent

  return sol, sol[-1][5]

def min_sigma(greatCircle, y_0):
  '''
    Parameters:
      greatCircle (double): range (in radians) from initial to final location
      y_0 (vector):         initial state vector
    Returns:
      sigma_1 (double):     optimal bank angle

  '''

  # performs a secant scheme to find the value for sigma that gets the vessel
  # as close to the desired location as necessary. Starts at the upper bound, 90 degrees.

  sigma_0, sigma_1 = pi / 4, pi / 4 - pi / 64
  _, z_0           = solve_for_sigma(sigma_0, y_0, greatCircle)
  _, z_1           = solve_for_sigma(sigma_1, y_0, greatCircle)

  epsilon = 1e-12

  # continues to search until the df/d_sigma is close enough to zero
  while abs(z_1 * (z_1 - z_0)) > epsilon:
    sigma_2 = sigma_1 - (z_1 / (z_1 - z_0)) * (sigma_1 - sigma_0) 
    sigma_0 = sigma_1
    sigma_1 = sigma_2

    z_0    = z_1
    _, z_1 = solve_for_sigma(sigma_1, y_0, greatCircle)

  return sigma_1

def find_time_descent(y_0, y_f):
  '''
    Parameters:
      y_0 (vector): initial state vector
      y_f (vector): final state vector from optimal path
    Returns:
      t (int):      descent time in seconds
      rList (list): radial distance at each time step
      sList (list): angle around great circle at each time step
  '''
  # This function finds the time of descent of the optimal path.
  # Continues to increment y until y is the same as y_f. 

  A     = 10000         # surface area of vessel
  Omega = 7.2921159e-5  # Earth's rotating rate

  y, t         = [y_0[0], y_0[1], y_0[2], v_0, y_0[3], y_0[4], y_0[5]], 0
  rList, sList = [y[0]], [y[6]]

  while y[0] > R_0 and abs(y[6] - y_f[5]) > 1e-5:
    len(y)
    y  = ode_t(y, Omega, sigma, m, A)
    t += 1

    rList.append(y[0])
    sList.append(y[6])

  return t, rList, sList

# greatCircle the the angle, in radians, in between the initial and
# final positions around the great circle. y_0 is the initial state 
# vector
greatCircle = .4*pi/180
y_0         = [8e4+R_0, 0, .2*pi/180, 10*pi/180, 0, greatCircle]

# sigma is the optimally found value for sigma, sol is a list of 
# statevectors at each step
sigma           = min_sigma(greatCircle, y_0)
sol, _          = solve_for_sigma(sigma, y_0, greatCircle) 
t, rList, sList = find_time_descent(y_0, sol[-1])
t=t*dt

print(t)

def find_correction():
    return 0

def printPath(sol, rList, sList):
  # plots the vessel's path, along with the Earth's surface

  xList0, yList0 = [], []
  for s in sol:
    r, s = s[0], s[5]
    x, y = r * sin(s), r * cos(s)

    xList0.append(x)
    yList0.append(y)

  xList1, yList1 = [], []
  for t in range(len(rList)):
    xList1.append(rList[t] * sin(sList[t]))
    yList1.append(rList[t] * cos(sList[t]))

  earthX, earthY = [], []
  sList = np.linspace(0, 1*pi/180, len(xList0))
  for d in range(len(sList)):
    x, y = R_0 * sin(sList[d]), R_0 * cos(sList[d])
    earthX.append(x)
    earthY.append(y)

  plt.figure(1)
  plt.plot(xList0, yList0, xList1, yList1, earthX, earthY)
  plt.savefig('flightPath.png')

printPath(sol, rList, sList)
