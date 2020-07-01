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
g_0   = 9.81
R_0   = 6378135
p_s   = 1.01325e5
rho_s = 1.2250
T_s   = 288.16
R     = 287

# energy-like variable used for path-finding
def E(r,v): return mu/r-0.5*v**2

def rho(h):
    a=1
    T=T_s+a*(h)
    return rho_s*(T/T_s)**-(g_0/(a*R)+1)
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
  D = 0.5 * rho(r - R_0) * V**2 * A * 0.5 / m
  L = 0.5 * rho(r - R_0) * V**2 * A * 0.25/m

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
  A     = 10000         # surface area of vessel
  Omega = 7.2921159e-5  # Earth's rotating rate

  # initial and final conditions for the vessel
  r_0 = y_0[0]
  v_0 = 2e2
  e_0 = E(r_0, v_0)
  r_f = (2e2+R_0)
  v_f = 0
  e_f = E(r_f,v_f)

  # finds the path that the vessel takes, returns the state vector
  # at each step and the final great circle range
  e   = np.linspace(e_0, e_f, 1000)
  sol = np.array(odeint(ode_e, y_0, e, args=(Omega, sigma, m, A)))

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

# greatCircle the the angle, in radians, in between the initial and
# final positions around the great circle. y_0 is the initial state 
# vector
greatCircle = .09*pi/180
y_0         = [2e4+R_0, 0, .2*pi/180, 10*pi/180, 0, greatCircle]

# sigma is the optimally found value for sigma, sol is a list of 
# statevectors at each step
sigma  = min_sigma(greatCircle, y_0)
sol, _ = solve_for_sigma(sigma, y_0, greatCircle) 

def printPath(sol):
  # plots the vessel's path, along with the Earth's surface

  xList0, yList0 = [], []
  for s in sol:
    r, s = s[0], s[5]
    x, y = r * sin(s), r * cos(s)

    xList0.append(x)
    yList0.append(y)

  earthX, earthY = [], []
  sList = np.linspace(0, .2*pi/180, len(xList0))
  for d in range(len(sList)):
    x, y = R_0 * sin(sList[d]), R_0 * cos(sList[d])
    earthX.append(x)
    earthY.append(y)

  plt.figure(1)
  plt.plot(xList0, yList0, earthX, earthY)
  plt.savefig('flightPath.png')

printPath(sol)