import time
import numpy as np
from math import tan, pi, exp, cos, log, factorial, sqrt, sin
import scipy
from scipy.integrate import odeint


def E(r,v):
    return 1/r-0.5*v**2
p_s=1.01325e5
rho_s=1.2250
g_0=9.81
R_0=6378135
T_s=288.16
R=287
def rho(h):
    a=1
    T=T_s+a*(h)
    return rho_s*(T/T_s)**-(g_0/(a*R)+1)
def CL(alpha,M):
    return 0.25
def CD(alpha,M):
    return 0.5
def ode_sys(y,e, Omega, sigma, m, A):
    r, theta, phi, gamma, psi, s = y
    V=(2*(1/r-e))**0.5
    D=0.5*rho(r)*V**2*A
    L=0.5*rho(r)*V**2*A
    dyde=[#r-dot
          V*sin(gamma),
          #theta-dot
          V*cos(gamma)*sin(psi)/(r*cos(phi)),
          #phi-dot
          V*cos(gamma)*cos(psi)/r,
          #V-dot
          #-D-(sin(gamma)/r**2)+Omega**2*r*cos(phi)*(sin(gamma)*cos(phi)- cos(gamma)*sin(phi)*cos(psi))
          #gamma-dot
          V**-1*(L*cos(sigma)+(V**2-1/r)*(cos(gamma)/r)+2*Omega*cos(phi)*sin(psi)+Omega**2*r*cos(phi)*(cos(gamma)*cos(phi)+sin(gamma)*cos(psi)*sin(phi))),
          #psi-dot
          V**-1*((L*sin(sigma)/cos(gamma)) + V**2/r*cos(gamma)*sin(psi)*tan(phi)-2*Omega*V*(tan(gamma)*cos(psi)*cos(phi)-sin(phi))+Omega**2*r/cos(gamma)*sin(psi)*sin(phi)*cos(phi)),
          #s-dot
          -V*cos(gamma)/r
          ]
    return dyde

A=1
m=100e3
sigma=10*pi/180
Omega=7.2921159e-5

r_0=200e3
v_0=2e3
e_0=E(r_0, v_0)
r_f=2e3
v_f=200
e_f = E(r_f,v_f)

s_f=0
y_0=[r_0, 86*pi/180, 39*pi/180, -10*pi/180,0,200e3]
e=np.linspace(e_0,e_f, 1000)

sol=np.array(odeint(ode_sys,y_0,e, args=(Omega, sigma, m, A)))
def z(s):
    s_e=odeint(ode_sys,y_0,e, args=(Omega, s, m,A))
    return s_e-s_f

def initialGuess(sigma_0):
    flag=False
    sigma_k=sigma
    #while(flag==False):
    #    sigma_k=sigma_k-lambda_k*(z(sigma_k))/(z(sigma_k))
    return 1
