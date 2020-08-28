from scipy.optimize import newton
from numpy import sin as s, cos as c, log, power as pow, sqrt, pi
import numpy as np

# read state
x0 = 35000; h0 = 75000; u0 = 500; v0 = 1600
m0 = 29000; md = 280*3; T = 800000*3; g = 9.8

def E(t): return log(1-md*t/m0)
def F(t): return (1-md*t/m0)*(E(t)-1)+1

def x(θ, t): return x0 + u0*t + m0*T*c(θ)*F(t)/pow(md,2)
def h(θ, t): return h0 + v0*t - 0.5*g*pow(t,2) + m0*T*s(θ)*F(t)/pow(md,2)
def u(θ, t): return u0 - T*c(θ)*E(t)/md
def v(θ, t): return v0 - g*t - T*s(θ)*E(t)/md

def hθ(θ, t): return m0*T*c(θ)*F(t)/pow(md,2)
def uθ(θ, t): return T*s(θ)*E(t)/md
def vθ(θ, t): return -T*c(θ)*E(t)/md

def func_θ(θc, tf):
    hf = h(θc, tf)
    uf = u(θc, tf)
    vf = v(θc, tf)
    return (
        uθ(θc,tf) * (vf*sqrt(pow(vf,2)+2*g*hf) + pow(vf,2) + 2*g*hf) +
        uf*(vθ(θc,tf)*(vf+sqrt(pow(vf,2)+2*g*hf)) + g*hθ(θc,tf))
    )

def x_max(θc, tf):
    hf = h(θc, tf)
    uf = u(θc, tf)
    vf = v(θc, tf)
    h_pdg = 0
    return -(uf/g) * (vf + sqrt(pow(vf,2) - 2*g*(h_pdg-hf)))


# secant search optimal t-cutoff
tf1 = 1; tf2 = 2
θc1 = newton(func_θ, pi, args=(tf1,))
θc2 = newton(func_θ, pi, args=(tf2,))
xm1 = x_max(θc1, tf1)
xm2 = x_max(θc2, tf2)
while abs(tf2 - tf1) > 0.001:
    tf3 = tf2 - ((tf2 - tf1) / (xm2 - xm1)) * (xm2 - x0)
    tf1 = tf2; tf2 = tf3
    θc1 = newton(func_θ, pi, args=(tf1,))
    θc2 = newton(func_θ, pi, args=(tf2,))
    xm1 = x_max(θc1, tf1)
    xm2 = x_max(θc2, tf2)

# optimal thrust angle and t-cutoff
θc = θc2
tf = tf2

