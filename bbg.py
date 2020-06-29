import numpy as np 
from math import cos, sin, sqrt

g      = 9.81
zeta   = .5 * (3-sqrt(5))
h_taem = 1000

def find_range(state0, t, theta_c):
	x0, h0, u0, v0, m0, T, dm, = state0

	E_t = math.log(1 - (dm * t) / m0)
	F_t = (1 - (dm * t) / m0) * (E_t - 1) + 1

	x_t = x0 + u0 * t + ((m0 * T * cos(theta_c)) / dm**2) * F_t
	h_t = h0 + v0 * t - (.5 * 9.81 * t**2) + ((m0 * T * sin(theta_c)) / dm**2) * F_t
	u_t = u0 - ((T * cos(theta_c)) / dm) * E_t
	v_t = v0 - g * t - ((T * sin(theta_c)) / dm) * E_t 

	return (-u_t / g) * (v_t + sqrt(v_t ** 2 - 2 * g * (h_taem - hf)))

def optimize_theta_c(state0, tf):
	thLow, thHigh = 120, 220

	th1 = int(thHigh - zeta * (thHigh - thLow))
    th2 = int(thLow + zeta * (thHigh - thLow))

    f1 = optimize_theta_c(state0, tf, th1)
    f2 = optimize_theta_c(state0, tf, th2)

    while th2 - th1 > 1:
    	if f1 < f2:
    		tLow, th1, f1, th2, f2 = 
    			th1, th2, f2, int(tLow + zeta * (tHigh - tLow)), optimize_theta_c(state0, tf, th2)
    	else:
    		tHigh, th2, f2, th1, f1 = 
    			th2, th1, f1, int(tHigh - zeta * (tHigh - tLow)), optimize_theta_c(state0, tf, th1)

    return abs(find_range(state0))

def optimize_tf(state0):
	tLow, tHigh = 0, 50

	t1 = int(tHigh - zeta * (tHigh - tLow))
    t2 = int(tLow + zeta * (tHigh - tLow))

    f1 = optimize_theta_c(state0, t1, theta_c)
    f2 = optimize_theta_c(state0, t2, theta_c)

    while t2 - t1 > 1:
    	if f1 > f2:
    		tLow, t1, f1, t2, f2 = 
    			t1, t2, f2, int(tLow + zeta * (tHigh - tLow)), optimize_theta_c(state0, t2)
    	else:
    		tHigh, t2, f2, t1, f1 = 
    			t2, t1, f1, int(tHigh - zeta * (tHigh - tLow)), optimize_theta_c(state0, t1)

    return optimize_theta_c(state0, .5 * (t2 - t1), theta_c)


    	








