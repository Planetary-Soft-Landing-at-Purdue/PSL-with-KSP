import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
import math
import sys
from decimal import *
from numpy import linalg
import time

start_time = time.time()
w = [0, 0, 0]  # planet angular velocity
m_f = 1505  # dry mass
m_0 = 1905  # wet mass
alpha = 4.53e-4  # fuel consumption rate kg/N
rho_1 = 4972  # lower bound on thrust
rho_2 = 13260  # upper bound on thrust
r_N = [0, 0, 0]  # final position, ECEF
rdot_N = [0, 0, 0]  # final velocity, ECEF
V_max = 0
theta = 180  # pointing constraint
dt = 0.1  # time step size
gamma = 4 * (math.pi) / 180

c = np.array([[1, 0, 0]]) / math.tan(gamma)
S = np.array([[0, 1, 0],
              [0, 0, 1]])
E = np.array([[1, 0, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0]])
F = np.array([[0, 0, 0, 0, 0, 0, 1]])
E_u = np.array([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0]])
E_v = np.array([[0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0]])
cE = np.matmul(c, E)
SE = np.matmul(S, E)

x_0 = np.array([1500, 500, 2000, -75, 0, 100, np.log(m_0)])
g = np.array([-3.71, 0, 0, 0])

s_W = np.array([[(w[2] ** 2 + w[1] ** 2), -1 * w[1] * w[0], -1 * w[2] * w[0], 0, 2 * w[2], -2 * w[1]],
                [-1 * w[1] * w[0], (w[2] ** 2 + w[0] ** 2), -1 * w[2] * w[1], -2 * w[2], 0, 2 * w[0]],
                [-1 * w[2] * w[0], -1 * w[2] * w[1], (w[1] ** 2 + w[0] ** 2), 2 * w[1], -2 * w[0], 0],
                ])

A = np.concatenate((np.zeros((3, 3)), np.identity(3)), axis=1)
A = np.concatenate((A, s_W, np.zeros((1, 6))))
A = np.concatenate((A, np.zeros((7, 1))), axis=1)

B = np.concatenate((np.zeros((3, 3)), np.identity(3), np.zeros((1, 3))), axis=0)
B = np.concatenate((B, np.array([[0, 0, 0, 0, 0, 0, -alpha]]).T), axis=1)

C = np.array([[0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0]])


def eMatrix(A, x):
    test = np.zeros((7, 7))
    for i in range(10):
        test = test + np.linalg.matrix_power(A, i) * x ** i / math.factorial(i)
    return np.array(test, dtype=float)


def Upsilon(n, i):
    test = np.zeros((4, 4 * n))
    test[0:, (4 * i):(4 * i + 4)] = np.identity(4)
    return test


def e(i, n):
    test = np.zeros((n, 1))
    test[i - 1] = 1
    return test


def matExp(A, x):  # finds matrix exponential using taylor polynomials, A -> matrix, x -> scalar that A is multiplied by
    expMat = np.zeros_like(A)
    for i in range(30):
        expMat = expMat + np.linalg.matrix_power(A, i) * x ** i / math.factorial(i)
    return expMat


accuracy = 500
x_r = matExp(A, dt)
# x_c = 0.5*dt*((matExp( A, 0) + x_r) @ B)
x_k = np.array(x_0)
x_c = 0.5 * dt / accuracy * np.sum(
    [(matExp(A, tau / accuracy * dt) + matExp(A, (tau / accuracy - dt / accuracy) * dt)) @ B for tau in
     range(1, accuracy)], axis=0)

for k in range(1, int(30 / dt)):

    # x_g = 0.5*dt*((matExp(A, k*dt) + matExp(A, (-k)*dt)) @ g)

    x_k = x_r @ x_k + x_c @ g  # state vector
    if (x_k[0] <= 1):
        print(str(k * dt) + " seconds to hit the ground")
        break

print("\nHeader compiled in %s seconds.\n" % (time.time() - start_time))