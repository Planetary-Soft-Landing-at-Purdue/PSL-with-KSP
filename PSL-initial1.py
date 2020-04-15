import numpy as np
import cvxpy as cp
import math
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d

steps, dimensions = 100, 3
X, Y, Z = 50, 2, -10
theta   = math.tan(math.pi / 4)

xyz = cp.Variable((steps, dimensions))
x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
parameters  = (cp.sum(x) - X)**2 + (cp.sum(y) - Y)**2 + (cp.sum(z) - Z)**2
objective   = cp.Minimize(parameters)
constraints = [z/theta <= x, x <= -z/theta,
               z/theta <= y, y <= -z/theta,
               sum(z) >= Z]

prob   = cp.Problem(objective, constraints)
result = prob.solve()

xyzSteps, xyzPos = np.zeros((steps, dimensions)), np.zeros((1, dimensions))
for s in range(steps):
    xyzPos += xyz.value[s]
    xyzSteps[s] = xyzPos

xSteps, ySteps, zSteps = xyzSteps[:, 0], xyzSteps[:, 1], xyzSteps[:, 2]
plt.figure(1)
plt.plot(xSteps, zSteps)
plt.xlabel("x-direction")
plt.ylabel("z-direction")
plt.figure(2)
plt.plot(ySteps, zSteps)
plt.xlabel("y-direction")
plt.ylabel("z-direction")