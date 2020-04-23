import cvxpy as cp

# Create scalar optimization variable.
x = cp.Variable()

# Create constraint.
constraints = [3<=x, x<=6]

# Form objective.
obj = cp.Minimize((x)**2)

# Form and solve problem.
prob = cp.Problem(obj, constraints)
prob.solve()  # Returns the optimal value.
print("status:", prob.status)
print("optimal value", prob.value)
print("optimal var", x.value)