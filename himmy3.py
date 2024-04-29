import cvxpy as cp
import numpy as np

C = np.array([0.5,0.5], [0.5, 0.5])

n = C.shape[0]
P = cp.Variable((n, n))

objective = cp.Minimize(cp.norm(C-P, "fro"))

constraints = [
    cp.sum(P, axis=0) == np.ones(n),  
    cp.sum(P, axis=1) == np.ones(n),  
    P >= 0  
]

problem = cp.Problem(objective, constraints)
problem.solve()

print("The optimal value is:", problem.value)
print("The optimal matrix P* is:")
print(P.value)