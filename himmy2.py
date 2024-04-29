import cvxpy as cp
import numpy as np

def optimize(vector):
    t = cp.Variable()
    objective = cp.Minimize(cp.sum_squares(t - vector))
    prob = cp.Problem(objective)
    prob.solve()
    print(f"(a) Optimal point t* for {vector}: {t.value:.3f}")

pid_vector = np.array([7,3,0,4,7,7,8,0,3])
one_two_three_four = np.array([1,2,3,4])
onehundred = np.array([1,2,3,4,100])

optimize(pid_vector)
optimize(one_two_three_four)
optimize(onehundred)

def optimize_again(vector):
    n = len(vector)
    x = cp.Variable(n)
    y = cp.Variable(n)
    t = cp.Variable()

    objective = cp.Minimize(cp.sum(x + y))
    constraints = [
        t * np.ones(n) - vector == x - y,
        x >= 0,
        y >= 0
    ]

    prob = cp.Problem(objective, constraints)
    prob.solve()
    print(f"(a) Optimal point t* for {vector}: {t.value:.3f}")

optimize_again(pid_vector)
optimize_again(one_two_three_four)
optimize_again(onehundred)

#K.P Not unique because the median can be the same from multiple vector


