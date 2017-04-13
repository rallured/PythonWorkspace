import numpy as np
import matplotlib.pyplot as plt
import cvxpy

def prob4_2():
    x1 = cvxpy.Variable()
    x2 = cvxpy.Variable()
    objective = cvxpy.Minimize(x1**2+9*x2**2)
    constraints = [2*x1+x2>=1, x1+3*x2>=1, x1 >=0, x2 >=0]
    prob = cvxpy.Problem(objective,constraints)
    prob.solve()
    print prob.value, x1.value, x2.value

    return

def prob4_4():
    n=100
    m=300
    A = np.random.uniform(size=(m,n))
    b = np.dot(A,np.ones(n)/2.)
    c = -np.random.uniform(size=n)

    x = cvxpy.Variable(n)
    objective = cvxpy.Minimize(c*x)
    constraints = [0 <= x, x <= 1, A*x<=b]

    prob = cvxpy.Problem(objective,constraints)

    prob.solve()
    print prob.value, x.value

    return
