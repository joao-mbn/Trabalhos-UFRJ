from scipy.optimize import fsolve
import numpy as np

a = 1.5e-4
ca2f = 1000
ca2i = 312.5

def equations(p):
    t2, q = p
    equation1 = 0.3/t2 - q
    equation2 = a * (ca2f - ca2i * np.exp(-a * t2)) / 2500 / (1 - np.exp(-a * t2)) - q
    return (equation1, equation2)

t2, q =  fsolve(equations, (500, 7))


