import numpy as np

# =====================================
# =====================================
# ============= Questão 1 =============
# =====================================
# =====================================

delta = 1e-4

y_delta1 = 0.750
y_delta2 = 0.2
y_delta3 = 0.005
y_delta4 = 0.045

y_0_1 = 0.728
y_0_2 = 0.137
y_0_3 = 0.080
y_0_4 = 0.055

d12 = 22e-6
d13 = 17e-6
d14 = 23e-6
d23 = 16e-6
d24 = 23e-6
d34 = 16e-6

c = 40.36

def dyi_dz(y_delta, y_0, delta): return (y_delta - y_0) / delta

dy1_dz = dyi_dz(y_delta1, y_0_1, delta)
dy2_dz = dyi_dz(y_delta2, y_0_2, delta)
dy3_dz = dyi_dz(y_delta3, y_0_3, delta)
dy4_dz = dyi_dz(y_delta4, y_0_4, delta)

y1 = np.average([y_delta1, y_0_1])
y2 = np.average([y_delta2, y_0_2])
y3 = np.average([y_delta3, y_0_3])
y4 = np.average([y_delta4, y_0_4])


a11, a12, a13 = round((y1/d12 + y3/d23 + y4/d24)/c, 4), round((-y2/d23)/c, 4), round((-y2/d24)/c, 4)
a21, a22, a23 = round((-y3/d23)/c, 4), round((y1/d13 + y2/d23 + y4/d34)/c, 4), round((-y3/d34)/c, 4)
a31, a32, a33 = round((-y4/d24)/c, 4), round((-y4/d34)/c, 4), round((y1/d14 + y2/d24 + y3/d34)/c, 4)

A = np.array([
    [a11, a12, a13],
    [a21, a22, a23],
    [a31, a32, a33]
])
b = np.array([-dy2_dz, -dy3_dz, -dy4_dz])

n2, n3, n4 = np.linalg.solve(A, b)
print(n2, n3, n4)