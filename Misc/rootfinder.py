
import scipy.optimize as sci
import numpy as np

c = 55600  # mol/m3
y = 0.7
P = 1e5  # Pa
R = 8.314  # J/mol/K
T = 298  # K

hg = 1e-3  # m2/s
hl = 1e-4  # m2/s

Ho2 = 4.63e9  # Pa
clInfO2 = 0.075  # mol/m3

sigma = 5e-5  # m
DeffO2 = 5e-11  # m2/s


def parenteses_denominador(no2):
    return y * P / R / T - no2 / hg


def denominador_ln(no2):
    return c - parenteses_denominador(no2) * R * T * c / Ho2


def numerador_ln(no2):
    return c - clInfO2 - no2 / hl


def ln(no2):
    return np.log(numerador_ln(no2) / denominador_ln(no2))


externo = DeffO2 * c / sigma


def lado_direito(no2):
    return ln(no2) * externo - no2


sci.root(lado_direito, 0).x[0]
