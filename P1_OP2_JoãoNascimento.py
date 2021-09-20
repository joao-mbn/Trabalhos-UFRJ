import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

def alfa(z):
    """volatilidade relativa das fases alfa(x) = a*x + b ou alfa(y) = a*y + b"""
    return alfa_a * z + alfa_b

# ==========================
# ========== MPD ===========
# ==========================

def h(y):
    """RELV X = h(Y)"""
    return y / (y + (1 - y) * alfa(y))

def op_MPD(x):

    for i in range(len(xy_emendas)):
        if x > xy_emendas[i][0]:
            a = coefs_eops[i][0]
            b = coefs_eops[i][1]
            return a * x + b
            
def estagios_MPD(x0, y1, xn, eficiencia = 1):

    ys = [None, y1]
    xs = [x0]

    y_plot = [y1, y1]
    x_plot = [x0]

    while xs[-1] > xn:
        x = (1-eficiencia) * xs[-1] + h(ys[-1])
        y = op_MPD(x)
        xs.append(x)
        ys.append(y)

        x_plot.append(x)
        x_plot.append(x)
        y_plot.append(y)
        y_plot.append(y)

    del y_plot[-1]

    return xs, ys, x_plot, y_plot

# ==========================
# ========== MPA ===========
# ==========================

def g(x):
    """RELV Y = g(X)"""
    return alfa(x) * x / (1 - x + alfa(x) * x)
    
def op_MPA(y):

    for i in range(len(xy_emendas)):
        if y < xy_emendas[i][1]:
            a = coefs_eops[i][0]
            b = coefs_eops[i][1]
            return a * y + b

def estagios_MPA(x1, y0, yn, eficiencia = 1):

    ys = [y0]
    xs = [None, x1]

    y_plot = [y0]
    x_plot = [x1, x1]

    while ys[-1] < yn:
        y = (1-eficiencia) * ys[-1] + g(xs[-1])
        x = op_MPA(y)
        xs.append(x)
        ys.append(y)

        x_plot.append(x)
        x_plot.append(x)
        y_plot.append(y)
        y_plot.append(y)

    del x_plot[-1]

    return xs, ys, x_plot, y_plot

def graficos_mct(x_plot, y_plot, relv_x, relv_y, eops):
    """Plot do Gráfico de Mc Cabe-Thiele, 
        os estágios,
        a reta de operação
        e a RELV"""
    plt.plot(x_plot, y_plot)
    plt.plot(relv_x, relv_y)
    for eop in eops:
        plt.plot(eop[0], eop[1])
    plt.show()

# ======================================================
# ======================================================
# ======================================================

def main(MPD):
    """Função principal: 
    gera os pontos e cria os gráficos"""

    if MPD:
        xs, ys, x_plot, y_plot = estagios_MPD(x0, y1, xn)

        relv_y = np.linspace(0, 1, 50)
        relv_x = h(relv_y)

        eops = []

        for i, coefs_eop in enumerate(coefs_eops):
            if i == 0:
                xs = np.linspace(xy_emendas[i][0], x0)
            else:
                xs = np.linspace(xy_emendas[i-1][0], xy_emendas[i][0])    
            
            a = coefs_eop[0]
            b = coefs_eop[1]
            ys = a * xs + b
            eops.append([xs, ys])
    
    else:
        xs, ys, x_plot, y_plot = estagios_MPA(x1, y0, yn)

        relv_x = np.linspace(0, 1, 50)
        relv_y = g(relv_x)

        eops = []

        for i, coefs_eop in enumerate(coefs_eops):
            if i == 0:
                ys = np.linspace(y0, xy_emendas[i][1])
            else:
                ys = np.linspace(xy_emendas[i-1][1], xy_emendas[i][1])     
            
            a = coefs_eop[0]
            b = coefs_eop[1]
            xs = a * ys + b
            eops.append([xs, ys])

    graficos_mct(x_plot, y_plot, relv_x, relv_y, eops)

# =========================================
# ======== Declaração de Variáveis ========
# =========================================

# Template MPD

# primeiro em cima, último embaixo

coefs_eops = [[0.8, 0.198], [1.58, -0.29145]]
xy_emendas = [[0.6275, 0.7], [0.5025, 0.5025]]
alfa_a, alfa_b = 2.183, 0.943

x0 = 0.99
y1 = 0.99
xn = xy_emendas[-1][0]

main(True)

# Template MPA

# primeiro embaixo, último em cima

coefs_eops = [[0.49395, 0.0050605], [0.9879, -0.08988]]
xy_emendas = [[0.1, 0.1922], [0.5971, 0.5971]]
alfa_a, alfa_b = 2.225, 3.139

x1 = 0.01
y0 = 0.01
yn = xy_emendas[-1][1]

main(False)