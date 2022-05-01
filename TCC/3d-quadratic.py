# Plot das curvas de investigação de Extratante e Ácido, usando modelo quadrático de regressão

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from matplotlib.gridspec import GridSpec

def get_features(x, y):
    features = {}
    features['x^0*y^0'] = np.matmul(x**0, y**0).flatten()
    features['x^0*y^1'] = np.matmul(x**0, y**1).flatten()
    features['x^0*y^2'] = np.matmul(x**0, y**2).flatten()
    features['x^1*y^0'] = np.matmul(x**1, y**0).flatten()
    features['x^1*y^1'] = np.matmul(x**1, y**1).flatten()
    features['x^1*y^2'] = np.matmul(x**1, y**2).flatten()
    features['x^2*y^0'] = np.matmul(x**2, y**0).flatten()
    features['x^2*y^1'] = np.matmul(x**2, y**1).flatten()
    return pd.DataFrame(features)

def get_experimental_points(x_values, y_values):
    dim = len(x_values)
    x = np.array(x_values).reshape(dim,1)
    y = np.array(y_values).reshape(1,dim)
    X,Y = np.meshgrid(x,y)
    return X, Y, get_features(x,y)

def get_smooth_points(x_values, y_values, regression):
    dim = 50
    x = np.linspace(x_values[0], x_values[-1], dim).reshape(dim,1)
    y = np.linspace(y_values[0], y_values[-1], dim).reshape(1,dim)
    x_smooth, y_smooth = np.meshgrid(x,y)
    dataset = get_features(x,y)
    z_model = regression.intercept_ + np.matmul(dataset.values, regression.coef_.reshape(-1,1)).reshape(dim,dim)
    return x_smooth, y_smooth, z_model

def _add_subplot(fig, gridspec, cmap, z, title, ph, extractant_concentration):

    z_experimental = np.array(z)
    x_experimental, y_experimental, dataset = get_experimental_points(ph, extractant_concentration)
    regression = LinearRegression().fit(dataset.values, z_experimental.flatten())
    x_smooth, y_smooth, z_model = get_smooth_points(ph, extractant_concentration, regression)

    ax = fig.add_subplot(gridspec, projection='3d')
    ax.plot_surface(x_smooth, y_smooth, z_model, cmap=cmap, label='prediction', rstride=1, cstride=1, alpha=None)
    ax.scatter(x_experimental, y_experimental, z_experimental, c = 'r', label = 'datapoints')
    ax.view_init(42, -136)
    ax.title.set_text(title)
    ax.set_xlabel('pH')
    ax.set_ylabel('Concentração de Extratante (%v/v)')

def main(chart_title, cmap, z_props, ph, extractant_concentration):
    fig = plt.figure()
    fig.suptitle(chart_title)
    gs = GridSpec(nrows=1, ncols=len(z_props))

    for index, z_prop in enumerate(z_props):
        z, z_title = z_prop
        _add_subplot(fig, gs[0, index], cmap, z, z_title, ph, extractant_concentration)

    plt.show()

# D2EHPA-HCl
e = [74,67,72,83,84,88,90,92,94]
beta = [2.59,2.7,2.6,2.58,2.6,2.58,2.34,2.36,2.32]
ph = [0.5, 1, 1.5]
ext = [17, 26, 35]
main('D2EHPA-HCl', 'plasma_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

# D2EHPA-HNO3
e = [56,64,68,74,77,73,84,93,90]
beta = [2.56,2.19,2.20,2.46,2.00,2.36,2.66,2.64,2.61]
ph = [1, 2, 3]
ext = [17, 26, 35]
main('D2EHPA-HNO3', 'plasma_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

# P507-HCl
e = [31,34,43,47,45,54,58,58,44]
beta = [2.78,2.12,2.14,3.07,2.77,2.48,3.51,3.13,2.37]
ph = [1, 2, 3]
ext = [17, 26, 35]
main('P507-HCl', 'cividis_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

# P507-HNO3
e = [28,36,35,37,44,44,52,53,53]
beta = [3.02,2.79,3.41,4.71,2.96,3.31,2.61,2.91,2.43]
ph = [1, 2, 3]
ext = [17, 26, 35]
main('P507-HNO3', 'cividis_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

# Cyanex 272-HCl
e = [7,10,9,7,12,14,11,14,12]
beta = [3.78,1.87,3.22,3.59,2.55,2.64,2.03,1.72,2.23]
ph = [1, 2, 3]
ext = [17, 26, 35]
main('Cyanex 272-HCl', 'viridis_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

# Cyanex 272-HNO3
e = [5,15,19,15,25,23,7,26,25]
beta = [1.35,2.36,1.23,0.97,1.25,1.18,0.96,1.34,1.5]
ph = [1, 2, 3]
ext = [17, 26, 35]
main('Cyanex 272-HNO3', 'viridis_r', [(e, '% Dy Extraído'), (beta, 'Beta Dy/Tb')], ph, ext)

