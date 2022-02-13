import matplotlib.pyplot as plt

def integral(f, a, b, h = 1e-1):
    """Calcula a integral numérica de uma função f no intervalo a, b.
    Args:
        f (function(arg: float)): função à ser integrada
        a (float): [limite inferior de integração]
        b (float): [limite superior de integração]
        h (float): [valor considerado infinitesimal]
    Returns:
        (float): [integral da função f no intervalo a, b]
    """
    n = int((b - a) / h)
    integral = 0
    for i in range(1, n + 1):
        integral += f(a + i * h) * h
    return integral

def gera_superficie_resposta(pontos, chart_type = 'economico'):

    temperaturas = [ponto[1] for ponto in pontos]
    carga = [ponto[2] for ponto in pontos]
    z = [ponto[3] for ponto in pontos]

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    color_map = 'inferno_r' if chart_type == 'economico' else 'viridis'
    scatter = ax.scatter(temperaturas, carga, z, c = z, cmap = color_map)
    ax.set_xlabel('Tempertura (ºC)')
    ax.set_ylabel('Carga (g/L)')
    ax.set_zlabel('{z}'.format(z = 'Conversão (%)' if chart_type == 'conversao' else 'Margem Bruta (%)'))
    plt.title("{catalisador}".format(catalisador = 'Dowex 50 - WX2' if pontos[0][0] else 'Amberlite IR - 122'))
    plt.colorbar(scatter)
    plt.show()
