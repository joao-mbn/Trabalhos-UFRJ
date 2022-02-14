import numpy as np
import helpers

def calcula_custos_ganhos(f1, T_reacao, concentracao_catalisador, isDowex, rendimento_reacao):

    # ======================================================================
    # Composições molares de uma mistura binária. O subscrito de refere à coluna. A fração explícita é sempre o IbAc

    # 1º Coluna
    xf1 = rendimento_reacao
    xb1 = 1
    xd1 = 0.11

    # 2º Coluna
    xf2 = xd1
    xb2 = 0
    xd2 = 0.5085

    # mols de IbAc puro no fundo da 1º coluna
    def b1(): return f1 * (xf1 - xd1) / (xb1 - xd1)
    def d1(): return f1 - b1()  # mols de IbAc + IbOH p/ a 2º coluna

    # mols de IbOH puro no fundo da 2º coluna
    def b2(): return d1() * (xf2 - xd2) / (xb2 - xd2)
    def d2(): return d1() - b2()  # mols da blena IbAc + IbOH

    # =========================================================================
    # cálculos ponderais
    # =========================================================================

    razao_g_ton = 1e-6
    razao_mol_kmol = 1e-3
    razao_cm3_l = 1e-3

    mm_ibac = 116.16  # g/mol
    mm_iboh = 74.122  # g/mol
    mm_hac = 60.052  # g/mol
    mm_mgso4 = 120.366 # g/mol
    mm_naoh = 39.997 # g/mol

    d_hac = 1.05  # g/cm3
    d_iboh = 0.802  # g/cm3
    pureza_iboh = 0.985
    pureza_hac = 0.99
    numero_hidratacao_mgso4 = 7

    vida_util_catalisador = 730 # dias de uso (2 anos)
    vida_util_mgso4 = vida_util_catalisador
    n_bateladas = 5 # por dia

    # Dados das correntes de entrada
    # massa das correntes de entrada, em g
    m_hac = f1 * mm_hac
    m_iboh = f1 * mm_iboh
    # volume das correntes de entrada, em L
    vol_hac = m_hac / d_hac * razao_cm3_l
    vol_iboh = m_iboh / d_iboh * razao_cm3_l

    """
    Assumiu-se que 10% do HAc não reagido não foi pego na decantação e não irá para a destilação na 3º coluna, mas será capturado por uma solução de NaOH.
    Assumiu-se que 10% da água produzida não foi pega na decantação e não irá para a destilação na 3º coluna, mas será capturado por MgSO4 seco.
    """
    hac_lavado = 0.1 * f1 * (1- rendimento_reacao) # mols.
    agua_nao_decantada = 0.1 * f1 * rendimento_reacao # mols.

    # quantidade dos sólidos: MgSO4 e catalisador são recuperados. NaOH é consumido.
    """espera-se que as espécies sejam recuperadas entre uma batelada e outra,
    logo elas precisam arcar com o volume diário das correntes"""
    quantidade_catalisador = concentracao_catalisador * (vol_hac + vol_iboh) / n_bateladas # g
    quantidade_mgso4 = agua_nao_decantada * mm_mgso4 / (numero_hidratacao_mgso4 * n_bateladas) # g
    quantidade_naoh = hac_lavado * mm_naoh # g

    # conversão das quantidades de matéria das correntes de mol p/ g.
    def m_corrente(x, mols_corrente): return (x * mm_ibac + (1 - x) * mm_iboh) * mols_corrente
    # massa das correntes das colunas, em g
    m_b1 = m_corrente(xb1, b1())
    m_b2 = m_corrente(xb2, b2())
    m_d2 = m_corrente(xd2, d2())
    m_b3 = m_hac * (1 - rendimento_reacao) - hac_lavado

    # =========================================================================
    # cálculos energéticos
    # =========================================================================

    razao_joule_kwh = 1 / 3.6e6

    T_ambiente = 25 #ºC
    T_max_efluente = 60 #ºC
    T_eb_hac = 118 #ºC
    T_eb_iboh = 108 #ºC
    T_eb_ibac = 117 #ºC

    # Entalpia de formação de gás ideal em J/kmol, 298K, 1bar
    entalpia_formacao_hac = -16.64e7
    entalpia_formacao_iboh = -27.51e7 # Valores do 1-Butanol
    entalpia_formacao_agua = -24.1814e7
    entalpia_formacao_ibac = -48.56e7 # Valores do acetato de butila

    def entalpia_vaporizacao(T, Tc, c1, c2, c3, c4):
        """T(K), resultado em J/kmol"""
        Tr = T/Tc
        return c1 * (1 - Tr) ** (c2 + c3 * Tr + c4 * Tr ** 2)

    def entalpia_reacao(T):
        """T(K), resultado em J/kmol"""
        entalpia_formacao_hac_liq = entalpia_formacao_hac - entalpia_vaporizacao(T, 591.95, 4.0179e7, 2.6037, -5.0031, 2.7069)
        entalpia_formacao_iboh_liq = entalpia_formacao_iboh - entalpia_vaporizacao(T, 563.1, 7.1274e7, 0.0483, 0.8966, -0.5116) # Valores do 1-Butanol
        entalpia_formacao_agua_liq = entalpia_formacao_agua - entalpia_vaporizacao(T, 647.093, 5.2053e7, 0.3199, -0.212, 0.25795)
        entalpia_formacao_ibac_liq = entalpia_formacao_ibac - entalpia_vaporizacao(T, 575.4, 5.8276e7, 0.38854, 0, 0) # Valores do acetato de butila
        return entalpia_formacao_ibac_liq + entalpia_formacao_agua_liq - entalpia_formacao_iboh_liq - entalpia_formacao_hac_liq

    def cp(T, c1, c2, c3, c4, c5):
        """T(K), resposta em J/(kmol K)"""
        return c1 + c2*T + c3*T**2 + c4*T**3 + c5*T**4

    def cp_massico(T, massa_molar, c1, c2, c3, c4, c5):
        """T(K), resposta em J/(g K)"""
        return cp(T, c1, c2, c3, c4, c5) * razao_mol_kmol / massa_molar

    def calor(cp, massa, T_inferior, T_superior): return massa * helpers.integral(cp, T_inferior, T_superior)

    def cp_hac(T): return cp_massico(T, mm_hac, 139640, -320.8, 0.8985, 0, 0) # J/(g K)
    def cp_iboh(T): return cp_massico(T, mm_iboh, 191200, -730.4, 2.2998, 0, 0) # J/(g K). Constantes do 1-Butanol
    def cp_ibac(T): return cp_massico(T, mm_ibac, 111850, 384.52, 0, 0, 0) # J/(g K). Constantes do acetato de butila
    cp_agua = 0.001 # kWh / kg ºC

    # calor consumido na reação, em J
    calor_reacao = entalpia_reacao(T_reacao + 273) * razao_mol_kmol * f1 * rendimento_reacao
    # calor para aquecer às correntes de entrada no reator, de T ambiente para T reação, em J
    calor_hac = calor(cp_hac, m_hac, T_ambiente + 273, T_reacao + 273)
    calor_iboh = calor(cp_iboh, m_iboh, T_ambiente + 273, T_reacao +273)
    # calor para aquecer às correntes das destilações, em J
    calor_hac_rec = calor(cp_hac, m_b3, T_ambiente + 273, T_eb_hac + 273)
    calor_iboh_rec = calor(cp_iboh, m_b2, T_ambiente + 273, T_eb_iboh + 273)
    calor_ibac = calor(cp_ibac, m_b1, T_ambiente + 273, T_eb_ibac + 273)

    calor_correntes = calor_hac + calor_iboh + calor_ibac + calor_hac_rec + calor_iboh_rec # J
    calor_total = calor_reacao + calor_correntes # J

    m_agua_resfriamento_reator = -calor_reacao * razao_joule_kwh / ((T_max_efluente - T_ambiente) * cp_agua) # Kg
    m_agua_resfriamento_correntes = calor_correntes * razao_joule_kwh / ((T_max_efluente - T_ambiente) * cp_agua) # Kg
    m_agua_troca_calor = m_agua_resfriamento_correntes + m_agua_resfriamento_reator # Kg

    # =========================================================================
    # cálculos econômicos
    # =========================================================================

    preco_ibac = 1500 # USD/ton
    preco_blenda = 0.2 * preco_ibac # USD/ton
    preco_iboh = 500 # USD/ton
    preco_hac = 620 # USD/ton
    preco_energia = 0.24 # USD/kWh
    preco_agua_troca_calor = 7.6e-4 # USD/kg
    preco_amberlite = 5 # USD/g
    preco_dowex = 7.5 # USD/g
    preco_naoh = 800 # USD/ton
    preco_mgso4 = 300 # USD/ton

    vendas = (preco_ibac * m_b1 + preco_blenda * m_d2) * razao_g_ton
    custo_reagentes = (preco_hac * m_hac / pureza_hac + preco_iboh * m_iboh / pureza_iboh) * razao_g_ton
    ganho_reagentes_recuperados =  (preco_hac * m_b3 + preco_iboh * m_b2) * razao_g_ton
    custo_naoh = quantidade_naoh * preco_naoh * razao_g_ton
    custo_mgso4 = quantidade_mgso4 * preco_mgso4 * razao_g_ton / vida_util_mgso4
    custo_catalisador = quantidade_catalisador * (preco_dowex if isDowex else preco_amberlite) / vida_util_catalisador
    custo_energia = (calor_total) * preco_energia * razao_joule_kwh
    custo_agua_troca_calor = m_agua_troca_calor * preco_agua_troca_calor

    margem_bruta_percentual = 100 * (
            vendas
            - custo_reagentes
            + ganho_reagentes_recuperados
            - custo_naoh
            - custo_mgso4
            - custo_catalisador
            - custo_energia
            - custo_agua_troca_calor
        ) / vendas

    return margem_bruta_percentual,\
        vendas,\
        custo_reagentes,\
        ganho_reagentes_recuperados,\
        custo_naoh,\
        custo_mgso4,\
        custo_catalisador,\
        custo_energia,\
        custo_agua_troca_calor

# =========================================================================

def acha_otimo(f1, esta_usando_Dowex, faixa_temperatura, faixa_catalisador):

    margem_bruta_percentual = 0
    vendas = 0
    custo_reagentes = 0
    ganho_reagentes_recuperados = 0
    custo_naoh = 0
    custo_mgso4 = 0
    custo_catalisador = 0
    custo_energia = 0
    custo_agua_troca_calor = 0
    isDowex = 0
    temperatura_reacao = 0
    concentracao_catalisador = 0

    pontos_dispersao = []

    for _isDowex in esta_usando_Dowex:
        for _temperatura_reacao in faixa_temperatura:
            for _concentracao_catalisador in faixa_catalisador:

                _rendimento_reacao = eficiencia_reacao(_isDowex, _temperatura_reacao, _concentracao_catalisador) / 100

                _margem_bruta_percentual,\
                _vendas,\
                _custo_reagentes,\
                _ganho_reagentes_recuperados,\
                _custo_naoh,\
                _custo_mgso4,\
                _custo_catalisador,\
                _custo_energia,\
                _custo_agua_troca_calor\
                    = calcula_custos_ganhos(f1, _temperatura_reacao, _concentracao_catalisador, _isDowex, _rendimento_reacao)

                pontos_dispersao.append([_isDowex, _temperatura_reacao, _concentracao_catalisador, _margem_bruta_percentual])

                if _margem_bruta_percentual > margem_bruta_percentual:
                    margem_bruta_percentual = _margem_bruta_percentual
                    vendas = _vendas
                    custo_reagentes = _custo_reagentes
                    ganho_reagentes_recuperados = _ganho_reagentes_recuperados
                    custo_catalisador = _custo_catalisador
                    custo_naoh = _custo_naoh
                    custo_mgso4 = _custo_mgso4
                    custo_energia = _custo_energia
                    custo_agua_troca_calor = _custo_agua_troca_calor
                    isDowex = _isDowex
                    temperatura_reacao = _temperatura_reacao
                    concentracao_catalisador = _concentracao_catalisador

    return pontos_dispersao,\
        {
            "Catalisador Utilizado": "Dowex 50 - WX2" if isDowex else "Amberlite IR - 122",
            "Concentração do Catalisador (g/L)": round(concentracao_catalisador, 2),
            "Temperatura de Reação (ºC)": round(temperatura_reacao, 2),
            "Margem Bruta (%)": round(margem_bruta_percentual, 2),
            "Vendas (USD)": round(vendas, 2),
            "Custo de Reagentes (USD)": round(custo_reagentes, 2),
            "Economia com Reagentes Recuperados (USD)": round(ganho_reagentes_recuperados, 2),
            "Custo do NaOH (USD)": round(custo_naoh, 2),
            "Custo do MgSO4 (USD)": round(custo_mgso4, 2),
            "Amortização do Catalisador (USD)": round(custo_catalisador, 2),
            "Custo de Energia (USD)": round(custo_energia, 2),
            "Custo de Água para Troca de Calor (USD)": round(custo_agua_troca_calor, 2)
        }

def eficiencia_reacao(isDowex, temperatura, concentracao_catalisador):
    if isDowex:
        return 59.18977 + 0.14273 * temperatura + 0.013917 * concentracao_catalisador + 1.7e-3 * temperatura * concentracao_catalisador
    else:
        return 40.34083 + 0.19633 * temperatura + 0.80342 * concentracao_catalisador + 1.7e-3 * temperatura * concentracao_catalisador

def gera_pontos_dispersao_conversao(temperaturas, faixa_catalisador):
    eficiencia_reacao_dowex = []
    eficiencia_reacao_amberlite = []

    for temperatura in temperaturas:
        for concentracao_catalisador in faixa_catalisador:
            eficiencia_reacao_dowex.append([True, temperatura, concentracao_catalisador, eficiencia_reacao(True, temperatura, concentracao_catalisador)])
            eficiencia_reacao_amberlite.append([False, temperatura, concentracao_catalisador, eficiencia_reacao(False, temperatura, concentracao_catalisador)])

    return eficiencia_reacao_dowex, eficiencia_reacao_amberlite

def main():
    """
    f1 é base de cálculo em mols/dia --> corrente de entrada de um reagente no reator e a corrente total na 1º coluna de destilação
    ela é o número arredondado para uma estimativa de produção do acetado de isobutila:

    35000 Ton/Ano de output --> 134 Ton/dia = 134e6 g/dia ~ 1.2e6 mols/dia
    """

    f1 = 1.2e6 # mols/dia

    limite_inferior_temperatura = 45
    limite_superior_temperatura = 100

    limite_inferior_catalisador = 5
    limite_superior_catalisador = 15

    esta_usando_Dowex = [True, False]
    faixa_temperatura = np.linspace(limite_inferior_temperatura, limite_superior_temperatura, limite_superior_temperatura - limite_inferior_temperatura + 1)
    faixa_catalisador = np.linspace(limite_inferior_catalisador, limite_superior_catalisador, (limite_superior_catalisador - limite_inferior_catalisador) * 10 + 1)
    eficiencia_reacao_dowex, eficiencia_reacao_amberlite = gera_pontos_dispersao_conversao(faixa_temperatura, faixa_catalisador)

    helpers.gera_superficie_resposta(eficiencia_reacao_dowex, 'conversao')
    helpers.gera_superficie_resposta(eficiencia_reacao_amberlite, 'conversao')

    pontos_dispersao, condicao_otima = acha_otimo(f1, esta_usando_Dowex, faixa_temperatura, faixa_catalisador)
    pontos_dispersao_dowex = list(filter(lambda ponto: ponto[0] == True, pontos_dispersao))
    pontos_dispersao_amberlite = list(filter(lambda ponto: ponto[0] == False, pontos_dispersao))

    helpers.gera_superficie_resposta(pontos_dispersao_dowex)
    helpers.gera_superficie_resposta(pontos_dispersao_amberlite)

    return condicao_otima

melhor_condicao = main()
# ========================================================================================================
# ========================================================================================================
# ========================================================================================================