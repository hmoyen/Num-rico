import numpy as np
from scipy.optimize import curve_fit

# Velocidade angular (rad/s)
omega = np.pi / 6

# Raio do círculo
R = 1  # Assumindo raio 1

# Tempos
tempos = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

# Calculando a posição no eixo x em função do tempo
velocidade_x = np.exp(-2*tempos)

# Calculando a aceleração linear no eixo x em função do tempo
aceleracao_x = 4*np.exp(-2*tempos)

# Imprimir os arrays
print("Posição no eixo x em função do tempo:")
print(velocidade_x)
print("\nAceleração linear no eixo x em função do tempo:")
print(aceleracao_x)

# Função que queremos ajustar
def func_to_fit(data, a, b, c):
    t, y = data
    return a * t + b * y + c

# Emparelhando os dados para o curve_fit
data = (velocidade_x, tempos)

# Ajustando a função aos dados usando curve_fit
popt, pcov = curve_fit(func_to_fit, data, aceleracao_x)

# Parâmetros ajustados
a, b, c = popt

# Calculando os valores previstos pela função ajustada
z_pred = func_to_fit(data, a, b, c)

# Imprimindo os parâmetros ajustados
print("Parâmetros Ajustados:")
print("a =", a)
print("b =", b)
print("c =", c)



