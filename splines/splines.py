import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class CubicSpline:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.n = len(x)
        self.a, self.b, self.c, self.d = self.compute_coefficients()

    def compute_coefficients(self):
        h = [self.x[i+1] - self.x[i] for i in range(self.n - 1)]
        alpha = [0] * self.n
        for i in range(1, self.n - 1):
            alpha[i] = 3 * ((self.y[i + 1] - self.y[i]) / h[i] - (self.y[i] - self.y[i - 1]) / h[i - 1])
        l = [1] + [0] * (self.n - 1)
        mu = [0] * self.n
        z = [0] * self.n
        for i in range(1, self.n - 1):
            l[i] = 2 * (self.x[i + 1] - self.x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / l[i]
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
        l[self.n - 1] = 1
        z[self.n - 1] = 0
        c = [0] * self.n
        b = [0] * self.n
        d = [0] * self.n
        for j in range(self.n - 2, -1, -1):
            c[j] = z[j] - mu[j] * c[j + 1]
            b[j] = (self.y[j + 1] - self.y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
            d[j] = (c[j + 1] - c[j]) / (3 * h[j])
        a = self.y[:-1]
        return a, b, c, d

    def __call__(self, x_eval):
        y_eval = []
        for x in x_eval:
            idx = 0
            for i in range(self.n - 1):
                if self.x[i] <= x <= self.x[i + 1]:
                    idx = i
                    break
            h = x - self.x[idx]
            y = self.a[idx] + self.b[idx] * h + self.c[idx] * h ** 2 + self.d[idx] * h ** 3
            y_eval.append(y)
        return y_eval

# Dados de entrada
x = np.linspace(0, 2*np.pi, 100)

# Definindo os valores de y1, y2 e y3 em função de x
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.exp(x)

# Criando splines cúbicos para cada conjunto de dados
spline1 = CubicSpline(x, y1)
spline2 = CubicSpline(x, y2)
spline3 = CubicSpline(x, y3)

# Avaliando os splines cúbicos em x
y1_interp = spline1(x)
y2_interp = spline2(x)
y3_interp = spline3(x)

# Configuração da figura 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plotando os pontos dos splines cúbicos em 3D com linhas contínuas
ax.plot(y1_interp, y2_interp, y3_interp, color='b', label='Linha dos Splines Cúbicos')

# Configuração dos eixos
ax.set_xlabel('X1')
ax.set_ylabel('X2')
ax.set_zlabel('X3')

# Adicionando legenda
plt.legend()

# Exibição do gráfico
plt.title('Linha dos Splines Cúbicos em 3D')
plt.show()
