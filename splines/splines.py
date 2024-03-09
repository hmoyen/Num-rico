import numpy as np
import matplotlib.pyplot as plt

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
x1 = np.linspace(0, 2*np.pi, 10)
x2 = np.linspace(0, 2, 10)
x3 = np.linspace(-1, 1, 10)
y1 = np.sin(x1)
y2 = np.cos(x2)
y3 = np.exp(x3)

# Criar splines cúbicos para cada conjunto de dados
spline1 = CubicSpline(x1, y1)
spline2 = CubicSpline(x2, y2)
spline3 = CubicSpline(x3, y3)

# Gerar novos pontos para plotagem
x1_new = np.linspace(0, 2*np.pi, 100)
x2_new = np.linspace(0, 2, 100)
x3_new = np.linspace(-1, 1, 100)
y1_new = spline1(x1_new)
y2_new = spline2(x2_new)
y3_new = spline3(x3_new)

# Expressões estimadas para cada spline
expression1 = f"y1(x) = {spline1.a[0]:.2f} + {spline1.b[0]:.2f} * (x - {spline1.x[0]:.2f}) + {spline1.c[0]:.2f} * (x - {spline1.x[0]:.2f})^2 + {spline1.d[0]:.2f} * (x - {spline1.x[0]:.2f})^3"
expression2 = f"y2(x) = {spline2.a[0]:.2f} + {spline2.b[0]:.2f} * (x - {spline2.x[0]:.2f}) + {spline2.c[0]:.2f} * (x - {spline2.x[0]:.2f})^2 + {spline2.d[0]:.2f} * (x - {spline2.x[0]:.2f})^3"
expression3 = f"y3(x) = {spline3.a[0]:.2f} + {spline3.b[0]:.2f} * (x - {spline3.x[0]:.2f}) + {spline3.c[0]:.2f} * (x - {spline3.x[0]:.2f})^2 + {spline3.d[0]:.2f} * (x - {spline3.x[0]:.2f})^3"

# Plotar os resultados
plt.figure(figsize=(10, 6))
plt.plot(x1_new, y1_new, label='Spline para y1(x1)')
plt.plot(x2_new, y2_new, label='Spline para y2(x2)')
plt.plot(x3_new, y3_new, label='Spline para y3(x3)')
plt.scatter(x1, y1, color='red', label='Pontos de y1(x1)')
plt.scatter(x2, y2, color='blue', label='Pontos de y2(x2)')
plt.scatter(x3, y3, color='green', label='Pontos de y3(x3)')

# Adicionando as expressões estimadas nos títulos
plt.title(f'Interpolação por Splines Cúbicos\n{expression1}\n{expression2}\n{expression3}')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()
