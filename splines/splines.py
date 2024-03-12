import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import matplotlib.pyplot as plt
from math import exp, sqrt 

__date__ = datetime(2019, 6, 6) # or version string or something
__author__ = "Joshua Simon"

"""Four_Step_Runge_Kutta_ODE1.py 

Implementation of the classic fourth-order method also refered as the
"original" Runge–Kutta method. This method is an implicit four step
Runge-Kutta method which solves an intial value problem numerically. 
"""
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
    
    def calc_func(self,t,d):

        # Derivada d nao aparece em funcoes que dependem apenas do tempo
        idx = 0
        for i in range(self.n - 1):
            if self.x[i] <= t <= self.x[i + 1]:
                idx = i
                break
        h = t - self.x[idx]
        print("H:", h)
        print("t:", t)
        print("self.x", self.x[idx])
        print("Length",len(self.x), "Index", idx)
        print("Parameters",self.a[idx], self.b[idx], self.c[idx], self.d[idx])
        y = self.a[idx] + self.b[idx] * h + self.c[idx] * h ** 2 + self.d[idx] * h ** 3

        return y
    

    def __call__(self, x_eval):
        y_eval = []
        for x in x_eval:
            idx = 0
            for i in range(self.n - 1):
                if self.x[i] <= x <= self.x[i + 1]:
                    idx = i
                    break
            h = x - self.x[idx]
            print("xzinho", self.x[idx])
            print("Hzinho", h)
            y = self.a[idx] + self.b[idx] * h + self.c[idx] * h ** 2 + self.d[idx] * h ** 3
            y_eval.append(y)
        return y_eval

def runge_kutta(f, x_0, y_0, h):
    """Four step Runge-Kutta method (RK4)
    Solves first order ODEs
    """
    k_0 = f(x_0, y_0)
    k_1 = f(x_0 + h/2, y_0 + h/2 * k_0)
    k_2 = f(x_0 + h/2, y_0 + h/2 * k_1)
    k_3 = f(x_0 + h, y_0 + h * k_2)

    k = 1/6 * (k_0 + 2.0*k_1 + 2.0*k_2 + k_3)

    x_1 = x_0 + h
    # y_1 = y_0 + h * k

    return x_1, k


def f(t, y):
    """Example first order ordinary differential equation (ODE)"""
    
    # Dados de entrada
    x = np.linspace(0, 2*np.pi, 100)

    y3 = np.cos(x)

    # Criando splines cúbicos para cada conjunto de dados
    spline3 = CubicSpline(x, y3)
    # y3_interp = spline3(x)

    # Avaliando os splines cúbicos em x
    y3_func = spline3.calc_func(t,y)
    print("Function;",y3_func)

    return y3_func

def v_func(values):

    x = np.linspace(0, 2*np.pi, 100)
    splineV = CubicSpline(x,values)
    # func= splineV.calc_func(t,y)
    return splineV.calc_func()


def pos_linear_model(p,v,a,delta):

    return (p+v*delta+(a*delta**2)*0.5)

def vel_linear_model(v,a,delta):

    return (v+a*delta)


if __name__=="__main__":
    # Initial values
    t_0 = 0.0
    v_0 = 0.0
    x_0 = -1

    # Step length 15Hz aproximadamente
    h = 2*np.pi/100 


    t_values = [t_0]
    v_values = [v_0]
    x_values = [x_0]

    # Calculate solution
    t = t_0
    v = v_0
    x = x_0

    contador = 0
    
    for _ in range(100):
        
        t, a = runge_kutta(f, t, v, h)
        # v_f = v_func(v_values)
        # t, x, v = runge_kutta(v_f, t, v, h)
        x = pos_linear_model(x,v,a,h)
        v = vel_linear_model(v,a,h)
        print(x)
        t_values.append(t)
        v_values.append(v)
        x_values.append(x)
        print(t, v)

        contador = contador + 1

    # Plot solution
    plt.plot(t_values, x_values) # Gráfico da posição 
    plt.show()

        # Dados de entrada
    x_sp = np.linspace(0, 2*np.pi, 100)

    # Definindo os valores de y1, y2 e y3 em função de x
    y1 = np.sin(x_sp)
    y2 = np.cos(x_sp)
    y3 = np.exp(x_sp)

    # Criando splines cúbicos para cada conjunto de dados
    spline1 = CubicSpline(x_sp, y1)
    spline2 = CubicSpline(x_sp, y2)
    spline3 = CubicSpline(x_sp, y3)

    # Avaliando os splines cúbicos em x
    y1_interp = spline1(x_sp)
    y2_interp = spline2(x_sp)
    y3_interp = spline3(x_sp)

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

    

