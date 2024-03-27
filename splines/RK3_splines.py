import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import matplotlib.pyplot as plt
from math import exp, sqrt

import pandas as pd 

__date__ = datetime(2024, 3, 24) # or version string or something
__author__ = "Helena Moyen"

global timestep
global x_sp
global division


# Conversions

# 5,3 Hz --> division = 50, interval = 3 * pi
# 16Hz --> division = 100, interval = 2 * pi
# 32Hz --> division = 200, interval = 2 * pi
# 64Hz --> division = 400, interval = 2 * pi : approximate frequency of Intel Realsense T265 accelerometer


interval =  np.pi
division = 512
timestep = interval/(division) 
x_sp = np.linspace(0, interval, division)

class CubicSpline:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.n = len(x)
        self.a, self.b, self.c, self.d = self.compute_coefficients()

    def compute_coefficients(self):
        # Calcula as diferenças entre os pontos adjacentes em x
        h = [self.x[i+1] - self.x[i] for i in range(self.n - 1)]
        
        # Inicializa a lista alpha com zeros
        alpha = [0] * self.n
        
        # Calcula os valores de alpha conforme a fórmula do pseudocódigo (Numerical Analysis)
        for i in range(1, self.n - 1):
            alpha[i] = 3 * ((self.y[i + 1] - self.y[i]) / h[i] - (self.y[i] - self.y[i - 1]) / h[i - 1])
        
        # Inicializa as listas l, mu e z
        l = [1] + [0] * (self.n - 1)
        mu = [0] * self.n
        z = [0] * self.n
        
        # Calcula os valores de l, mu e z 
        for i in range(1, self.n - 1):
            l[i] = 2 * (self.x[i + 1] - self.x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / l[i]
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
        
        # Define os valores finais de l e z
        l[self.n - 1] = 1
        z[self.n - 1] = 0
        
        # Inicializa as listas c, b e d
        c = [0] * self.n
        b = [0] * self.n
        d = [0] * self.n
        
        # Calcula os valores de c, b e d 
        for j in range(self.n - 2, -1, -1):
            c[j] = z[j] - mu[j] * c[j + 1]
            b[j] = (self.y[j + 1] - self.y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
            d[j] = (c[j + 1] - c[j]) / (3 * h[j])
        
        # Define os valores de a como os valores de y, exceto o último ponto
        a = self.y[:-1]
        
        # Retorna os coeficientes calculados
        return a, b, c, d

    
    def calc_func(self,t,d):

        # Derivada d nao aparece em funcoes que dependem apenas do tempo
        idx = 0
        for i in range(self.n - 1):
            if self.x[i] <= t <= self.x[i + 1]:
                idx = i
                break
        h = t - self.x[idx]
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
            y = self.a[idx] + self.b[idx] * h + self.c[idx] * h ** 2 + self.d[idx] * h ** 3
            y_eval.append(y)
        return y_eval

def runge_kutta(f, x_0, y_0, h):
    """Three step Runge-Kutta method (RK3)
    Solves first order ODEs
    """
    k_0 = f.calc_func(x_0, y_0)
    k_1 = f.calc_func(x_0 + h/2, y_0 + h/2 * k_0)
    k_2 = f.calc_func(x_0 + h, y_0 - h*k_0 + 2*h*k_1)

    k = 1/6 * (k_0 + 4.0*k_1 + k_2)

    x_1 = x_0 + h

    return x_1, k


def pos_linear_model(p,v,a,delta):

    # p -> posição anterior
    # v --> velocidade estimada com runge kutta
    # a --> medida ponderada das acelerações dado por runge kutta
    # delta --> tempo

    return (p+v*delta+(a*delta**2)*0.5)

def vel_linear_model(v,a,delta):

    return (v+a*delta)

def exact_solution(index, x):

    #Calculating exact solutions for f1,f2 and f3

    if index == 'x1':
        return (-np.cos(x))
    elif index == 'v1':
        return (np.sin(x))
    elif index == 'x2':
        return (-np.sin(x))
    elif index == 'v2':
        return (-np.cos(x))
    elif index == 'x3':
        return (np.exp(x))
    elif index == 'v3':
        return (np.exp(x))

def calculate_error(approximate, exact):
    """Calculate error between approximate and exact solutions"""
    approximate_array = np.array(approximate)
    exact_array = np.array(exact)
    return np.abs(approximate_array - exact_array)

def plot_variable(t_values, exact_values, approx_values, variable_name, plot_func, t_func):
    # Plot the exact and approximate solutions
    plt.plot(t_func, plot_func, label=f'Exact Solution for {variable_name}', color='green')
    plt.plot(t_values, approx_values, label=f'Approximate Solution for {variable_name}', color='red', linestyle=(0,(1,1,3,1)))
    plt.xlabel('time t (in seconds)')
    plt.ylabel(f'{variable_name} (in units)')
    plt.title(f'Numerical Approximation of State Variable {variable_name}')
    plt.legend()
    plt.show()

    # Calculate the error
    error = calculate_error(exact_values, approx_values)
    print("Approx solution", approx_values[int(division/2)])
    print("Time solution", t_values[int(division/2)])
    print("Error", error[int(division/2)])
    print("h", timestep)
    print("n", division)

    # Create DataFrame for error
    error_df = pd.DataFrame({
        't': t_values,
        'Exact Solution': exact_values,
        'Approximate Solution (Runge-Kutta)': approx_values,
        'e': error
    })

    # Slice the DataFrame to select only 10 rows in the middle of the interval
    start_index = len(t_values) // 2 - 5
    end_index = len(t_values) // 2 + 5
    error_df_slice = error_df.iloc[start_index:end_index]

    # Convert DataFrame to LaTeX table format
    latex_table = error_df_slice.to_latex(index=False)
    print(f"{variable_name} global error")
    # print(latex_table)

    # Display the sliced table in a figure
    plt.figure(figsize=(10, 6))
    table = plt.table(cellText=error_df_slice.values,
                      colLabels=error_df_slice.columns,
                      cellLoc='center',
                      loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)
    plt.axis('off')
    plt.show()


if __name__=="__main__":

    # Function

    y1 = np.cos(x_sp)
    y2 = np.sin(x_sp)
    y3 = np.exp(x_sp) 
    

    # Criando splines cúbicos para cada conjunto de dados
    spline1 = CubicSpline(x_sp, y1)
    spline2 = CubicSpline(x_sp, y2)
    spline3 = CubicSpline(x_sp, y3)


    # Initial values
    t_0 = 0.0

    v_0_1 = 0.0
    x_0_1 = -1

    v_0_2 = -1
    x_0_2 = 0.0

    v_0_3 = 1
    x_0_3 = 1


    h = timestep
    

    t_values = [t_0]

    v1_values = [v_0_1]
    x1_values = [x_0_1]

    v2_values = [v_0_2]
    x2_values = [x_0_2]

    v3_values = [v_0_3]
    x3_values = [x_0_3]

    # Calculate solution
    t = t_0

    v1 = v_0_1
    x1 = x_0_1

    v2= v_0_2
    x2 = x_0_2 

    v3 = v_0_3
    x3 = x_0_3

    #Exact solution
    v1_exact_values = [v_0_1]
    x1_exact_values = [x_0_1]

    v2_exact_values = [v_0_2]
    x2_exact_values = [x_0_2]

    v3_exact_values = [v_0_3]
    x3_exact_values = [x_0_3]

    v1_exact = v_0_1
    x1_exact = x_0_1

    v2_exact = v_0_2
    x2_exact = x_0_2

    v3_exact = v_0_3
    x3_exact = x_0_3

    contador = 0
    
    for _ in range(division):
        
        t_last = t
        
        t, a1 = runge_kutta(spline1, t_last, v1, h)
        t, a2 = runge_kutta(spline2, t_last, v2, h)
        t, a3 = runge_kutta(spline3, t_last, v3, h)

        x1 = pos_linear_model(x1,v1,a1,h)
        v1 = vel_linear_model(v1,a1,h)

        x2 = pos_linear_model(x2,v2,a2,h)
        v2 = vel_linear_model(v2,a2,h)

        x3 = pos_linear_model(x3,v3,a3,h)
        v3 = vel_linear_model(v3,a3,h)

        v1_exact = exact_solution('v1',t)
        x1_exact = exact_solution('x1',t)
        v1_exact_values.append(v1_exact)
        x1_exact_values.append(x1_exact)

        v2_exact = exact_solution('v2',t)
        x2_exact = exact_solution('x2',t)
        v2_exact_values.append(v2_exact)
        x2_exact_values.append(x2_exact)

        v3_exact = exact_solution('v3',t)
        x3_exact = exact_solution('x3',t)
        v3_exact_values.append(v3_exact)
        x3_exact_values.append(x3_exact)
        # print(x)
        t_values.append(t)
        v1_values.append(v1)
        x1_values.append(x1)

        v2_values.append(v2)
        x2_values.append(x2)

        v3_values.append(v3)
        x3_values.append(x3)
        # print(t, v)

        contador = contador + 1

    t_complete = np.linspace(0, interval, 200)
    x1_complete = []
    x2_complete = []
    x3_complete = []
    v1_complete = []
    v2_complete = []
    v3_complete = []
    
    for i in t_complete:
        x1 = exact_solution('x1',i)
        x2 = exact_solution('x2',i)
        x3 = exact_solution('x3',i)
        v1 = exact_solution('v1',i)
        v2 = exact_solution('v2',i)
        v3 = exact_solution('v3',i)
        v1_complete.append(v1)
        x1_complete.append(x1)

        v2_complete.append(v2)
        x2_complete.append(x2)

        v3_complete.append(v3)
        x3_complete.append(x3)

    #Plot and analyze each variable
    variables = ['x1', 'x2', 'x3', 'v1', 'v2', 'v3']
    for variable in variables:
        exact_values = globals()[f'{variable}_exact_values']
        approx_values = globals()[f'{variable}_values']
        func = globals()[f'{variable}_complete']
        plot_variable(t_values, exact_values, approx_values, variable, func, t_complete)
    

    # Configuração da figura 3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plotando os pontos dos splines cúbicos em 3D com linhas contínuas
    ax.plot(x1_values, x2_values, x3_values, color='b', label='Linha da Trajetória Aproximada')
    ax.plot(x1_complete, x2_complete, x3_complete, color='r', label='Linha da Trajetória Real')

    # Configuração dos eixos
    ax.set_xlabel('X1')
    ax.set_ylabel('X2')
    ax.set_zlabel('X3')

    # Adicionando legenda
    plt.legend()

    # Exibição do gráfico
    plt.title('Posição Linear')
    plt.show()

    

