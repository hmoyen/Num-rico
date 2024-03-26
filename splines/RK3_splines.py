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

interval = 4*np.pi
division = 100
timestep = interval/(division) 
x_sp = np.linspace(0, interval, division)

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
        # print("H:", h)
        # print("t:", t)
        # print("self.x", self.x[idx])
        # print("Length",len(self.x), "Index", idx)
        # print("Parameters",self.a[idx], self.b[idx], self.c[idx], self.d[idx])
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
            # print("xzinho", self.x[idx])
            # print("Hzinho", h)
            y = self.a[idx] + self.b[idx] * h + self.c[idx] * h ** 2 + self.d[idx] * h ** 3
            y_eval.append(y)
        return y_eval

def runge_kutta(f, x_0, y_0, h):
    """Three step Runge-Kutta method (RK3)
    Solves first order ODEs
    """
    k_0 = f(x_0, y_0)
    k_1 = f(x_0 + h/2, y_0 + h/2 * k_0)
    k_2 = f(x_0 + h, y_0 - h*k_0 + 2*h*k_1)

    k = 1/6 * (k_0 + 4.0*k_1 + k_2)

    x_1 = x_0 + h
    # y_1 = y_0 + h * k

    return x_1, k


def f1(t, y):
    """Example first order ordinary differential equation (ODE)"""
    
    # Dados de entrada
    x = x_sp

    y = np.cos(x) # Function

    # Criando splines cúbicos para cada conjunto de dados
    spline = CubicSpline(x, y)
    # y3_interp = spline3(x)

    # Avaliando os splines cúbicos em x
    y_func = spline.calc_func(t,y)
    # print("Function;",y3_func)

    return y_func

def f2(t, y):
    """Example first order ordinary differential equation (ODE)"""
    
    # Dados de entrada
    x = x_sp

    y = np.sin(x) # Function

    # Criando splines cúbicos para cada conjunto de dados
    spline = CubicSpline(x, y)
    # y3_interp = spline3(x)

    # Avaliando os splines cúbicos em x
    y_func = spline.calc_func(t,y)
    # print("Function;",y3_func)

    return y_func

def f3(t, y):
    """Example first order ordinary differential equation (ODE)"""
    
    # Dados de entrada
    x = x_sp

    y = np.exp(x) # Function

    # Criando splines cúbicos para cada conjunto de dados
    spline = CubicSpline(x, y)
    # y3_interp = spline3(x)

    # Avaliando os splines cúbicos em x
    y_func = spline.calc_func(t,y)
    # print("Function;",y3_func)

    return y_func

def v_func(values):

    x = x_sp
    splineV = CubicSpline(x,values)
    # func= splineV.calc_func(t,y)
    return splineV.calc_func()


def pos_linear_model(p,v,a,delta):

    return (p+v*delta+(a*delta**2)*0.5)

def vel_linear_model(v,a,delta):

    return (v+a*delta)

def exact_solution(index, derivative, x):

    #Calculating exact solutions for f1,f2 and f3

    if index == 1 and derivative == 1:
        return (-np.cos(x))
    elif index == 1 and derivative ==2:
        return (np.sin(x))
    elif index == 2 and derivative == 1:
        return (-np.sin(x))
    elif index == 2 and derivative ==2:
        return (-np.cos(x))
    elif index == 3 and derivative ==1:
        return (np.exp(x))
    elif index == 3 and derivative ==2:
        return (np.exp(x))

def calculate_error(approximate, exact):
    """Calculate error between approximate and exact solutions"""
    approximate_array = np.array(approximate)
    exact_array = np.array(exact)
    return np.abs(approximate_array - exact_array)

def plot_variable(t_values, exact_values, approx_values, variable_name):
    # Plot the exact and approximate solutions
    plt.plot(t_values, exact_values, label=f'Exact Solution for {variable_name}', color='green')
    plt.plot(t_values, approx_values, label=f'Approximate Solution for {variable_name}', color='red', linestyle=(0,(1,1,3,1)))
    plt.xlabel('time t (in seconds)')
    plt.ylabel(f'{variable_name} (in units)')
    plt.title(f'Numerical Approximation of State Variable {variable_name}')
    plt.legend()
    plt.show()

    # Calculate the error
    error = calculate_error(exact_values, approx_values)

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
    print(latex_table)

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
    # Initial values
    t_0 = 0.0

    v_0_1 = 0.0
    x_0_1 = -1

    v_0_2 = -1
    x_0_2 = 0.0

    v_0_3 = 1
    x_0_3 = 1

    # Step length 15Hz aproximadamente
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
        
        t, a1 = runge_kutta(f1, t_last, v1, h)
        t, a2 = runge_kutta(f2, t_last, v2, h)
        t, a3 = runge_kutta(f3, t_last, v3, h)

        x1 = pos_linear_model(x1,v1,a1,h)
        v1 = vel_linear_model(v1,a1,h)

        x2 = pos_linear_model(x2,v2,a2,h)
        v2 = vel_linear_model(v2,a2,h)

        x3 = pos_linear_model(x3,v3,a3,h)
        v3 = vel_linear_model(v3,a3,h)

        v1_exact = exact_solution(1,2,t)
        x1_exact = exact_solution(1,1,t)
        v1_exact_values.append(v1_exact)
        x1_exact_values.append(x1_exact)

        v2_exact = exact_solution(2,2,t)
        x2_exact = exact_solution(2,1,t)
        v2_exact_values.append(v2_exact)
        x2_exact_values.append(x2_exact)

        v3_exact = exact_solution(3,2,t)
        x3_exact = exact_solution(3,1,t)
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

    # Plot and analyze each variable
    variables = ['x1', 'x2', 'x3', 'v1', 'v2', 'v3']
    for variable in variables:
        exact_values = globals()[f'{variable}_exact_values']
        approx_values = globals()[f'{variable}_values']
        plot_variable(t_values, exact_values, approx_values, variable)
    

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

    

