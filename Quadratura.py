import math
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import matplotlib.pyplot as plt
from math import exp, sqrt

import pandas as pd 


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
    
      
class max_seg_deriv:
    # Function

    y1 = np.cos(x_sp)
    y2 = np.sin(x_sp)
    y3 = np.exp(x_sp) 
    

    # Criando splines cúbicos para cada conjunto de dados
    spline1 = CubicSpline(x_sp, y1)
    spline2 = CubicSpline(x_sp, y2)
    spline3 = CubicSpline(x_sp, y3)
    

      

if __name__== "__main__":
# f_t0 = input("Digite valor inicial: ")
# t = int(input("Digite o tamanho do intervalo: "))
# n = input("Digite a quantidade de medicoes: ")
# a = 0
# b = float(t) 

    f_t0 = 1
    n = 8 #quantidade de medições
    a = 0
    b = float(n) 

# n = int(n)
# val = []
# f_t0 = float(f_t0)
# val.append(f_t0)
# for k in range(0, n):
#     val.append(float(input("Digite o proximo valor: ")))

n = int(n)
val = []
f_t0 = float(f_t0)
val.append(f_t0)
val.append(2.65)
val.append(3.61)
val.append(4.36)
val.append(5.0)
val.append(5.57)
val.append(6.08)
val.append(6.56)
val.append(7.0)
# val.append(1)
# val.append(f_t0)
# val.append(1)
# val.append(f_t0)
# val.append(1)
# val.append(f_t0)
# val.append(1)



 

h = (b-a)/n 
soma = 0
for k in range(1, n):
    soma += (val[k])
    
        
soma = soma * 2
soma = soma +  val[a] + val[n]
    
res = soma * h * 0.5

print (res)



erro_max = ((b-a)^3/12*(n)^2) * max_seg_deriv


    