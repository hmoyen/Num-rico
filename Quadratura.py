import math
import sympy as sp


  

      


f_t0 = input("Digite valor inicial: ")
t = int(input("Digite o tamanho do intervalo: "))
n = input("Digite a quantidade de medicoes: ")
a = 0
b = float(t) 

n = int(n)
val = []
f_t0 = float(f_t0)
val.append(f_t0)
for k in range(0, n):
    val.append(float(input("Digite o proximo valor: ")))
 

h = (b-a)/n 
soma = 0
for k in range(1, n):
    soma += (val[k])
    
        
soma = soma * 2
soma = soma +  val[a] + val[n]
    
res = soma * h * 0.5

print (res)




# Definindo os símbolos
t, t0, n, tau = sp.symbols('t t0 n tau')

# Definindo a função de segunda derivada omega''(tau)
omega_double_prime = sp.Function('omega')(tau).diff(tau, 2)

# Definindo a fórmula do erro
erro = sp.Abs((t - t0)**3 / (12 * n**2) * sp.Max(omega_double_prime, (tau, t, t0)))

# Mostrando a fórmula
print("Fórmula do erro:")
sp.pprint(erro)
    

# Definindo o símbolo e a função
x = sp.symbols('x')


# Encontrando a segunda derivada
f_double_prime = sp.diff(f, x, 2)

# Mostrando a segunda derivada
print("Segunda derivada de f(x):")
print(f_double_prime)
    