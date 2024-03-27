import math


  

      


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

    
    