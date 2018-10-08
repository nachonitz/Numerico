#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TP1 ANALISIS NUMERICO
Problemas de busqueda de raices
Comparacion de metodos

2 Cuatrimestre 2018

Por:
Ignacio

"""

#Imports
from scipy.optimize import brentq
import timeit #Para calcular tiempo de corrida
import numpy as np #Manejo de arrays
import matplotlib.pylab as plt #Rutinas gráficas

#Funciones f1, f2, f3 y derivadas
def f1(x): 
    return x**2-2
def df1(x): 
    return 2*x
def ddf1(x): 
    return 2

def f2(x): 
    return 0
def df2(x): 
    return 0
def ddf2(x): 
    return 0

def f3(t): 
    return 0
def df3(t): 
    return 0
def ddf3(t): 
    return 0

#Funciones busqueda de raices
def bisec(f, a, b, a_tol, n_max):
    """
    Devolver (x0, delta), raiz y cota de error por metodo de la biseccion
    Datos deben cumplir f(a)*f(b) > 0
    """
    x = a+(b-a)/2    #mejor que (a+b)/2 segun Burden
    delta = (b-a)/2
    
    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'a', 'b'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x, a, b))
    
    for i in range(n_max):
        if f(a) * f(x) > 0:
            a = x
        else:
            b = x
        x_old = x
        x = a+(b-a)/2 #(a+b)/2
        delta = np.abs(x - x_old)
        
        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x, a, b))
        
        if delta <= a_tol: #Hubo convergencia
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1
    
    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x, delta, i+1

def secante(f, x0, x1, a_tol, n_max):
    """
    Devolver (x, delta), raiz y cota de error por metodo de la secante
    """
    delta = 0

    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'x_-1', 'delta'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x1, x0, delta))

    for i in range(n_max):
        x = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        delta = np.abs(x - x1)
        x0 = x1
        x1 = x
        
        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x1, x0, delta))
        
        #Chequear convergencia
        if delta <= a_tol: #Hubo convergencia
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    #Si todavia no salio es que no hubo convergencia:
    raise ValueError('No hubo convergencia')
    return x. delta, i+1

#Intervalo para buscar raiz
a = 0.0
b = 2.0

#Parametros para el algoritmo
a_tol1 = 1e-5
a_tol2 = 1e-13
n_max = 100

#Grafica de las funciones
#Ver https://matplotlib.org
xx = np.linspace(a, b, 256+1)
yy = f1(xx)
nombre = 'f1'
plt.figure(figsize=(10,7))
plt.plot(xx, yy, lw=2)
#plt.legend(loc=best)
plt.xlabel('x')
plt.ylabel(nombre +'(x)')
plt.title('Funcion '+ nombre)
plt.grid(True)
plt.savefig(nombre + '.png')
plt.show()

print('----------------')
print('Metodo biseccion')
print('----------------')
print('')
print('Funcion f1, a_tol = '+str(a_tol1))
r, delta, n_iter = bisec(f1, a, b, a_tol1, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')
print('Funcion f1, a_tol = '+str(a_tol2))
r, delta, n_iter = bisec(f1, a, b, a_tol2, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')

print('----------------')
print('Metodo secante')
print('----------------')
print('')
print('Funcion f1, a_tol = '+str(a_tol1))
r, delta, n_iter = secante(f1, a, b, a_tol1, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')
print('Funcion f1, a_tol = '+str(a_tol2))
r, delta, n_iter = secante(f1, a, b, a_tol2, n_max)
print('raiz = ' +str(r))
print('delta= ' +str(delta))
print('n_ite= ' +str(n_iter))
print('')

print('----------------')
print('Metodo brent')
print('----------------')
print('')
print('Funcion f1, a_tol por defecto para la funcion')
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html
r, results = brentq(f1, a, b, full_output=True)
print('raiz = ' +str(r))
print('Resultados: ')
print(results)


