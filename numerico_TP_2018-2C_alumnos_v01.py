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
import math

#Funciones f1, f2, f3 y derivadas
def f1(x): 
    return x**2-2
def df1(x): 
    return 2*x
def ddf1(x): 
    return 2

def f2(x): 
    return x**5 - 6.6*x**4 + 5.12*x**3 + 21.312*x**2 - 38.016*x + 17.28
def df2(x): 
    return 5*x**4 - 26.4*x**3 + 15.36*x**2 +42.624*x - 38.016
def ddf2(x): 
    return 20*x**3 - 79.2*x**2 + 30.72*x + 42.624


def f3(t): 
    return (t - 1.5) * math.exp(-4 * (t - 1.5)**2)
def df3(t): 
    return (-8*t + 12.0)*(t - 1.5) * math.exp(-4*(t - 1.5)**2) + math.exp(-4*(t - 1.5)**2)
def ddf3(t): 
    return (-24*t + (t - 1.5)*(8*t - 12.0)**2) * math.exp(-4* (t - 1.5)**2)

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
    print('No hubo convergencia')
    return x, delta, i+1

def secante(f, x0, x1, a_tol, n_max):
    """
    Devolver (x, delta), raiz y cota de error por metodo de la secante
    """
    delta = 0

    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'x_-1', 'delta'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x1, x0, delta))

    f_x0 = f(x0)
    f_x1 = f(x1)

    for i in range(n_max):
        x = x1 - f_x1*(x1-x0)/(f_x1-f_x0)
        delta = np.abs(x - x1)
        x0 = x1
        x1 = x

        f_x0 = f_x1
        f_x1 = f(x1)
        
        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x1, x0, delta))
        
        #Chequear convergencia
        if delta <= a_tol: #Hubo convergencia
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    #Si todavia no salio es que no hubo convergencia:
    print('No hubo convergencia')
    return x, delta, i+1

def newton_raphson(f, df, a_tol, n_max):

    delta = 0
    x = 1.0


    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'x_-1', 'delta'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x, 0, delta))

    for i in range(n_max):
        x_old = x
        x = x_old - f(x_old)/df(x_old)
        delta = np.abs(x - x_old)

        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x, x_old, delta))

        if delta <= a_tol:
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    print('No hubo convergencia')
    return x, delta, i+1

def newton_raphson_modificado(f, df, ddf, a_tol, n_max):

    delta = 0
    x = 1.0

    print('{0:^4} {1:^17} {2:^17} {3:^17}'.format('i', 'x', 'x_-1', 'delta'))
    print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(0, x, 0, delta))

    for i in range(n_max):
        x_old = x
        x = x_old - (f(x_old)*df(x_old))/(df(x_old)*df(x_old) - f(x_old)*ddf(x_old))
        delta = np.abs(x - x_old)

        print('{0:4} {1: .14f} {2: .14f} {3: .14f}'.format(i+1, x, x_old, delta))

        if delta <= a_tol:
            print('Hubo convergencia, n_iter = ' + str(i+1))
            return x, delta, i+1

    print('No hubo convergencia')
    return x, delta, i+1

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

funciones = ((f1,df1,ddf1),(f2,df2,ddf2),(f3,df3,ddf3))
metodos = ("biseccion","secante","newton_raphson","newton_raphson_modificado","brent")
tolerancias = (a_tol1,a_tol2)
for funcion in funciones:
    if funcion == (f1,df1,ddf1):
        funcion_actual = "f1"
    if funcion == (f2,df2,ddf2):
        funcion_actual = "f2"
    if funcion == (f3,df3,ddf3):
        funcion_actual = "f3"
    for metodo in metodos:
        
        
            
        if metodo != "brent":
            print('----------------')
            print('Metodo {}'.format(metodo))
            print('----------------')
            print('')
            print('Funcion {}, a_tol = '.format(funcion_actual)+str(a_tol1))
            if metodo == "biseccion":
                r, delta, n_iter = bisec(funcion[0], a, b, a_tol1, n_max)
                
            if metodo == "secante":
                r, delta, n_iter = secante(funcion[0], a, b, a_tol1, n_max)
                
            if metodo ==  "newton_raphson":
                r, delta, n_iter = newton_raphson(funcion[0], funcion[1], a_tol1, n_max)
                
            if metodo == "newton_raphson_modificado":
                r, delta, n_iter = newton_raphson_modificado(funcion[0], funcion[1], funcion[2], a_tol1, n_max)
            print('raiz = ' +str(r))
            print('delta= ' +str(delta))
            print('n_ite= ' +str(n_iter))
            print('')
            print('Funcion {}, a_tol = '.format(funcion_actual)+str(a_tol2))
        
            if metodo == "biseccion":
                metodo_actual = "biseccion"
                r, delta, n_iter = bisec(funcion[0], a, b, a_tol2, n_max)
                
            if metodo == "secante":
                r, delta, n_iter = secante(funcion[0], a, b, a_tol2, n_max)
                
            if metodo ==  "newton_raphson":
                r, delta, n_iter = newton_raphson(funcion[0], funcion[1], a_tol2, n_max)
                
            if metodo == "newton_raphson_modificado":
                r, delta, n_iter = newton_raphson_modificado(funcion[0], funcion[1], funcion[2], a_tol2, n_max)
            
            print('raiz = ' +str(r))
            print('delta= ' +str(delta))
            print('n_ite= ' +str(n_iter))
            print('')

        if metodo == "brent":
            print('----------------')
            print('Metodo brent')
            print('----------------')
            print('')
            print('Funcion {}, a_tol por defecto para la funcion'.format(funcion_actual))
            #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html
            r, results = brentq(funcion[0], a, b, full_output=True)
            print('raiz = ' +str(r))
            print('Resultados: ')
            print(results)
