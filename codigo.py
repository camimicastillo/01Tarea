#! /usr/bin/env python

'''
Este es un script que integra el espectro del Sol
e integra la funcion de Planck
'''


import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as cons
from astropy import units as u
from scipy import integrate as integrate

#Lectura del archivo que contiene el espectro del Sol
D = np.loadtxt("sun_AM0.dat")
x = D[:,0]
y = D[:,1]

#Conversion de los datos a Flujo [cgs] y Longitud de onda [um]
xdato= x*u.nm
x_um= xdato.to('um')
ydato= y*u.W*(u.m**-2)*(u.nm**-1)
y_cgs= ydato.to('erg/(s cm2 um)')



#Plotear Flujo vs Longitud de onda
plt.clf()
plt.plot(x_um, y_cgs)
plt.xlim(0,7)
plt.xlabel('Longitud de onda [$ \mu m $]')
plt.ylabel('Flujo [$erg / s * cm^2 * \mu m$]')
plt.title('Grafico de Flujo versus Longitud de onda')
plt.savefig('Figura.png', bbox_inches='tight')
plt.show()


#Integracion del espectro en longitud de onda usando el metodo del trapecio
#"sum" corresponde al valor de la integral del espectro (Constante solar)
#y "Ls" al valor de la Luminosidad total del sol
Ks=0
for i in range(1696):
    valor= ((y_cgs[i]+y_cgs[i+1])*(x_um[i+1] - x_um[i])) / 2
    Ks += valor
print ('Integral del espectro del Sol: CONSTANTE SOLAR')
print (Ks)

r= cons.au #Distancia entre el Sol y la Tierra en m
rcm= r.to('cm')
Ls= 4*np.pi*(rcm*rcm)*Ks
print ('Luminosidad total del Sol')
print (Ls)

#Integracion de la funcion de Planck
#A partir de la ecuacion dada en la tarea, se realiza el cambio de
#variable y= arctan(x). Ahora calcularemos esta integral con el metodo de
#Simpson, pero usando una formula vista en clases que ocupa el metodo del
#Trapecio para llegar al metodo de Simpson

#Definir vectores x e y(x)
xin=0.01
xfi=3.14 / 2
x_pl= np.arange(xin, xfi, 0.01)
y_pl= (np.tan(x_pl)**3)/((np.cos(x_pl)**2)*(np.exp(np.tan(x_pl))-1))
elementos=len(x_pl)

#Iteracion Trapecio para n puntos (n=elementos del vector x_pl)
Tn=0
for j in range(elementos-1):
    ar= ((y_pl[j]+y_pl[j+1])*(x_pl[j+1] - x_pl[j])) / 2
    Tn += ar
print('Valor integral con Trapecio con n puntos')
print (Tn)

#Definir nuevos vectores x e y(x) con el doble de elementos
xinn=0.01
xfin=3.14 / 2
x_pln= np.arange(xin, xfi, 0.005)
y_pln= (np.tan(x_pln)**3)/((np.cos(x_pln)**2)*(np.exp(np.tan(x_pln))-1))
elementoss=len(x_pln)

#Iteracion Trapecio para 2n puntos
T_2n=0
for k in range(elementoss-1):
    inte= ((y_pln[k]+y_pln[k+1])*(x_pln[k+1] - x_pln[k])) / 2
    T_2n += inte
print ('Valor integral con Trapecio 2n puntos')
print (T_2n)

#Aplicar formula vista en clases para Simpson
S= (4*T_2n / 3) - (Tn)/3
print ('Valor integral con metodo de Simpson (formula vista en clases que utiliza Trapecio)')
print (S)

#Usamos el metodo del punto medio para calcular la integral cuando la
#funcion diverge, en este caso en 0
yin_medio= (np.tan(xin/2)**3)/((np.cos(xin/2)**2)*(np.exp(np.tan(xin/2))-1))
area_0=xin*yin_medio
print ('Valor de integral con metodo punto medio en x=0, debido a que es una divergencia')
print (area_0)

#Valor final de la integral
Sfinal= S+area_0
print ('Valor de la integral total: Simpson + Punto medio')
print (Sfinal)

#Valor de la integral de la funcion de Planck: que es constantes * Integral
t= 5778 * u.K
cplanck= ((2*np.pi*cons.h)* (((cons.k_B*t)/(cons.h))**4))/((cons.c)**2)
planck= cplanck * Sfinal
print ('Valor de integral de funcion de Planck: integral calculada * constantes . FLUJO DE ENERGIA ')
print(planck)

#Conversion de la constante solar a J/m^2 * s
KsJ=Ks.to('J / (m2 s)')

#Calcular el radio del Sol con Formula que relaciona Constante Solar, el flujo
#de energia y la unidad astronomica
rSOL= ((KsJ*(cons.au**2))/planck)**0.5
print ('Valor del radio del Sol')
print (rSOL)


#Recalcular las integrales anteriores con los metodos predeterminados de python
#Calculo integral del espectro con metodo de trapecio predeterminado
esp_trapz= np.trapz(y_cgs, x_um)
print ('Valor integral del espectro con metodo trapecio predeterminado')
print (esp_trapz)

#Calculo de integral de funcion de Planck con metodo de trapecio predeterminado
pln_trapz= np.trapz(y_pl, x_pl)
plntrapz= pln_trapz * cplanck
print ('Valor integral de Planck con metodo trapecio predeterminado')
print (plntrapz)

#Calculo de la integral de funcion de Planck con quad
fu= lambda x: (x**3)/(np.exp(x) - 1)
pln_quad= integrate.quad( fu , 0, np.inf)
print ('Valor de la integral desde 0 a infinito de (x^3 / exp(x) - 1)')
print (pln_quad)

pln_quad1= pln_quad[0]
plnquad= pln_quad1 * cplanck
print ('Valor de la integral de la funcion de Planck con metodo predeterminado quad')
print (plnquad)

#Calculo del tiempo de ejecucion
