#! /usr/bin/env python

'''
Este es un script que integra el espectro del Sol
e integra la funcion de Planck
'''


import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as cons
from astropy import units as u

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
