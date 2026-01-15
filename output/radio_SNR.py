# ==============================================================================

# Este programa estima y grafica el radio del RSN como función del tiempo.

# ==============================================================================

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np

# ==============================================================================

# Constantes astrofísicas y fundamentales
AMU  = 1.660538782e-24
KB   = 1.380650400e-16
PC   = 3.085677588e+18
AU   = 1.495978707e+13
YR   = 3.155673600e+7
KYR  = 3.155673600e+10
MSUN = 1.988920000e+33
KPS  = 1.0e5

# Parámetros del problema
X_SNR, Y_SNR = 0, 0
time_step = 1  # Intervalo de tiempo entre archivos (en kyr)

# ==============================================================================

# Rango de outputs a graficar
out_min = 0
out_max = 50

# Guardar gráficas a disco en lugar de mostrarlas?
save_plots = True

# Estructura de los nombres de archivos
fname_template = "output_{:04d}.dat"

# Límites de las escalas de colores
# Poner None para usar el mínimo/máximo de los datos
dens_min = 1e-25
dens_max = 1e-23
dens_log = True
vel_min = 0
vel_max = 500
vel_log = False
pres_min = 1e-10
pres_max = 1e-9
pres_log = True
temp_min = 1e4
temp_max = 1e8
temp_log = True
mu0_amb = 1.3
mui_amb = 0.61

figsize=(21,10)
# figsize=(12,12)

tiempos = []
radios_expansion = []

# ==============================================================================

# Función para calcular el radio de expansión en cada iteración
#Sobre el eje x y el eje y suponiendo una velocidad puramente radial, determinamos la posicion
# donde se halle el maximo cambio de velocidad. A partir de esa posicion se determina el maximo respecto al centro
# de la explosion de supernova. Se toma un promedio de las medidas de ambas direcciones
def calcular_radio_expansion(x, y, vel_mag):
  # Calcular la diferencia en la magnitud de velocidad a lo largo del eje x
    dif_vel_x = np.diff(vel_mag[:, 0]) 
    #dif_vel_y = np.diff(vel_mag[125, :])
    # Encontrar el índice donde la diferencia es mayor que el umbral
    indice_discontinuidad_x = np.argmax(dif_vel_x)
    #indice_discontinuidad_y = np.argmax(dif_vel_y)
    # Calcular el radio de expansión
    radio_expansion = abs(x[indice_discontinuidad_x] - X_SNR)/PC
    #radio_expansion_y = abs(y[indice_discontinuidad_y] - Y_SNR)/PC

    #radio_expansion = (radio_expansion_x + radio_expansion_y)/2 
    #print(radio_expansion_x, radio_expansion_y, radio_expansion)    
    return radio_expansion

# Función de Sedov-Taylor para graficar de fondo    
def sedov_taylor(t):
#r en pc
# t en kyr
    return 1 + 5.07 * t**(2/5)
    
    
    

# Procesar cada archivo de salida y calcular el radio de expansión en cada iteración
for nout in range(out_min, out_max+1):

  # Generar nombre de archivo
  fname = fname_template.format(nout)

  # Leer datos del archivo
  t, = np.fromfile(fname, dtype="d", count=1)
  NEQ, NX, NY = np.fromfile(fname, dtype="i", count=3, offset=8)
  X1, X2, Y1, Y2 = np.fromfile(fname, dtype="d", count=4, offset=20)
  data = np.fromfile(fname, dtype="d", count=NEQ*NX*NY, offset=52).reshape((NEQ, NX, NY))

  # Extensión espacial de los eje x y y
  extent = [X1/PC, X2/PC, Y1/PC, Y2/PC]
  
  # Para mayor comodidad, separamos las primitivas individuales
  dens = data[0, :, :]
  velx = data[1, :, :]
  vely = data[2, :, :]
  pres = data[3, :, :]

  DX = (X2-X1)/NX      # Espaciamiento de la malla en X 
  DY = (Y2-Y1)/NY      # Espaciamiento de la malla en Y

  
  # Generar listas x e y con progresión aritmética
  x = X1 + np.arange(NX + 1) * DX
  y = Y1 + np.arange(NY + 1) * DY
  
  # ----------------------------------------------------------------------------

  # Calculamos la magnitud de la velocidad y la temperatura
  # Nótese que como estas operaciones se hacen con arreglos de numpy, se aplican
  # elemento por elemento en todo el arreglo
  vel_mag = np.sqrt(velx**2 + vely**2)
  temp = mui_amb*AMU/KB*pres/dens
  temp[temp < 1e4] = (mu0_amb*AMU/KB*pres/dens)[temp < 1e4]
    
  tiempo = nout # Supongo que el tiempo es simplemente el número de iteración
  tiempos.append(tiempo)
  radio = calcular_radio_expansion(x, y, vel_mag)
  radios_expansion.append(radio)

#print(tiempos, radios_expansion)
# Graficar el radio de expansión como función del tiempo
plt.plot(tiempos, radios_expansion, marker='o', label='Radios simulados')
t_sedov = np.linspace(0, 50, 100)  # Ajusta según el rango de tus datos
r_sedov = sedov_taylor(t_sedov)
plt.plot(t_sedov, r_sedov, color='red', linestyle='-', label='Sedov-Taylor')

plt.xlabel('Tiempo (kyr)')
plt.ylabel('Radio de expansión (pc)')
plt.title('Radio de expansión del remanente de supernova')
plt.grid(True)
plt.legend()
plt.savefig('Sedov-Taylor.png')
plt.show()

    
    
    
    
    
    

  
