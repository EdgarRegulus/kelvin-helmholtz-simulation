# ==============================================================================

# Este programa grafica la salida del código que resuelve las ecuaciones
# de Euler en 2D, en formato binario.

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

# ==============================================================================

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

  # ----------------------------------------------------------------------------

  # Calculamos la magnitud de la velocidad y la temperatura
  # Nótese que como estas operaciones se hacen con arreglos de numpy, se aplican
  # elemento por elemento en todo el arreglo
  vel_mag = np.sqrt(velx**2 + vely**2)
  temp = mui_amb*AMU/KB*pres/dens
  temp[temp < 1e4] = (mu0_amb*AMU/KB*pres/dens)[temp < 1e4]

  plt.figure(figsize=figsize)

  # Densidad
  plt.subplot(2,2,1)
  if dens_log:
    im = plt.imshow(dens.T, cmap="jet", norm=LogNorm(vmin=dens_min, vmax=dens_max), origin='lower', extent=extent)
  else:
    im = plt.imshow(dens.T, cmap="jet", vmin=dens_min, vmax=dens_max, origin='lower', extent=extent)
  plt.title("Densidad [g/cm³]")
  plt.xlabel("x [pc]")
  plt.ylabel("y [pc]")
  ax_divider = make_axes_locatable(plt.gca())
  cax = ax_divider.append_axes("right", size="3%", pad="2%")
  plt.colorbar(im, cax=cax)

  # Velocidad
  plt.subplot(2,2,2)
  if vel_log:
    im = plt.imshow(vel_mag.T/KPS, cmap="rainbow", norm=LogNorm(vmin=vel_min, vmax=vel_max), origin='lower', extent=extent)
  else:
    im = plt.imshow(vel_mag.T/KPS, cmap="rainbow", vmin=vel_min, vmax=vel_max, origin='lower', extent=extent)
  plt.title("Velocidad (magnitud) [km/s]")
  plt.xlabel("x [pc]")
  plt.ylabel("y [pc]")  
  ax_divider = make_axes_locatable(plt.gca())
  cax = ax_divider.append_axes("right", size="3%", pad="2%")
  plt.colorbar(im, cax=cax)

  # Presión
  plt.subplot(2,2,3)
  if pres_log:
    im = plt.imshow(pres.T, cmap="viridis", norm=LogNorm(vmin=pres_min, vmax=pres_max), origin='lower', extent=extent)
  else:
    im = plt.imshow(pres.T, cmap="viridis", vmin=pres_min, vmax=pres_max, origin='lower', extent=extent)
  plt.title("Presión [erg/cm^3]")
  plt.xlabel("x [pc]")
  plt.ylabel("y [pc]")  
  ax_divider = make_axes_locatable(plt.gca())
  cax = ax_divider.append_axes("right", size="3%", pad="2%")
  plt.colorbar(im, cax=cax)

  # Temperatura
  plt.subplot(2,2,4)
  if temp_log:
    im = plt.imshow(temp.T, cmap="plasma", norm=LogNorm(vmin=temp_min, vmax=temp_max), origin='lower', extent=extent)
  else:
    im = plt.imshow(temp.T, cmap="plasma", vmin=temp_min, vmax=temp_max, origin='lower', extent=extent)
  plt.title("Temperatura [K]")
  plt.xlabel("x [pc]")
  plt.ylabel("y [pc]")  
  ax_divider = make_axes_locatable(plt.gca())
  cax = ax_divider.append_axes("right", size="3%", pad="2%")
  plt.colorbar(im, cax=cax)  

  plt.suptitle("t = {:,.0f} yr".format(t/YR))

  plt.tight_layout()
  # plt.subplots_adjust(top=0.934, bottom=0.032, left=0.024, right=0.992, hspace=0.164, wspace=0.012)

  if save_plots:
    out_fname = os.path.splitext(fname)[0] + ".png"
    plt.savefig(out_fname)
    print(out_fname)
    plt.close()
  else:
    plt.show()
