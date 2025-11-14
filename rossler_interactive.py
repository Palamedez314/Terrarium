import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Rössler-System
def Rossler(x, y, z, a, b, c):
    x_dot = - y - z
    y_dot = x + a * y
    z_dot = b + x * z - c * z
    return x_dot, y_dot, z_dot

# Initialisierung der Parameter
a, b = 0.2, 0.2
c_init = 2.9
c_min = 2
c_max = 7
starting_point = (1.0, 1.0, 1.0)
dt = 0.01
step_count = 50000

# Erstellen der Arrays
xs=np.empty((step_count + 1,))
ys=np.empty((step_count + 1,))
zs=np.empty((step_count + 1,))

# Erstellen von Figure- und Axes-Objekten
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
fig.subplots_adjust(bottom=0.2, left=0.15)
l, = ax.plot([], [], [], lw=0.5)

# Setup der Axes-Eigenschaften
ax.set(xlim3d=(-6, 6), xlabel='x')
ax.set(ylim3d=(-6, 6), ylabel='y')
ax.set(zlim3d=(0, 6), zlabel='z')

# Update-Funktion für c-Parameter
def c_update(c):
    xs[0], ys[0], zs[0] = starting_point
    for i in range(step_count):
        x_dot, y_dot, z_dot = Rossler(xs[i], ys[i], zs[i], 0.2, 0.2, c)
        xs[i+1] = xs[i] + (x_dot*dt)
        ys[i+1] = ys[i] + (y_dot*dt)
        zs[i+1] = zs[i] + (z_dot*dt)
    l.set_data_3d(xs, ys, zs)
    fig.canvas.draw_idle()

# Initialisierung von ax mit Trajektorie des Rösler Systems für c = c_init
c_update(c_init)

# Update-Funktion für Zoom
def zoom_update(exp):
    bound = 6 * (10 ** -exp)
    ax.set(xlim3d=(-bound, bound))
    ax.set(ylim3d=(-bound, bound))
    ax.set(zlim3d=(0, 2 * bound))
    fig.canvas.draw_idle()

# Erstellen der Slider
hslideraxis = fig.add_axes([0.35, 0.1, 0.55, 0.03])
hslider = Slider(hslideraxis, label='c-Parameter',
                valmin=c_min, valmax=c_max, valinit=c_init)
hslider.on_changed(c_update)

vslideraxis = fig.add_axes([0.08, 0.25, 0.03, 0.65])
vslider = Slider(vslideraxis, label='Zoom',
                valmin=-1, valmax=1, valinit=0, orientation="vertical")
vslider.on_changed(zoom_update)

# Anzeigen des Plots
plt.show()