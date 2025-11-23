import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from time import time

# Rössler-System
def Rossler(x, y, z, a, b, c):
    x_dot = - y - z
    y_dot = x + a * y
    z_dot = b + x * z - c * z
    return x_dot, y_dot, z_dot

# Parameter des Systems
a, b = 0.2, 0.2
c_init = 2.9
c_min = 2
c_max = 7

# Startpunkt der Trajektorie
starting_point = (1.0, 1.0, 1.0)

# Schrittweite und -zahl der Simulation
dt = 0.01
step_count = 30000

# Koordinatenschranke
plot_limit = 6.0

# Erstellen der Arrays
xs=np.empty((step_count + 1,))
ys=np.empty((step_count + 1,))
zs=np.empty((step_count + 1,))

# Erstellen von Figure- und Axes-Objekten
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
#print(ax)
fig.subplots_adjust(bottom=0.2, left=0.15)
l, = ax.plot([], [], [], lw=0.5)

plt.title("Rössler-Attraktor")

from mpl_toolkits.mplot3d.art3d import Line3DCollection


def colored_3d_line(ax, points, c):
    # Bau der Segmente
    points_rs = points.T.reshape(-1, 1, 3)
    #print(points)
    c = c[:len(points)]
    segs = np.concatenate([points_rs[:-1], points_rs[1:]], axis=1)
    #print(segs)
    # Normalize mappt die segmente linear auf [0,1]; für colormap
    lc = Line3DCollection(segs, cmap="plasma", norm=plt.Normalize(0,1))
    lc.set_array(c)
    lc.set_linewidth(0.8)
    ax.add_collection3d(lc)

color = np.linspace(0, 2, step_count + 1)

# Update-Funktion für c-Parameter
def c_update(c):

    t1 = time()

    xs[0], ys[0], zs[0] = starting_point
    # derivative = []

    for i in range(step_count):
        der = Rossler(xs[i], ys[i], zs[i], 0.2, 0.2, c)
        # derivative.append(der)

        xs[i+1] = xs[i] + (der[0]*dt)
        ys[i+1] = ys[i] + (der[1]*dt)
        zs[i+1] = zs[i] + (der[2]*dt)

    t15 = time()
    
    
    # derivative_arr = np.array(derivative)
    # print(derivative_arr)
    # # derivative_arr = 25*np.ones((step_count,3))
    # derivative_norm = np.absolute(derivative_arr.T[0]) + np.absolute(derivative_arr.T[0]) + np.absolute(derivative_arr.T[0])

    t2 = time()

    # points_eff = []
    # print(len(derivative_norm))
    # sum = 0
    # acc = 1
    # for i in range(step_count):
    #     sum += derivative_norm[i]
    #     if sum > acc:
    #         points_eff.append([xs[i], ys[i], zs[i]])
    #         sum -= acc * (sum // acc)
    # print(len(points_eff))

    # # a = [points_eff[j][2] == zs[j] for j in range(len(points_eff))]
    # # print(a)

    # points_eff_arr = np.array(points_eff).T
    # # b = (points_eff_arr == np.array([xs, ys, zs])[:len(points_eff)])
    # # print(b)
    # # print(points_eff_arr)
    # # print(len(points_eff_arr))


    for coll in ax.collections:
        coll.remove()

    t3 = time()
    print(np.array([xs,ys,zs])[:100])

    colored_3d_line(ax, np.array([xs,ys,zs])[:100], color)
    fig.canvas.draw_idle()

    t4 = time()

    print(t15-t1, t2-t15, t3-t2, t4-t3)

# Initialisierung von ax mit Trajektorie des Rösler Systems für c = c_init
c_update(c_init)

# Update-Funktion für Zoom
def zoom_update(exponent):
    bound = plot_limit * (10 ** -exponent)
    ax.set(xlim3d=(-bound, bound))
    ax.set(ylim3d=(-bound, bound))
    ax.set(zlim3d=(0, 2 * bound))
    fig.canvas.draw_idle()

# Setup der anfänglichen Koordinatenschranken (-plot_limit, plot_limit)^2 x (0, 2*plot_limit)
zoom_update(0)

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