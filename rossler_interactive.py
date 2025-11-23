import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, Button

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
c_current = c_init

# Startpunkt der Trajektorie
starting_point = [1.0, 1.0, 1.0]

# Schrittweite und -zahl der Simulation
dt = 0.05
step_count = 5000


# Koordinatenschranke
plot_limit = 6.0

# Erstellen der Arrays
xs=np.empty((step_count + 1,))
ys=np.empty((step_count + 1,))
zs=np.empty((step_count + 1,))

# Erstellen von Figure- und Axes-Objekten
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
print(ax)
fig.subplots_adjust(bottom=0.2, left=0.15)
l, = ax.plot([], [], [], lw=0.5)

plt.title("Rössler-Attraktor")

import warnings #wieder entfernen
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def colored_line(x, y, z, c, ax, **lc_kwargs):
    
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Kp was das macht, war schon da
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)


    #Die "50" sind nur zum testen
    c=c[:5]

    x = np.asarray(x)[:5]
    y = np.asarray(y)[:5]
    z = np.asarray(z)[:5]
    #print(x)
    #print(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))
    z_midpts = np.hstack((z[0], 0.5 * (z[1:] + z[:-1]), z[-1]))

    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1], z_midpts[:-1]))[:, np.newaxis, :]
    #print("coord_start:\n", coord_start)
    coord_mid = np.column_stack((x, y, z))[:, np.newaxis, :]
    #print("coord_mid:\n", coord_mid)
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:], z_midpts[1:]))[:, np.newaxis, :]
    #print("coord_end:\n", coord_end)
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)
    #print("coord_segments:\n", segments)

    lc = Line3DCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment
    #print(lc)
    
    return ax.add_collection3d(lc)

def colored_3d_line(x, y, z, c, ax, **lc_kwargs):
    # Bau der Segmente
    points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    #print(points)
    segs = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segs)
    # Normalize mappt die segmente linear auf [0,1]; für colormap
    lc = Line3DCollection(segs, norm=plt.Normalize(0,1),**lc_kwargs)
    lc.set_array(c)
    lc.set_linewidth(0.8)
    ax.add_collection3d(lc)
    return lc

color = np.linspace(0, 2, 5001)


# Update-Funktion für c-Parameter
def c_update(c):
    print(starting_point)
    c_current = c
    xs[0], ys[0], zs[0] = starting_point
    for i in range(step_count):
        x_dot, y_dot, z_dot = Rossler(xs[i], ys[i], zs[i], 0.2, 0.2, c)
        xs[i+1] = xs[i] + (x_dot*dt)
        ys[i+1] = ys[i] + (y_dot*dt)
        zs[i+1] = zs[i] + (z_dot*dt)


    for coll in ax.collections:
        coll.remove() 
    _lines = colored_3d_line(xs, ys, zs, color, ax, cmap="hsv")
    fig.canvas.draw_idle()

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

# Erstellen der TextBoxen
textbox_x_ax = fig.add_axes([0.85,0.9,0.1,0.05])
textbox_x = TextBox(textbox_x_ax,"X-Wert","1.0")
textbox_x.on_submit(lambda x_str : (starting_point.pop(0),starting_point.insert(0,float(x_str))))

textbox_y_ax = fig.add_axes([0.85,0.85,0.1,0.05])
textbox_y = TextBox(textbox_y_ax,"Y-Wert","1.0")
textbox_y.on_submit(lambda y_str : (starting_point.pop(1),starting_point.insert(1,float(y_str))))

textbox_z_ax = fig.add_axes([0.85,0.8,0.1,0.05])
textbox_z = TextBox(textbox_z_ax,"Z-Wert","1.0")
textbox_z.on_submit(lambda z_str : (starting_point.pop(2),starting_point.insert(2,float(z_str))))

submit_button_ax = fig.add_axes([0.86,0.75,0.08,0.04])
submit_button = Button(submit_button_ax,"Submit")
submit_button.on_clicked(lambda _clicked : c_update(c_current))


# Anzeigen des Plots
plt.show()