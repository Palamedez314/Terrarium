import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox, Button, CheckButtons
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# fix für komische 3D-Rotation für neue mathlib-Versionen:
# import matplotlib as mpl
# mpl.rcParams['axes3d.mouserotationstyle'] = 'azel' 

# Lorenz-System
def Lorenz(x, y, z, sigma, b, r):
    x_dot = sigma * (y - x)
    y_dot = r * x - y - x * z
    z_dot = x*y - b * z
    return x_dot, y_dot, z_dot

# Parameter des Systems
sigma, b = 10, 8/3
r_init = 28
r_min = 0
r_max = 35

# Startpunkt der Trajektorie zu Beginn
starting_point = [1.0,1.0,1.0]

# Schrittweite und -zahl der Simulation
dt = 0.01
step_count = 50000

# Anzahl verschiedene Farben im Farbverlauf
line_len = 1000

# Nur den Grenzzyklus zu Beginn anzeigen
only_limit_cycle = False

# Koordinatenschranke
plot_limit = 25.0

# Erstellen der Arrays
xs=np.empty((step_count + 1,))
ys=np.empty((step_count + 1,))
zs=np.empty((step_count + 1,))

# Erstellen von Figure- und Axes-Objekten
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
fig.subplots_adjust(left=0.1, bottom=0.2)
l, = ax.plot([], [], [], lw=0.5)
plt.title("Lorenz-Attraktor")

# Erstellen und Plotten verschiedenfarbiger Liniensegmente der Länge line_len
def colored_3d_line(x, y, z, ax):

    # Zusammenbauen der Segmente
    points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    segs = np.concatenate(
        [points[i::(line_len-1)]
            [:(len(points)//(line_len))]
         for i in range(line_len) ], axis=1)

    # Farbverlauf Blau -> Rot
    rgba = np.array(
        [(i / len(segs),
          0,
          1 - i / len(segs),
          1)
          for i in range(len(segs))])

    # Abschneidepunkt des Arrays für den Grenzzyklus
    devider = int(0.7*(step_count/line_len))

    # Erstellen und übergeben der 3D-LineCollections
    lc = Line3DCollection(segs, colors=rgba)

    lc_only_limitCycle = Line3DCollection(
        segs[devider:], colors=rgba[devider:])
    
    lc.set_linewidth(0.8)
    lc_only_limitCycle.set_linewidth(0.8)
    ax.add_collection3d(lc)
    ax.add_collection3d(lc_only_limitCycle)

    # Sichtbarkeit der zwei Line3DCollections 
    lc.set_visible(not only_limit_cycle)
    lc_only_limitCycle.set_visible(only_limit_cycle)

# Update-Funktion für r-Parameter
def r_update(r):
    xs[0], ys[0], zs[0] = starting_point

    # Iterierte Simulation des Rössler-Systems
    for i in range(step_count):
        x_dot, y_dot, z_dot = Lorenz(
            xs[i], ys[i], zs[i], sigma, b, r)
        
        xs[i+1] = xs[i] + (x_dot*dt)
        ys[i+1] = ys[i] + (y_dot*dt)
        zs[i+1] = zs[i] + (z_dot*dt)

    # Plotten der neuen Trajektorie
    for coll in ax.collections:
        coll.remove() 
    colored_3d_line(xs, ys, zs, ax)
    fig.canvas.draw_idle()

# Initialisierung von ax mit Trajektorie des Rösler Systems für r = r_init
r_update(r_init)

# Update-Funktion für Zoom
def zoom_update(exponent):
    bound = plot_limit * (10 ** -exponent)
    ax.set(xlim3d=(-bound, bound))
    ax.set(ylim3d=(-bound, bound))
    ax.set(zlim3d=(0, 2 * bound))
    fig.canvas.draw_idle()

# Setup der anfänglichen Koordinatenschranken (-plot_limit, plot_limit)^2 x (0, 2*plot_limit)
zoom_update(0)

# Update-Funktion für Checkbox
def callback(_label):
    global only_limit_cycle
    coll1, coll2 = ax.collections[0], ax.collections[1]
    only_limit_cycle = not only_limit_cycle
    coll1.set_visible(not only_limit_cycle)
    coll2.set_visible(only_limit_cycle)
    fig.canvas.draw_idle()

# Erstellen der Slider
hslideraxis = fig.add_axes([0.35, 0.1, 0.55, 0.03])
hslider = Slider(hslideraxis, label='r-Parameter',
                valmin=r_min, valmax=r_max,
                valinit=r_init)
hslider.on_changed(r_update)

vslideraxis = fig.add_axes([0.08, 0.25, 0.03, 0.65])
vslider = Slider(vslideraxis, label='Zoom',
                valmin=-0.5, valmax=2,
                valinit=0, orientation="vertical")
vslider.on_changed(zoom_update)

# Erstellen von TextBoxen und Button
error_text = plt.text(1.97,-0.5,"",horizontalalignment="center",color="tab:red")
def change_starting_point(text,index):
    try:
        starting_point[index] = float(text)
    except ValueError:
        error_text.set(text="Bitte Zahlen eingeben")
    else:
        error_text.set(text="")

textbox_x_ax = fig.add_axes([0.85,0.6,0.1,0.05])
textbox_x = TextBox(textbox_x_ax,"X",
                    str(starting_point[0]))
textbox_x.on_submit(lambda x_str: 
                    change_starting_point(x_str,0))

textbox_y_ax = fig.add_axes([0.85,0.55,0.1,0.05])
textbox_y = TextBox(textbox_y_ax,"Y",
                    str(starting_point[1]))
textbox_y.on_submit(lambda y_str: 
                    change_starting_point(y_str,1))

textbox_z_ax = fig.add_axes([0.85,0.5,0.1,0.05])
textbox_z = TextBox(textbox_z_ax,"Z",
                    str(starting_point[2]))
textbox_z.on_submit(lambda z_str: 
                    change_starting_point(z_str,2))

submit_button_ax = fig.add_axes([0.86,0.45,0.08,0.04])
submit_button = Button(submit_button_ax,"Submit")
submit_button.on_clicked(lambda _clicked : 
                         r_update(hslider.val))

# Erstellen der Checkbox
rax = ax.inset_axes([0.95,0.08,0.5,0.12])
check = CheckButtons(
    ax=rax,
    labels=["nur Grenzzyklus\n anzeigen"],
    actives=[only_limit_cycle],
    label_props={'color': ["red"]},
    frame_props={'edgecolor': ["red"]},
    check_props={'facecolor': ["red"]})
check.on_clicked(callback)


# Anzeigen des Plots
plt.show()