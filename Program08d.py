# Programs 08d: Animation of a Chua circuit bifurcation.
# You can watch a YouTube animation on the web.
# Search for Chua circuit AND oscilloscope.
from matplotlib import pyplot as plt
from matplotlib.animation import ArtistAnimation
import numpy as np
from scipy.integrate import odeint

fig=plt.figure()
mo = -1/7
m1 = 2/7
tmax = 100

def chua(x, t):
    return [a * (x[1] - (m1 * x[0] + (mo - m1) / 2 *
    (np.abs(x[0] + 1) - \
    np.abs(x[0] - 1)))), x[0] - x[1] + x[2], -15 * x[1]]

time = np.arange(0, tmax, 0.01)
x0=[1.96, -0.0519, -3.077]
myimages = []

for a in np.arange(8, 11, 0.1):
    xs = odeint(chua, x0, time)
    imgplot = plt.plot(xs[:, 0], xs[:, 1], "r-")
    myimages.append(imgplot)
my_anim=ArtistAnimation(fig, myimages, interval=500, blit=False, repeat_delay=500)
plt.show()
