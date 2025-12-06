# Program 08b: The Lorenz attractor. See Figure 8.11.
# In this case, the odeint numerical solver was used to solve the ODE.
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Lorenz parameters and initial conditions
sigma, beta, rho = 10, 2.667, 28
x0, y0, z0 = 0, 1, 1.05
# Maximum time point and total number of time points
tmax, n = 100, 10000

def Lorenz(X, t, sigma, beta, rho):
    """The Lorenz equations."""
    x, y, z = X
    dx = -sigma * (x - y)
    dy = rho * x - y - x * z
    dz = -beta * z + x * y
    return dx, dy, dz

t = np.linspace(0, tmax, n)
f = odeint(Lorenz, (x0, y0, z0), t, args=(sigma, beta, rho))
x, y, z = f.T

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, 'b-', lw=0.5)
ax.set_xlabel('x', fontsize=15)
ax.set_ylabel('y', fontsize=15)
ax.set_zlabel('z', fontsize=15)
plt.tick_params(labelsize=15)
ax.set_title('Lorenz Attractor', fontsize=15)
plt.show()