from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as ppl

# parameters
u0 = 0.3 # m/s
ug = 0.04 # m/s
gamma = (86400*2.89)**-1 # s^-1
gammag = (86400*28.9)**-1 # s^-1
f = 1.05e-4 # s^-1
x0 = 0.0
y0 = 0.0

t = np.linspace(0, 10*86400, 1000)


def traj(t, x0, y0, u0, ug, gammag, gamma, f):
    return \
        x0 + ug/gammag * (1 - np.exp(-gammag * t)) + \
        (u0 -ug) * f / (f*f + gamma*gamma) * \
        (gamma/f + np.exp(-gamma*t) * (np.sin(f*t) - gamma/f*np.cos(f*t))), \
        y0 - (u0 - ug) * f / (f*f + gamma*gamma) *  \
        (1 - np.exp(-gamma*t) * (np.cos(f*t) + gamma/f*np.sin(f*t)))

# grid resolution
dx = 250.
dy = dx
trms = np.genfromtxt("results/theoretical/19880101-0000/test_run.asc")
x = trms[:, 2] * dx
y = trms[:, 3] * dy
x -= x[0]
y -= y[0]

xt, yt = traj(t, x0, y0, u0, ug, gammag, gamma, f)
xt -= xt[0]
yt -= yt[0]

fig = ppl.figure()

ax = fig.add_subplot(111)
ax.set_aspect("equal")

ax.plot(xt,yt)
ax.plot(x, y)

fig.show()
