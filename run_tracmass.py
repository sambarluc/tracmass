from tmpyui import tmsim
import numpy as np
d=tmsim("/home/cimatori/Work/tracmass",
        "MITgcm_leman_mp", "template_input", "MITgcm_leman.in",
        "/home/cimatori/Work/mitgcm/lac_leman/run_2016")

# we seed along a transect off the Rhone outflow
npts = 5
p0 = (553319., 137896.) # km CH03
p1 = (555803., 139517.)
z0 = 2.5
z1 = 100.
dz = 5.0
zs = np.arange(z0, z1, dz)
xs = np.linspace(p0[0], p1[0], npts)
ys = np.linspace(p0[1], p1[1], npts)
X, Z = np.meshgrid(xs, zs)
Y, Z = np.meshgrid(ys, zs)

d.ijklims(ii=X.ravel(), jj=Y.ravel(), kk=Z.ravel(), inds=False,)
d.run(3, ni0=15)

from numpy import savetxt

savetxt("exit_codes_rtrm.txt", d.exit)
