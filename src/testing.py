#Toolbox Testing

import astro_toolbox as toolbox
import numpy as np
import math as m
import matplotlib.pyplot as plt

mu = 1
ae = 1

r = np.array([0, 1, 1])
v = np.array([-1, 0, -0.05])

params = toolbox.rv_orbparams(r[0],r[1],r[2],v[0],v[1],v[2],1)
pqw = toolbox.param_to_pqw(params[1],params[2],100)
IJK = toolbox.pqw_to_ijk(pqw, params[3], params[4], params[5])


#Break IJK coords into individual coord lists for plotting
icoord = []
jcoord = []
kcoord = []
for point in IJK:
    icoord.append(point[0])
    jcoord.append(point[1])
    kcoord.append(point[2])

#generate date for wireframe Earth:
u,v =  np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    
px = ae*np.cos(u)*np.sin(v)
py = ae*np.sin(u)*np.sin(v)
pz = ae*np.cos(v)

#Frame sizer
ijk_max_list = (max(icoord), max(jcoord), max(kcoord), abs(min(icoord)), abs(min(jcoord)), abs(min(kcoord)))
max_ijk = max(ijk_max_list)

#Plot orbit in ijk frame
fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(111, projection = '3d')
ax2.set_xlim(-max_ijk, max_ijk)
ax2.set_ylim(-max_ijk, max_ijk)
ax2.set_zlim(-max_ijk, max_ijk)
ax2.plot_wireframe(px, py, pz, rstride=1, cstride=1, color='green', linewidth=0.5)
ax2.scatter(icoord, jcoord, kcoord)
ax2.title.set_text('IJK Frame Plot of Orbit \n Type: %s' %params[7])

plt.show()