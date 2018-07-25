import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt

os.system("rm -r notch-out/")

#Waveguide Math
a = 4 # length of waveguide
h = 0.2 # height of Waveguide
e = 0.4 # etch fraction [0 -> 1]
w = 0.1 # width of the notch
H = 10 * h

h2 = h * (1 - e) #height of the etched region
d = 1 #spacing of the waveguides
posd = d / 2 #pos half d
negd = -d / 2 #neg half d

ang = 80 #angle of the sidewall
angrad = ang * pi / 180
bottomoffset = h / tan(angrad)
c = bottomoffset / 2
r = sqrt(bottomoffset ** 2 + h ** 2)

wavelength = 0.6372
fcen = 1 / wavelength
df = 0.05

n_e = 2.2012 # refractive index inside waveguide (at wavelength 0.6372)
n_c = 1.4569 # refractive index outside waveguide (at wavelength 0.6372)

#Waveguide Geometry
cell = mp.Vector3(a, H)

default_material = mp.Medium(epsilon = n_c ** 2)

geometry = [mp.Block(mp.Vector3(a + 2 * h, h),
            center = mp.Vector3(0,0),
            material = mp.Medium(epsilon = n_e ** 2)),
            mp.Block(mp.Vector3(w, e * h),
            center = mp.Vector3(0, h2 / 2),
            material = mp.Medium())]

sources = [mp.EigenModeSource(mp.GaussianSource(frequency = fcen, fwidth = df),
                              size = mp.Vector3(0,H),
                              center = mp.Vector3(-0.75 * a / 2, 0))]

pml_layers = [mp.PML(0.2)]

resolution = 50

sim = mp.Simulation(cell_size = cell,
                    boundary_layers = pml_layers,
                    geometry = geometry,
                    sources = sources,
                    resolution = resolution)
'''
#--------------------------------------------------
#FOR DISPLATYING THE GEOMETRY

sim.run(until = 200)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
plt.figure(dpi=100)
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
#plt.axis('off')
plt.show()

quit()
#----------------------------------------------------
'''
#------------------------------------------------------
#FOR GENERATING THE ELECTRIC FIELD GIF
#Note: After running this program, write the following commands in Terminal:
    # $ source deactivate mp
    # $ cd notch-out/
    # $ python ../NotchIP.py

sim.use_output_directory()
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
        until = 200)
#sim.run(mp.at_every(0.6 , mp.output_png(mp.Ez, "-Zc dkbluered")), until=200)

quit()
#---------------------------------------------------------
