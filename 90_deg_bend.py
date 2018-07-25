import meep as mp
import numpy as np
import matplotlib.pyplot as plt

cell = mp.Vector3(16, 16, 0)
geometry = [mp.Block(mp.Vector3(12,1,1e20), center = mp.Vector3(-2, 3.5), material = mp.Medium(epsilon = 12)), mp.Block(mp.Vector3(1, 12, 1e20), center = mp.Vector3(3.5, -2), material = mp.Medium(epsilon = 12))]
pml_layers = [mp.PML(1.0)]
resolution = 10

sources = [mp.Source(mp.ContinuousSource(wavelength=2*(11**0.5), width = 20), component=mp.Ez, center = mp.Vector3(-7, 3.5), size = mp.Vector3(0,1))]

sim = mp.Simulation(cell_size = cell, boundary_layers = pml_layers, geometry = geometry, sources = sources, resolution = resolution)

vals = []

def get_slice(sim): 
	center = mp.Vector3(0,3.5)
	size = mp.Vector3(16,0)
	vals.append(sim.get_array(center = center, size = size, component = mp.Ez))

#For GIF output:

#sim.use_output_directory()

#HDF5 output for images: 

#sim.run(mp.at_beginning(mp.output_epsilon), mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)), until = 200)

#PNG output for images: 

#sim.run(mp.at_every(0.6, mp.output_png(mp.Ez, "-Zc dkbluered")), until = 200)

#For y = 3.5 slice output: 

#Entire computational cell: 

#sim.run(mp.at_beginning(mp.output_epsilon), mp.at_every(0.6, get_slice), until = 200)

#Subset of computational cell: 

sim.run(mp.in_volume(mp.Volume(mp.Vector3(0,3.5), size=mp.Vector3(16,0)), mp.to_appended("ez-slice", mp.output_efield_z)), until=200)   

plt.figure(dpi = 100)
plt.imshow(vals, interpolation = 'spline36', cmap = 'RdBu')
plt.axis('off')
plt.show()
