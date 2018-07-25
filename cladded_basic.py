import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

os.system("rm -r cladded_basic-out/")

cell = mp.Vector3(20, 10)

geometry = [mp.Block(mp.Vector3(1e20, 1),
            center = mp.Vector3(0,0),
            material = mp.Medium(epsilon = 12)),
            mp.Block(mp.Vector3(4, 4),
            center = mp.Vector3(0,0),
            e1 = mp.Vector3(3,2),
            e2 = mp.Vector3(3,-2),
            material = mp.Medium(epsilon = 24))]

fcen = 0.15  # pulse center frequency
df = 0.1     # pulse width (in frequency)

sources = [mp.EigenModeSource(mp.GaussianSource(frequency = fcen, fwidth = df),
                              size = mp.Vector3(0,1),
                              center = mp.Vector3(-7,0))]

pml_layers = [mp.PML(1.0)]

resolution = 10

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
'''
#------------------------------------------------------
#FOR GENERATING THE ELECTRIC FIELD GIF
#Note: After running this program, write the following commands in Terminal:
    # $ source deactivate mp
    # $ cd cladded_basic-out/
    # $ python ../ClaBasIP.py

sim.use_output_directory()
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
        until = 200)
#sim.run(mp.at_every(0.6 , mp.output_png(mp.Ez, "-Zc dkbluered")), until=200)

quit()
#---------------------------------------------------------
'''

#---------------------------------------------------------
#FOR GENERATING THE TRANSMITTANCE SPECTRUM

nfreq = 100  # number of frequencies at which to compute flux

# reflected flux 1
refl_fr1 = mp.FluxRegion(center=mp.Vector3(-9,0), size=mp.Vector3(0,2))
refl1 = sim.add_flux(fcen, df, nfreq, refl_fr1)

#reflected flux 2
refl_fr2 = mp.FluxRegion(center=mp.Vector3(-5,0), size=mp.Vector3(0,2))
refl2 = sim.add_flux(fcen, df, nfreq, refl_fr2)

#transmitted flux
tran_fr = mp.FluxRegion(center = mp.Vector3(5,0), size = mp.Vector3(0,2))
tran = sim.add_flux(fcen, df, nfreq, tran_fr)

pt = mp.Vector3(9.75,0)

sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

# save incident power for reflection planes
straight_refl1_flux = mp.get_fluxes(refl1)
straight_refl2_flux = mp.get_fluxes(refl2)

# save incident power for transmission plane
straight_tran_flux = mp.get_fluxes(tran)

wl = [] #list of wavelengths
Rs = [] #reflectance spectrum
Ts = [] #transmittance spectrum
Is = [] #spectrum of incident light

refl1_flux = mp.get_fluxes(refl1)
refl2_flux = mp.get_fluxes(refl2)
tran_flux = mp.get_fluxes(tran)

flux_freqs = mp.get_flux_freqs(refl1)

for i in range(nfreq):
    wl = np.append(wl, 1/flux_freqs[i])
    Rs = np.append(Rs, -refl1_flux[i])
    Ts = np.append(Ts, tran_flux[i])
    Is = np.append(Is, refl2_flux[i] - refl1_flux[i])

plt.plot(wl,Rs,'bo-',label='reflectance')
plt.plot(wl,Ts,'ro-',label='transmittance')
plt.plot(wl,Is-Rs-Ts,'go-',label='loss')
plt.axis([5.0, 10.0, 0.0, 6.5])
plt.xlabel("wavelength (μm)")
plt.legend(loc="upper right")
plt.show()

quit()
#-------------------------------------------------------------
