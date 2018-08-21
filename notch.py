import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt
from scipy.optimize import fsolve

#os.system("rm -r notch-out/")

f1 = open('notch.out', 'w')
f2 = open('notch.txt', 'w')

ws = []
Rs = []
Ts = []
Ss = []
NET_LIST = []
Sus = []
Sds = []
NET_LOSS_LIST = []
norm_Sus = []

def notch(w):

	print("#----------------------------------------")
	print("NOTCH WIDTH: %s nanometers" % (w * 1000))
	print("#----------------------------------------")

	# w is the width of the notch in the waveguide

	#Waveguide Math
	a = 4 # length of cell
	h = 0.2 # height of Waveguide
	e = 0.4 # etch fraction [0 -> 1]
	H = 10 * h #height of cell

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

	#---------------------------------------------------------
	'''

	#---------------------------------------------------------
	#FOR GENERATING THE TRANSMITTANCE SPECTRUM

	nfreq = 20  # number of frequencies at which to compute flux

	# reflected flux 1
	refl_fr1 = mp.FluxRegion(center=mp.Vector3(-0.9 * a/2,0), size=mp.Vector3(0,3*h))
	refl1 = sim.add_flux(fcen, df, nfreq, refl_fr1)

	#reflected flux 2
	refl_fr2 = mp.FluxRegion(center=mp.Vector3(-0.6 * a/2,0), size=mp.Vector3(0,3*h))
	refl2 = sim.add_flux(fcen, df, nfreq, refl_fr2)

	#transmitted flux
	tran_fr = mp.FluxRegion(center = mp.Vector3(0.6 * a/2,0), size = mp.Vector3(0,3*h))
	tran = sim.add_flux(fcen, df, nfreq, tran_fr)

	#flux loss above the waveguide
	su_fr = mp.FluxRegion(center = mp.Vector3(0, 3*h), size = mp.Vector3(a,0))
	su = sim.add_flux(fcen, df, nfreq, su_fr)

	#flux loss below the waveguide
	sd_fr = mp.FluxRegion(center = mp.Vector3(0, -3*h), size = mp.Vector3(a,0))
	sd = sim.add_flux(fcen, df, nfreq, sd_fr)

	refl_vals = []
	tran_vals = []

	def get_refl_slice(sim):
	    center = mp.Vector3(-0.9 * a/2,0)
	    size = mp.Vector3(0,3*h)
	    refl_vals.append(sim.get_array(center=center, size=size, component=mp.Ez))

	def get_tran_slice(sim): 
		center = mp.Vector3(0.6 * a/2,0)
		size = mp.Vector3(0,3*h)
		tran_vals.append(sim.get_array(center=center, size=size, component=mp.Ez))

	pt = mp.Vector3(9.75,0)

	sim.run(mp.at_beginning(mp.output_epsilon), 
			mp.at_time(100, get_refl_slice), 
			mp.at_time(100, get_tran_slice), 
			until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

	refl_val = refl_vals[0]
	tran_val = tran_vals[0]

	fund_func = lambda n_eff: sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2))
	first_order_func = lambda n_eff: sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2) - pi / 2)

	initial_guess = n_e

	n_eff0 = fsolve(fund_func, initial_guess)
	n_eff1 = fsolve(first_order_func, initial_guess)

	ky0 = np.absolute(2 * pi / wavelength * sqrt(n_e**2 - n_eff**2))
	ky1 = np.absolute(2 * pi / wavelength * sqrt(n_c**2 - n_eff**2))

	E_fund = lambda y : cos(ky0 * y) if np.absolute(y) < h / 2 else exp(-ky1 * (np.absolute(y) - h / 2))
	E_first_order = lambda y : sin(ky0 * y) if np.absolute(y) < h / 2 else exp(-ky1 * (np.absolute(y) - h / 2))

	# save incident power for reflection planes
	straight_refl1_flux = mp.get_fluxes(refl1)
	straight_refl2_flux = mp.get_fluxes(refl2)

	# save incident power for transmission plane
	straight_tran_flux = mp.get_fluxes(tran)

	# save incident power for loss planes
	straight_su_flux = mp.get_fluxes(su)
	straight_sd_flux = mp.get_fluxes(sd)

	wl = [] #list of wavelengths

	refl1_flux = mp.get_fluxes(refl1)
	refl2_flux = mp.get_fluxes(refl2)
	tran_flux = mp.get_fluxes(tran)
	su_flux = mp.get_fluxes(su)
	sd_flux = mp.get_fluxes(sd)

	flux_freqs = mp.get_flux_freqs(refl1)

	for i in range(nfreq):
	    wl = np.append(wl, 1/flux_freqs[i])

	for ind, elt in enumerate(wl):
	    #print(round(elt, 4))
	    if round(elt, 3) == 0.637:
	        #print("ALERT: MATCH FOUND")
	        index = ind

	R = -refl1_flux[index] / (refl2_flux[index] - refl1_flux[index])
	T = tran_flux[index] / (refl2_flux[index] - refl1_flux[index])
	S = (refl2_flux[index] - tran_flux[index]) / (refl2_flux[index] - refl1_flux[index])
	Su = su_flux[index] / (refl2_flux[index] - refl1_flux[index])
	Sd = -sd_flux[index] / (refl2_flux[index] - refl1_flux[index])

	norm_Su = S * Su / (Su + Sd)

	NET = round((R + T + S) * 100, 0)
	if NET > 100.0:
	    NET = 100.0

	NET_LOSS = round((Su + Sd) / S * 100, 0)
	if NET_LOSS > 100.0:
		NET_LOSS = 100.0

	'''
	np.append(ws, [w * 1000])
	np.append(Rs, [R * 100])
	np.append(Ts, [T * 100])
	np.append(Ss, [S * 100])
	np.append(NET_LIST, [NET])
	np.append(Sus, [Su * 100])
	np.append(Sds, [Sd * 100])
	np.append(NET_LOSS_LIST, [NET_LOSS])
	np.append(norm_Sus, [norm_Su * 100])
	'''

	ws.append(w*1000)
	Rs.append(R*100)
	Ts.append(T*100)
	Ss.append(S*100)
	NET_LIST.append(NET)
	Sus.append(Su*100)
	Sds.append(Sd*100)
	NET_LOSS_LIST.append(NET_LOSS)
	norm_Sus.append(norm_Su*100)

	f1.write("--------------------------------------------------- \n")
	f1.write("Notch Width: %s nanometers \n" % (w * 1000))
	f1.write("Reflection Percentage: %s \n" % (R * 100))
	f1.write("Transmission Percentage: %s \n" % (T * 100))
	f1.write("Total Loss Percentage: %s \n" % (S * 100))
	f1.write("Percentage of Light Accounted For: %s \n" % (NET))
	f1.write("Upper Loss Percentage: %s \n" % (Su * 100))
	f1.write("Lower Loss Percentage: %s \n" % (Sd * 100))
	f1.write("Percentage of Total Loss Accounted For: %s \n" % (NET_LOSS))
	f1.write("Normalized Upper Loss Percentage: %s \n" % (norm_Su * 100))
	f1.write("--------------------------------------------------- \n")

	sim.reset_meep()

	#-------------------------------------------------------------

for notch_index in range(4, 32, 2):
	notch_width = notch_index / 100
	notch(notch_width)

f2.write("%s \n" % (ws))
f2.write("%s \n" % (Rs))
f2.write("%s \n" % (Ts))
f2.write("%s \n" % (Ss))
f2.write("%s \n" % (Sus))
f2.write("%s \n" % (Sds))
f2.write("%s \n" % (norm_Sus))

f1.close()
f2.close()

plt.plot(ws, Rs,'bo-',label='reflectance')
plt.plot(ws, Ts,'ro-',label='transmittance')
plt.plot(ws, Ss,'go-',label='net loss')
plt.plot(ws, Sus, 'co-', label='upper loss')
plt.plot(ws, Sds, 'mo-', label='lower loss')
plt.plot(ws, norm_Sus, 'ko-', label='normalized upper loss')
plt.axis([40.0, 300.0, 0.0, 100.0])
plt.xlabel("wavelength (nm)")
plt.ylabel("percentage (%)")
plt.legend(loc="center right")
plt.show()
