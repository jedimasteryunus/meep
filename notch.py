import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, tan, sqrt
from scipy.optimize import fsolve

os.system("rm -r notch-out/")

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

n_eff_funds = []
n_eff_firsts = []

r00s = []
r01s = []
r10s = []
r11s = []

t00s = []
t01s = []
t10s = []
t11s = []

su0s = []
sd0s = []
su1s = []
sd1s = []

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

	pml_layers = [mp.PML(0.2)]

	resolution = 50

	r00 = None
	r01 = None
	r10 = None
	r11 = None

	t00 = None
	t01 = None
	t10 = None
	t11 = None

	su0 = None
	sd0 = None
	su1 = None
	sd1 = None

	for mode in [0, 1]:

		if mode == 0:
			eig_parity = mp.EVEN_Y	# Fundamental
			print("-----------")
			print("MODE TYPE: FUNDAMENTAL")

		else:
			eig_parity = mp.ODD_Y #First Order
			print("-----------")
			print("MODE TYPE: FIRST ORDER")

		sources = [	mp.EigenModeSource(
						mp.GaussianSource(	frequency = fcen,
											fwidth = df),
											size = mp.Vector3(0,H),
											center = mp.Vector3(-0.75 * a / 2, 0),
					eig_parity = eig_parity) ]

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

		# ------------------------ CODE FOR SEPARATING FUND AND FIRST ORDER MODE STARTS HERE ------------------------

		refl_vals = []
		tran_vals = []

		def get_refl_slice(sim):
		    center = mp.Vector3(-0.9 * a/2,0)
		    size = mp.Vector3(0,H)
		    refl_vals.append(sim.get_array(center=center, size=size, component=mp.Ez, cmplx = True))

		def get_tran_slice(sim):
			center = mp.Vector3(0.6 * a/2,0)
			size = mp.Vector3(0,H)
			tran_vals.append(sim.get_array(center=center, size=size, component=mp.Ez, cmplx = True))

		pt = mp.Vector3(9.75,0)

		gif = False

		if gif and w == 0.1:
			sim.use_output_directory()
			sim.run(mp.at_beginning(mp.output_epsilon),
					mp.at_end(get_refl_slice),
					mp.at_end(get_tran_slice), until=100)
			sim.run(mp.at_every(wavelength / 20, mp.output_efield_z), until=wavelength)
			sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))
		else:
			sim.run(mp.at_beginning(mp.output_epsilon),
					mp.at_time(100, get_refl_slice),
					mp.at_time(100, get_tran_slice),
					until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

		refl_val = refl_vals[0]
		tran_val = tran_vals[0]

		#In the fsolve we need to somehow ensure that n_eff satisfies n_e <= n_eff <= n_c

		def fund_func(n_eff):
			if n_eff > n_c and n_eff < n_e:
				return sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2))

		def first_order_func(n_eff):
			if n_eff > n_c and n_eff < n_e:
				return sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2) - pi / 2)

		initial_guess = (n_c + n_e) / 2

		n_eff_fund = 	fsolve(fund_func, initial_guess)
		n_eff_first = 	fsolve(first_order_func, initial_guess)

		if len(n_eff_funds) == 0:
			n_eff_funds.append(n_eff_fund[0])

		if len(n_eff_firsts) == 0:
			n_eff_firsts.append(n_eff_first[0])

		ky0_fund =  	np.absolute(2 * pi / wavelength * sqrt(n_e**2 - n_eff_fund **2))
		ky0_first = 	np.absolute(2 * pi / wavelength * sqrt(n_e**2 - n_eff_first**2))

		ky1_fund =  	np.absolute(2 * pi / wavelength * sqrt(n_eff_fund **2 - n_c**2))
		ky1_first = 	np.absolute(2 * pi / wavelength * sqrt(n_eff_first**2 - n_c**2))

		E_fund = 		lambda y : cos(ky0_fund  * y) if np.absolute(y) < h / 2 else cos(ky0_fund  * h/2) * np.exp(-ky1_fund  * (np.absolute(y) - h / 2))
		E_first_order = lambda y : sin(ky0_first * y) if np.absolute(y) < h / 2 else sin(ky0_first * h/2) * np.exp(-ky1_first * (np.absolute(y) - h / 2))

		y_list = np.arange(-H/2, H/2, 1/50)

		#print("Y LIST: ", y_list)
		#print("SIZE OF Y LIST: ", y_list.size)

		E_fund_vec = np.zeros(y_list.size)
		E_first_order_vec = np.zeros(y_list.size)

		for index in range(y_list.size):
			y = y_list[index]
			E_fund_vec[index] = E_fund(y)
			E_first_order_vec[index] = E_first_order(y)

		#print("E VECTOR: ", refl_val)
		#print("E0 VECTOR: ", E_fund_vec)

		fund_refl_amp = 		np.dot(refl_val, E_fund_vec) 		/ np.dot(E_fund_vec, E_fund_vec)
		first_order_refl_amp = 	np.dot(refl_val, E_first_order_vec) / np.dot(E_first_order_vec, E_first_order_vec)
		fund_tran_amp = 		np.dot(tran_val, E_fund_vec) 		/ np.dot(E_fund_vec, E_fund_vec)
		first_order_tran_amp = 	np.dot(tran_val, E_first_order_vec) / np.dot(E_first_order_vec, E_first_order_vec)

		fund_refl_power = 			np.abs(fund_tran_amp)			** 2
		first_order_refl_power = 	np.abs(first_order_refl_amp) 	** 2
		fund_tran_power = 			np.abs(fund_tran_amp) 			** 2
		first_order_tran_power = 	np.abs(first_order_tran_amp) 	** 2

		fund_refl_ratio = 			fund_refl_power 		/ (fund_refl_power + first_order_refl_power)
		first_order_refl_ratio = 	first_order_refl_power 	/ (fund_refl_power + first_order_refl_power)
		fund_tran_ratio = 			fund_tran_power 		/ (fund_tran_power + first_order_tran_power)
		first_order_tran_ratio = 	first_order_tran_power 	/ (fund_tran_power + first_order_tran_power)

		print("Percentage of reflected light in fundamental mode: ", 	fund_refl_ratio * 100)
		print("Percentage of reflected light in first order mode: ", 	first_order_refl_ratio * 100)
		print("Percentage of transmitted light in fundamental mode: ", 	fund_tran_ratio * 100)
		print("Percentage of transmitted light in first order mode: ", 	first_order_tran_ratio * 100)

		# ------------------------ CODE FOR SEPARATING FUND AND FIRST ORDER MODE ENDS HERE ------------------------

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

		R = 	-refl1_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		T = 	 tran_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		S = 	(refl2_flux[index] - tran_flux[index]) 	/ (refl2_flux[index] - refl1_flux[index])
		Su = 	 su_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		Sd = 	-sd_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])

		LT = -w/2 + 0.9 * a/2	# Distance from the end of the notch to the transmission monitor
		LR = 0.6 * a/2 - w/2	# Distance from the beginning of the notch to the reflection monitor

		r = sqrt(R)

		r_fund = 	(r  * fund_refl_ratio) * (fund_refl_amp / np.abs(fund_refl_amp)) * (np.exp(2j*pi * LR * n_eff_fund  / wavelength))
		# The amplitude ... times the phase ... accounting for the distance to the detector.
		#print("Re(r_fund): ", np.real(r_fund))
		#print("Im(r_fund): ", np.imag(r_fund))

		r_first = 	(r * first_order_refl_ratio) * (first_order_refl_amp / np.abs(first_order_refl_amp)) * (np.exp(2j*pi * LR * n_eff_first / wavelength))
		# The amplitude ... times the phase ... accounting for the distance to the detector.

		t = sqrt(T)

		t_fund =    (t * fund_tran_ratio) * (fund_tran_amp / np.abs(fund_tran_amp))	* (np.exp(-2j*pi * LT * n_eff_fund  / wavelength))
		# The amplitude ... times the phase ... accounting for the distance to the detector.

		t_first =   (t * first_order_tran_ratio) * (first_order_tran_amp / np.abs(first_order_tran_amp)) * (np.exp(-2j*pi * LT * n_eff_first  / wavelength))
		#The amplitude ... times the phase ... accounting for the distance to the detector.

		if mode == 0:
			r00 = r_fund
			r01 = r_first

			t00 = t_fund
			t01 = t_first

			su0 = sqrt(Su)
			sd0 = sqrt(Sd)

			r00s.append(r00[0])
			r01s.append(r01[0])
			t00s.append(t00[0])
			t01s.append(t01[0])
			su0s.append(su0)
			sd0s.append(sd0)
		else:
			r10 = r_fund
			r11 = r_first

			t10 = t_fund
			t11 = t_first

			su1 = sqrt(Su)
			sd1 = sqrt(Sd)

			r10s.append(r10[0])
			r11s.append(r11[0])
			t10s.append(t10[0])
			t11s.append(t11[0])
			su1s.append(su1)
			sd1s.append(sd1)

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

		if mode == 0:
			ws.append(w*1000)
			Rs.append(R*100)
			Ts.append(T*100)
			Ss.append(S*100)
			NET_LIST.append(NET)
			Sus.append(Su*100)
			Sds.append(Sd*100)
			NET_LOSS_LIST.append(NET_LOSS)
			norm_Sus.append(norm_Su*100)

		if mode == 0:

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
			f1.write("\n \n")

			f1.write("FUNDAMENTAL MODE \n")
			f1.write("n_eff:   %s \n" % (n_eff_fund[0]))
			f1.write("Re(r00): %s \n" % (np.real(r00))[0])
			f1.write("Im(r00): %s \n" % (np.imag(r00))[0])
			f1.write("Re(r01): %s \n" % (np.real(r01))[0])
			f1.write("Im(r01): %s \n" % (np.imag(r01))[0])
			f1.write("Re(t00): %s \n" % (np.real(t00))[0])
			f1.write("Im(t00): %s \n" % (np.imag(t00))[0])
			f1.write("Re(t01): %s \n" % (np.real(t01))[0])
			f1.write("Im(t01): %s \n" % (np.imag(t01))[0])
			f1.write("Re(su0): %s \n" % (np.real(su0)))
			f1.write("Im(su0): %s \n" % (np.imag(su0)))
			f1.write("Re(sd0): %s \n" % (np.real(sd0)))
			f1.write("Im(sd0): %s \n" % (np.imag(sd0)))
			f1.write("\n")

		else:

			f1.write("FIRST ORDER MODE \n")
			f1.write("n_eff:   %s \n" % (n_eff_first[0]))
			f1.write("Re(r10): %s \n" % (np.real(r10))[0])
			f1.write("Im(r10): %s \n" % (np.imag(r10))[0])
			f1.write("Re(r11): %s \n" % (np.real(r11))[0])
			f1.write("Im(r11): %s \n" % (np.imag(r11))[0])
			f1.write("Re(t10): %s \n" % (np.real(t10))[0])
			f1.write("Im(t10): %s \n" % (np.imag(t10))[0])
			f1.write("Re(t11): %s \n" % (np.real(t11))[0])
			f1.write("Im(t11): %s \n" % (np.imag(t11))[0])
			f1.write("Re(su1): %s \n" % (np.real(su1)))
			f1.write("Im(su1): %s \n" % (np.imag(su1)))
			f1.write("Re(sd1): %s \n" % (np.real(sd1)))
			f1.write("Im(sd1): %s \n" % (np.imag(sd1)))

			f1.write("--------------------------------------------------- \n")

		sim.reset_meep()

	#-------------------------------------------------------------

notch(0)

old_widths = range(4, 32, 2)

widths = [elt/100 for elt in old_widths]

p0 = t00s[0] * np.exp(n_eff_funds[0]  * 2j * pi * (Ls + widths/2) / wavelength)
p1 = t11s[0] * np.exp(n_eff_firsts[0] * 2j * pi * (Ls + widths/2) / wavelength)

ws = []
Rs = []
Ts = []
Ss = []
NET_LIST = []
Sus = []
Sds = []
NET_LOSS_LIST = []
norm_Sus = []

r00s = []
r01s = []
r10s = []
r11s = []

t00s = []
t01s = []
t10s = []
t11s = []

su0s = []
sd0s = []
su1s = []
sd1s = []

for notch_index in range(4, 12, 2):
	notch_width = notch_index / 100
	notch(notch_width)

def writeme(value, phasecompensation):
	f2.write("%s \n" % np.real(np.divide(value, phasecompensation)))		# 7
	f2.write("%s \n" % np.imag(np.divide(value, phasecompensation)))		# 8

f2.write("%s \n" % (ws))			# 0
f2.write("%s \n" % (Rs))			# 1
f2.write("%s \n" % (Ts))			# 2
f2.write("%s \n" % (Ss))			# 3
f2.write("%s \n" % (Sus))			# 4
f2.write("%s \n" % (Sds))			# 5
f2.write("%s \n" % (norm_Sus))		# 6

# f2.write("%s \n" % (real_r00s))		# 7
writeme(r00s, p0);
writeme(r01s, p0);
writeme(t00s, p0);
writeme(t01s, p0);
writeme(sd0s, p0);
writeme(su0s, p0);
# f2.write("%s \n" % np.real(np.divide(r00s, p0)))		# 7
# f2.write("%s \n" % np.imag(np.divide(r00s, p0)))		# 8
# f2.write("%s \n" % (real_r01s))		# 9
# f2.write("%s \n" % (imag_r01s))		# 10
# f2.write("%s \n" % (real_t00s))		# 11
# f2.write("%s \n" % (imag_t00s))		# 12
# f2.write("%s \n" % (real_t01s))		# 13
# f2.write("%s \n" % (imag_t01s))		# 14
# f2.write("%s \n" % (real_su0s))		# 15
# f2.write("%s \n" % (imag_su0s))		# 16
# f2.write("%s \n" % (real_sd0s))		# 17
# f2.write("%s \n" % (imag_sd0s))		# 18

writeme(r10s, p1);
writeme(r11s, p1);
writeme(t10s, p1);
writeme(t11s, p1);
writeme(sd1s, p1);
writeme(su1s, p1);
# f2.write("%s \n" % (real_r10s))		# 19
# f2.write("%s \n" % (imag_r10s))		# 20
# f2.write("%s \n" % (real_r11s))		# 21
# f2.write("%s \n" % (imag_r11s))		# 22
# f2.write("%s \n" % (real_t10s))		# 23
# f2.write("%s \n" % (imag_t10s))		# 24
# f2.write("%s \n" % (real_t11s))		# 25
# f2.write("%s \n" % (imag_t11s))		# 26
# f2.write("%s \n" % (real_su1s))		# 27
# f2.write("%s \n" % (imag_su1s))		# 28
# f2.write("%s \n" % (real_sd1s))		# 29
# f2.write("%s \n" % (imag_sd1s))		# 30

f2.write("%s \n" % (n_eff_funds))	# 31
f2.write("%s \n" % (n_eff_firsts))	# 32

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
