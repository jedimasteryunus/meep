import os
from meep import mpb
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, tan, sqrt
from scipy.optimize import fsolve

os.system("rm -r notch-out/")

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

a = 8 						# length of cell
h = 0.4 					# height of Waveguide
dpml = .5;
# monitorheight = .8
monitorheight = .6
H = monitorheight + 2*dpml 	#height of cell
resolution = 50

Lr2 = 	2.00-a/2	# Position of reflection monitor2
Ls = 	1.75-a/2	# Position of source
Lr1 = 	1.50-a/2	# Position of reflection monitor1
Lt = 	a/2-1.50	# Position of transmisison monitor
# monitorheight = 4*h

only_fund = True

if True:
	case = "AlN420"
	wavelength = 0.4203
	widths = [0, .09, .1, .11, .2]
	# widths = np.arange(0, .12, .01)

	print(widths)

	n_e = 2.2415
	n_o = 2.1888
	n_x = n_o
	n_y = n_e
	n_z = n_o

	n_u = 1.4681
	n_l = 1.7800
	# n_l = 1.4681

	e = .5

	h = .2
	hu = .2
	hl = 100

	# only_fund = False

	ang = 90
else:
	case = "LN637"
	wavelength = 0.6372
	widths = [0, .09, .1, .11]

	n_e = 2.2012
	n_o = 2.3048
	n_x = n_o
	n_y = n_o
	n_z = n_e

	n_u = 1.4569
	n_l = 1.4569

	e = .4

	h = .2
	hu = 100
	hl = 100

	ang = 80

f1 = open('notch-' + case + '.out', 'w')
f2 = open('notch-' + case + '.txt', 'w')

hu = min(hu, (H-h)/2)
hl = min(hl, (H-h)/2)

n_m = max(n_x, n_y, n_z)
n_c = max(n_l, n_u)

fcen = 1 / wavelength
df = 0.05

default_material = 	mp.Medium(epsilon 		= 1)
upper_material = 	mp.Medium(epsilon 		= n_u ** 2)
core_material = 	mp.Medium(epsilon_diag 	= mp.Vector3(n_x ** 2, n_y ** 2, n_z ** 2))
lower_material = 	mp.Medium(epsilon 		= n_l ** 2)

geometry_lattice = mp.Lattice(size=mp.Vector3(0, monitorheight, 0))

# bandNum = 4

# geometry = [mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material)]

geometrympb = [	mp.Block(mp.Vector3(mp.inf, monitorheight, 	mp.inf), center = mp.Vector3(0,0), 		material = default_material),
				mp.Block(mp.Vector3(mp.inf, hu+h/2, 		mp.inf), center = mp.Vector3(0,  (hu+h/2)/2), 	material = upper_material),
				mp.Block(mp.Vector3(mp.inf, hl+h/2, 		mp.inf), center = mp.Vector3(0, -(hl+h/2)/2), 	material = lower_material),
				mp.Block(mp.Vector3(mp.inf, h, 				mp.inf), center = mp.Vector3(0, 0), 			material = core_material)]

ms = mpb.ModeSolver(
	geometry_lattice=geometry_lattice,
	resolution=resolution,
	# ensure_periodicity=False,
	geometry=geometrympb
)

band = 1;

ms.find_k(
	# p              = mp.EVEN_Z,               # Polarization
	p              = mp.ODD_Z,               # Polarization
	omega          = fcen,              # Omega to find corresponding k
	band_min       = band,                  # Minimum band index
	band_max       = band,                  # Max band index
	korig_and_kdir = mp.Vector3(1,0,0),      # K vector orientation
	tol            = 1e-7,               # solver tolerance
	kmag_guess     = fcen * n_m,       # initial guess
	kmag_min       = fcen * n_c,        # Minimum
	kmag_max       = fcen * n_m)         # Maximum

E1 = np.squeeze(ms.get_efield(band, bloch_phase=False))
H1 = np.squeeze(ms.get_hfield(band, bloch_phase=False))
# E2 = np.squeeze(ms.get_hfield(3, bloch_phase=False))
# H2 = np.squeeze(ms.get_efield(3, bloch_phase=False))

# Ie_true = 	abs(E1[:,0])**2 + abs(E1[:,1])**2 + abs(E1[:,2])**2
# Ih_true = 	abs(H1[:,0])**2 + abs(H1[:,1])**2 + abs(H1[:,2])**2
# Ie_true2 = 	abs(E2[:,0])**2 + abs(E2[:,1])**2 + abs(E2[:,2])**2
# Ih_true2 = 	abs(H2[:,0])**2 + abs(H2[:,1])**2 + abs(H2[:,2])**2

E0 = np.abs(E1[:,2])
H0 = np.abs(H1[:,2])

# y_list = np.arange(-monitorheight/2-.5/resolution, monitorheight/2+1.5/resolution, 1/resolution)
y_list = np.arange(-monitorheight/2+.5/resolution, monitorheight/2+.5/resolution, 1/resolution)

# plt.plot(y_list, np.real(E1[:,2]), 'b-', label='E1(2)')
# plt.plot(y_list, np.real(H1[:,1]), 'r-', label='H1(1)')
# plt.plot(y_list, np.real(E1[:,1]), 'go-', label='E1(1)')
# plt.plot(y_list, np.real(H1[:,2]), 'co-', label='H1(2)')
# # plt.axis([40.0, 300.0, 0.0, 100.0])
# plt.xlabel("y (um)")
# plt.ylabel("Field (a.u.)")
# plt.legend(loc="center right")
# plt.show()
#
# print(y_list)
# print(y_list.size)
#
# quit()

def notch(w):

	print("#----------------------------------------")
	print("NOTCH WIDTH: %s nanometers" % (w * 1000))
	print("#----------------------------------------")

	# w is the width of the notch in the waveguide

	angrad = ang * pi / 180
	bottomoffset = e * h / tan(angrad)

	vertices = [mp.Vector3( w/2 + bottomoffset, (.5 + e)*h),
				mp.Vector3(-w/2 - bottomoffset, (.5 + e)*h),
				mp.Vector3(-w/2 + bottomoffset, (.5 - e)*h),
				mp.Vector3( w/2 - bottomoffset, (.5 - e)*h)]

	if bottomoffset > w/2:
		ratio = (w/2)/bottomoffset
		vertices = [mp.Vector3( w/2, h/2),
					mp.Vector3(-w/2, h/2),
					mp.Vector3(0, (.5 - e*ratio)*h)]

	print(vertices)

	#Waveguide Geometry
	cell = mp.Vector3(a, H)

	geometry = [mp.Block(cell, 					center = mp.Vector3(0,0), 				material = default_material),
				mp.Block(mp.Vector3(a, hu+h/2), center = mp.Vector3(0,  (hu+h/2)/2), 	material = upper_material),
				mp.Block(mp.Vector3(a, hl+h/2), center = mp.Vector3(0, -(hl+h/2)/2), 	material = lower_material),
				mp.Block(mp.Vector3(a, h), 		center = mp.Vector3(0,0), 				material = core_material)]

	if w > 0:
		geometry.append(mp.Prism(vertices, height=1, material = upper_material))

	pml_layers = [mp.Absorber(thickness=dpml)]

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

	modes = [0, 1]

	if only_fund:
		modes = [0];

	# eig_parity_fund = 	mp.EVEN_Z+mp.EVEN_Y;
	# eig_parity_first = 	mp.EVEN_Z+mp.ODD_Y;

	eig_parity_fund = 	mp.EVEN_Y;
	eig_parity_first = 	mp.ODD_Y;

	# for mode in [0, 1]:
	for mode in modes:

		if mode == 0:
			eig_parity = eig_parity_fund	# Fundamental
			print("-----------")
			print("MODE TYPE: FUNDAMENTAL")

		else:
			eig_parity = eig_parity_first 	# First Order
			print("-----------")
			print("MODE TYPE: FIRST ORDER")

		sources = [	mp.EigenModeSource(
						mp.GaussianSource(	frequency = fcen,
											fwidth = df),
											size = mp.Vector3(0,H),
											center = mp.Vector3(Ls, 0),
											eig_parity = eig_parity) ]

		sim = mp.Simulation(cell_size = cell,
							boundary_layers = pml_layers,
							geometry = geometry,
							sources = sources,
							resolution = resolution,
							force_complex_fields = True)
		'''
		#--------------------------------------------------
		#FOR DISPLAYING THE GEOMETRY

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
		# FOR GENERATING THE TRANSMITTANCE SPECTRUM

		nfreq = 1  # number of frequencies at which to compute flux

		refl_fr1 = 	mp.FluxRegion(center=mp.Vector3(Lr1,0), 	size=mp.Vector3(0,monitorheight))	# Reflected flux 1
		refl_fr2 = 	mp.FluxRegion(center=mp.Vector3(Lr2,0), 	size=mp.Vector3(0,monitorheight))	# Reflected flux 2
		tran_fr = 	mp.FluxRegion(center=mp.Vector3(Lt,0), 		size=mp.Vector3(0,monitorheight))	# Transmitted flux
		su_fr = 	mp.FluxRegion(center=mp.Vector3(0, monitorheight/2),	size=mp.Vector3(a,0))	# Flux loss above the waveguide
		sd_fr = 	mp.FluxRegion(center=mp.Vector3(0,-monitorheight/2), 	size=mp.Vector3(a,0))	# Flux loss below the waveguide

		refl1 = sim.add_flux(fcen, df, nfreq, refl_fr1)
		refl2 = sim.add_flux(fcen, df, nfreq, refl_fr2)
		tran = 	sim.add_flux(fcen, df, nfreq, tran_fr)
		su = 	sim.add_flux(fcen, df, nfreq, su_fr)
		sd = 	sim.add_flux(fcen, df, nfreq, sd_fr)

		# ------------------------ CODE FOR SEPARATING FUND AND FIRST ORDER MODE STARTS HERE ------------------------


		refl_vals = []
		tran_vals = []

		def get_refl_slice(sim):
			# print(sim.get_array(center=mp.Vector3(Lr1,0), size=mp.Vector3(0,H), component=mp.Ez, cmplx=True))
			# refl_val = sim.get_array(center=mp.Vector3(Lr1,0), size=mp.Vector3(0,H), component=mp.Ez, cmplx=True)
			refl_vals.append(sim.get_array(center=mp.Vector3(Lr1,0), size=mp.Vector3(0,monitorheight-2/resolution), component=mp.Ez, cmplx=True))


		def get_tran_slice(sim):
			# print(sim.get_array(center=mp.Vector3(Lt, 0), size=mp.Vector3(0,H), component=mp.Ez, cmplx=True))
			# tran_val = sim.get_array(center=mp.Vector3(Lt, 0), size=mp.Vector3(0,H), component=mp.Ez, cmplx=True)
			tran_vals.append(sim.get_array(center=mp.Vector3(Lt, 0), size=mp.Vector3(0,monitorheight-2/resolution), component=mp.Ez, cmplx=True))

		gif = True

		if gif: # and w == 0.1:
			sim.use_output_directory()
			sim.run(mp.at_beginning(mp.output_epsilon),
					mp.at_end(get_refl_slice),
					mp.at_end(get_tran_slice),
					until=100)
			# sim.run(mp.at_every(wavelength / 20, mp.output_efield_z), until=wavelength)
			# sim.run(mp.at_every(wavelength/20 , mp.output_png(mp.Ez, "-RZc bluered -A notch-out/notch-eps-000000000.h5 -a gray:.2")), until=19*wavelength/20)
			sim.run(mp.at_every(wavelength/20 , mp.output_png(mp.Ez, "-RZc bluered -A notch-out/notch-eps-000000.00.h5 -a gray:.2")), until=19*wavelength/20)
			sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(), 1e-5))
		else:
			sim.run(mp.at_beginning(mp.output_epsilon),
					mp.at_time(100, get_refl_slice),
					mp.at_time(100, get_tran_slice),
					until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(), 1e-5))

		os.system("h5topng notch-out/notch-eps-000000.00.h5; mv notch-out/notch-eps-000000.00.png " 	+ case + "-" + str(int(w*1000)) + "-" + str(mode) + "-eps.png")
		os.system("cp notch-out/notch-ez-000100.00.png " 	+ case + "-" + str(int(w*1000)) + "-" + str(mode) + ".png")
		os.system("convert notch-out/notch-ez-*.png    " 	+ case + "-" + str(int(w*1000)) + "-" + str(mode) + ".gif")

		# get_eigenmode(fcen, mp., refl_fr1, 1, kpoint)
		# v = mp.volume(mp.vec(Lr1, -monitorheight/2), mp.vec(Lr1, monitorheight/2))
		# mode = get_eigenmode(fcen, mp.X, v, v, 1, mp.vec(0, 0, 0), True, 0, 0, 1e-7, True)
		# print(mode.amplitude)

		# coef_refl_fund = 	mp.get_eigenmode_coefficients_and_kpoints(refl_fr1, 1, eig_parity=eig_parity_fund, 	eig_resolution=resolution, eig_tolerance=1e-7)
		# coef_refl_first = 	mp.get_eigenmode_coefficients_and_kpoints(refl_fr1, 1, eig_parity=eig_parity_first, eig_resolution=resolution, eig_tolerance=1e-7)
		#
		# coef_tran_fund = 	mp.get_eigenmode_coefficients_and_kpoints(tran_fr, 	1, eig_parity=eig_parity_fund, 	eig_resolution=resolution, eig_tolerance=1e-7)
		# coef_tran_first = 	mp.get_eigenmode_coefficients_and_kpoints(tran_fr, 	1, eig_parity=eig_parity_first, eig_resolution=resolution, eig_tolerance=1e-7)

		ep = mp.ODD_Z

		coef_refl_fund, 	vgrp, kpoints_fund 	= 	sim.get_eigenmode_coefficients(refl1, [1], eig_parity=ep, 	eig_resolution=resolution, eig_tolerance=1e-7)
		coef_refl_first, 	vgrp, kpoints_first =	sim.get_eigenmode_coefficients(refl1, [2], eig_parity=ep,	eig_resolution=resolution, eig_tolerance=1e-7)

		coef_tran_fund, 	vgrp, kpoints_fund 	=	sim.get_eigenmode_coefficients(tran,  [1], eig_parity=ep, 	eig_resolution=resolution, eig_tolerance=1e-7)
		coef_tran_first,	vgrp, kpoints_first = 	sim.get_eigenmode_coefficients(tran,  [2], eig_parity=ep, 	eig_resolution=resolution, eig_tolerance=1e-7)

		# print(kpoints_fund)
		# print(kpoints_first)
		# print(kpoints_fund[0])
		# print(type(kpoints_fund[0]))
		# print(dir(kpoints_fund[0]))

		n_eff_fund = 	wavelength*kpoints_fund[0].x
		n_eff_first = 	wavelength*kpoints_first[0].x

		print(n_eff_fund)
		print(n_eff_first)

		# print(coef_refl_fund)
		# print(type(coef_refl_fund))
		# print(dir(coef_refl_fund))

		# print(coef_refl_fund[0])
		# print(coef_refl_fund[0][0,0,:])
		#
		# fund_refl_amp = 		coef_refl_fund[0][0,0,1];
		# first_order_refl_amp =	coef_refl_first[0][0,0,1];
		# fund_tran_amp =			coef_tran_fund[0][0,0,0];
		# first_order_tran_amp =	coef_tran_first[0][0,0,0];

		print("get_eigenmode_coefficients:\n")

		print(coef_refl_fund)
		print(coef_refl_first)
		print(coef_tran_fund)
		print(coef_tran_first)

		print("\n")
		# print(coef_refl_fund[0,0,:])

		fund_refl_amp = 		coef_refl_fund[0,0,1];
		first_order_refl_amp =	coef_refl_first[0,0,1];
		fund_tran_amp =			coef_tran_fund[0,0,0];
		first_order_tran_amp =	coef_tran_first[0,0,0];

		refl_val = refl_vals[0]
		tran_val = tran_vals[0]

		# n_eff must satisfy n_e <= n_eff <= n_c for the mode to be bound.

		def fund_func(n_eff):
			if n_eff >= n_c and n_eff <= n_e:
				return sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2))

		def first_order_func(n_eff):
			if n_eff >= n_c and n_eff <= n_e:
				return sqrt(n_eff**2 - n_c**2) - sqrt(n_e**2 - n_eff**2) * tan(pi * h / wavelength * sqrt(n_e**2 - n_eff**2) - pi / 2)

		initial_guess = (n_c + n_e) / 2

		# n_eff_fund = 	fsolve(fund_func, initial_guess)
		# n_eff_first = 	fsolve(first_order_func, n_c)

		print(n_eff_fund, n_eff_first)

		assert(n_eff_fund > n_eff_first)

		if len(n_eff_funds) == 0:
			n_eff_funds.append(n_eff_fund)

		if len(n_eff_firsts) == 0:
			n_eff_firsts.append(n_eff_first)

		ky0_fund =  	np.abs(2 * pi / wavelength * sqrt(n_e**2 - n_eff_fund **2))
		ky0_first = 	np.abs(2 * pi / wavelength * sqrt(n_e**2 - n_eff_first**2))

		ky1_fund =  	ky0_fund # np.abs(2 * pi / wavelength * sqrt(n_eff_fund **2 - n_c**2))
		ky1_first = 	ky0_first # np.abs(2 * pi / wavelength * sqrt(n_eff_first**2 - n_c**2))

		E_fund = 		lambda y : cos(ky0_fund  * y) if np.abs(y) < h / 2 else cos(ky0_fund  * h/2) * np.exp(-ky1_fund  * (np.abs(y) - h / 2))
		E_first_order = lambda y : sin(ky0_first * y) if np.abs(y) < h / 2 else sin(ky0_first * h/2) * np.exp(-ky1_first * (np.abs(y) - h / 2)) * np.sign(y)
		# y_list = np.arange(-H/2+.5/resolution, H/2-.5/resolution, 1/resolution)

		#print("Y LIST: ", y_list)
		#print("SIZE OF Y LIST: ", y_list.size)

		E_fund_vec = np.zeros(y_list.size)
		E_first_order_vec = np.zeros(y_list.size)

		for index in range(y_list.size):
			y = y_list[index]
			E_fund_vec[index] = E_fund(y)
			E_first_order_vec[index] = E_first_order(y)

		# print(dir(sim))
		# print(type(sim.get_eigenmode_coefficients))
		# print(dir(sim.get_eigenmode_coefficients))
		# print(type(sim.get_eigenmode))
		# print(dir(sim.get_eigenmode))
		# print(sim.get_eigenmode.__code__.co_varnames)
		# print(sim.get_eigenmode.__defaults__)

		# E1 = sim.get_eigenmode(fcen, mp.X, refl1.where, 1, None)
		# E2 = sim.get_eigenmode(fcen, mp.X, refl1.where, 2, None)
		# E3 = sim.get_eigenmode(fcen, mp.X, refl1.where, 3, None)
		# E4 = sim.get_eigenmode(fcen, mp.X, refl1.where, 4, None)

		# E1 = sim.get_eigenmode(fcen, mp.X, refl1.where, 1, mp.Vector3(0, 0, 0))
		# E2 = sim.get_eigenmode(fcen, mp.X, refl1.where, 2, mp.Vector3(0, 0, 0))
		# E3 = sim.get_eigenmode(fcen, mp.X, refl1.where, 3, mp.Vector3(0, 0, 0))
		# E4 = sim.get_eigenmode(fcen, mp.X, refl1.where, 4, mp.Vector3(0, 0, 0))

		# print(refl1.where)
		# numEA = 0
		# E1 = sim.fields.get_eigenmode(fcen, mp.X, refl1.where, refl1.where, 1, mp.vec(0,0,0), True, 0, numEA, 1e-7, True)

		# print(type(E1))
		# print(dir(E1))
		#
		# print(refl1.where)
		# # print(E1, E2, E3, E4)
		# print(E1.amplitude, E1.band_num, E1.group_velocity, E1.k)
		# print(type(E1.amplitude))
		# print(dir(E1.amplitude))
		# print(doc(E1.amplitude))
		# print(self(E1.amplitude))
		# print(E1.amplitude(y_list))
		# print(type(E1.amplitude(y_list)))
		# print(dir(E1.amplitude(y_list)))
		# print(E1.amplitude, E2.amplitude, E3.amplitude, E4.amplitude)

		# E1 = sim.get_eigenmode(fcen, mp.X, refl1.where, refl1.where, 1, mp.vec(0, 0, 0))
		# E2 = sim.get_eigenmode(fcen, mp.X, refl1.where, refl1.where, 2, mp.vec(0, 0, 0))
		# E3 = sim.get_eigenmode(fcen, mp.X, refl1.where, refl1.where, 3, mp.vec(0, 0, 0))
		# E4 = sim.get_eigenmode(fcen, mp.X, refl1.where, refl1.where, 4, mp.vec(0, 0, 0))

		# print(y_list)
		#
		# E1a = np.zeros(y_list.size, dtype='Complex128');
		# E2a = np.zeros(y_list.size, dtype='Complex128');
		# E3a = np.zeros(y_list.size, dtype='Complex128');
		# E4a = np.zeros(y_list.size, dtype='Complex128');
		#
		# # print(mp.eigenmode_amplitude.__code__.co_varnames)
		# # print(mp.eigenmode_amplitude.__defaults__)
		#
		# for i in range(y_list.size):
		# 	# print(E1)
		# 	# print(E1.swigobj)
		# 	# print(E1.amplitude)
		# 	# print(mp.vec(Lr1, y_list[i], 0))
		# 	# print(mp.Vector3(Lr1, y_list[i], 0))
		# 	# print(mp.Ez)
		# 	# E1a[i] = mp.eigenmode_amplitude(E1.swigobj, mp.vec(Lr1, y_list[i], 0), mp.Ez);
		# 	eval_point = mp.vec(Lr1, y_list[i], 0)
		# 	# print(eval_point)
		# 	E1a[i] = mp.eigenmode_amplitude(E1.swigobj, eval_point, mp.Ex)
		# 	E2a[i] = mp.eigenmode_amplitude(E2.swigobj, eval_point, mp.Ex)
		# 	E3a[i] = mp.eigenmode_amplitude(E3.swigobj, eval_point, mp.Ex)
		# 	E4a[i] = mp.eigenmode_amplitude(E4.swigobj, eval_point, mp.Ex)
		# 	# E1a[i] = mp.eigenmode_amplitude(E1.amplitude, mp.Vector3(Lr1, y_list[i], 0), mp.Ez);
		#
		# plt.plot(y_list, np.abs(E1a)**2, 'bo-', label='E1')
		# plt.plot(y_list, np.abs(E2a)**2, 'ro-', label='E2')
		# plt.plot(y_list, np.abs(E3a)**2, 'go-', label='E3')
		# plt.plot(y_list, np.abs(E4a)**2, 'co-', label='E4')
		# # plt.axis([40.0, 300.0, 0.0, 100.0])
		# plt.xlabel("y (um)")
		# plt.ylabel("Field (a.u.)")
		# plt.legend(loc="center right")
		# plt.show()


		# print("r VECTOR: ", 	refl_val)
		# print("t VECTOR: ", 	tran_val)
		# print("E0 VECTOR: ", 	E_fund_vec)
		# print("E1 VECTOR: ", 	E_first_order_vec)
		#

		# fund_refl_amp_2 = 			np.conj(np.dot(refl_val, E_fund_vec) 		/ np.dot(E_fund_vec, E_fund_vec))					# Conjugate becasue MEEP uses physics exp(kz-wt) rather than engineering exp(wt-kz)
		# first_order_refl_amp_2 = 	np.conj(np.dot(refl_val, E_first_order_vec) / np.dot(E_first_order_vec, E_first_order_vec))
		# fund_tran_amp_2 = 			np.conj(np.dot(tran_val, E_fund_vec) 		/ np.dot(E_fund_vec, E_fund_vec))
		# first_order_tran_amp_2 = 	np.conj(np.dot(tran_val, E_first_order_vec) / np.dot(E_first_order_vec, E_first_order_vec))

		fund_refl_amp = 			np.conj(np.dot(refl_val, E0) / np.dot(E0, E0))					# Conjugate becasue MEEP uses physics exp(kz-wt) rather than engineering exp(wt-kz)
		first_order_refl_amp = 		0 # np.conj(np.dot(refl_val, E1[:,2]) / np.dot(E1[:,2], E1[:,2]))
		fund_tran_amp = 			np.conj(np.dot(tran_val, E0) / np.dot(E0, E0))
		first_order_tran_amp = 		0 # np.conj(np.dot(tran_val, E1[:,2]) / np.dot(E1[:,2], E1[:,2]))

		fund_refl = 				np.conj(fund_refl_amp)*E0
		not_fund_refl = 			refl_val - fund_refl
		fund_tran = 				np.conj(fund_tran_amp)*E0
		not_fund_tran = 			tran_val - fund_tran

		fund_refl_power =			np.dot(np.conj(fund_refl),		fund_refl)
		not_fund_refl_power = 		np.dot(np.conj(not_fund_refl), 	not_fund_refl)
		fund_tran_power =			np.dot(np.conj(fund_tran),		fund_tran)
		not_fund_tran_power = 		np.dot(np.conj(not_fund_tran), 	not_fund_tran)

		fund_refl_ratio = 			np.abs(fund_refl_power 		/ (fund_refl_power + not_fund_refl_power))
		first_order_refl_ratio = 	0
		fund_tran_ratio = 			np.abs(fund_tran_power 		/ (fund_tran_power + not_fund_tran_power))
		first_order_tran_ratio = 	0

		# plt.plot(y_list, np.abs(refl_val),	'bo-',label='reflectance')
		# plt.plot(y_list, np.abs(tran_val),	'ro-',label='transmittance')
		# plt.plot(y_list, E0,	'go-',label='E0')
		# plt.plot(y_list, fund_refl_amp*E0,	'co-',label='over')
		# # plt.axis([40.0, 300.0, 0.0, 100.0])
		# plt.xlabel("y (um)")
		# plt.ylabel("Field")
		# plt.legend(loc="center right")
		# plt.show()
		#
		# print("\n")
		#
		# print(tran_val.size, refl_val.size)
		#
		# print("\n")
		#
		# print(tran_val, refl_val, E0)
		#
		# print("\n")
		# # print(np.conj(tran_val), tran_val, E1[:,2])
		# #
		# # print(np.conj(refl_val), refl_val, E1[:,2])
		#
		# refl_tot_power = np.abs(np.dot(np.conj(refl_val), refl_val))
		# tran_tot_power = np.abs(np.dot(np.conj(tran_val), tran_val))
		#
		# print(fund_refl_amp, refl_tot_power, fund_tran_amp, tran_tot_power)

		# print(fund_refl_amp 		, fund_refl_amp_2 		)
		# print(first_order_refl_amp 	, first_order_refl_amp_2)
		# print(fund_tran_amp 		, fund_tran_amp_2 		)
		# print(first_order_tran_amp 	, first_order_tran_amp_2)
		#
		# print(np.angle(fund_refl_amp), 			np.angle(fund_refl_amp_2))
		# print(np.angle(first_order_refl_amp), 	np.angle(first_order_refl_amp_2))
		# print(np.angle(fund_tran_amp), 			np.angle(fund_tran_amp_2))
		# print(np.angle(first_order_tran_amp), 	np.angle(first_order_tran_amp_2))

		# fund_refl_power = 			np.abs(fund_tran_amp)			** 2
		# fund_refl_power = 			np.abs(fund_refl_amp)			** 2
		# first_order_refl_power = 	np.abs(first_order_refl_amp) 	** 2
		# fund_tran_power = 			np.abs(fund_tran_amp) 			** 2
		# first_order_tran_power = 	np.abs(first_order_tran_amp) 	** 2

		# print(fund_refl_power, first_order_refl_power, fund_tran_power, first_order_tran_power)

		# fund_refl_ratio = 			fund_refl_power 		/ (fund_refl_power + first_order_refl_power)
		# first_order_refl_ratio = 	first_order_refl_power 	/ (fund_refl_power + first_order_refl_power)
		# fund_tran_ratio = 			fund_tran_power 		/ (fund_tran_power + first_order_tran_power)
		# first_order_tran_ratio = 	first_order_tran_power 	/ (fund_tran_power + first_order_tran_power)

		# fund_refl_power = 			np.abs(fund_refl_amp)			** 2
		# first_order_refl_power = 	np.abs(first_order_refl_amp)	** 2
		# fund_tran_power = 			np.abs(fund_tran_amp)			** 2
		# first_order_tran_power = 	np.abs(first_order_tran_amp)	** 2
		#
		# fund_refl_ratio = 			fund_refl_power 		/ refl_tot_power
		# first_order_refl_ratio = 	first_order_refl_power 	/ refl_tot_power
		# fund_tran_ratio = 			fund_tran_power 		/ tran_tot_power
		# first_order_tran_ratio = 	first_order_tran_power 	/ tran_tot_power
		#
		# fund_refl_ratio = 			1
		# first_order_refl_ratio = 	0
		# fund_tran_ratio = 			1
		# first_order_tran_ratio = 	0

		print("Percentage of reflected light in fundamental mode: ", 	fund_refl_ratio * 100)
		print("Percentage of reflected light in first order mode: ", 	first_order_refl_ratio * 100)
		print("Percentage of transmitted light in fundamental mode: ", 	fund_tran_ratio * 100)
		print("Percentage of transmitted light in first order mode: ", 	first_order_tran_ratio * 100)

		# ------------------------ CODE FOR SEPARATING FUND AND FIRST ORDER MODE ENDS HERE ------------------------

		wl = [] #list of wavelengths

		refl1_flux = 	mp.get_fluxes(refl1)
		refl2_flux = 	mp.get_fluxes(refl2)
		tran_flux = 	mp.get_fluxes(tran)
		su_flux = 		mp.get_fluxes(su)
		sd_flux = 		mp.get_fluxes(sd)

		flux_freqs = mp.get_flux_freqs(refl1)

		for i in range(nfreq):
			wl = np.append(wl, 1/flux_freqs[i])
			print(1/flux_freqs[i])

		# for ind, elt in enumerate(wl):
		#     #print(round(elt, 4))
		#     if round(elt, 3) == 0.637:
		#         #print("ALERT: MATCH FOUND")
		#         index = ind

		index = 0

		rp = refl1_flux[index]
		tp = tran_flux[index]

		# print("rp/tp:\n")
		# print(rp)
		# print(tp)
		# print("\n")

		R = 	-refl1_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		T = 	 tran_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		S = 	(refl2_flux[index] - tran_flux[index]) 	/ (refl2_flux[index] - refl1_flux[index])
		Su = 	 su_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
		Sd = 	-sd_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])

		S_correction = (1 - R - T) / (Su + Sd)

		# print(R, T, S, Su, Sd)

		r = sqrt(R)
		t = sqrt(T)

		# The amplitude ... times the phase ... accounting for the distance to the detector (Reverse exp(-kz) of phase).
		# r_fund = 	(r * fund_refl_ratio) 			* (fund_refl_amp 		/ np.abs(fund_refl_amp)) 		* (np.exp( 2j*pi * (-Lr1 - w/2) * n_eff_fund  / wavelength))	# Lr1 is negative because it is the (negative) position of the monitor, not the distace from the monitor to the center.
		# r_first = 	(r * first_order_refl_ratio) 	* (first_order_refl_amp / np.abs(first_order_refl_amp)) * (np.exp( 2j*pi * (-Lr1 - w/2) * n_eff_first / wavelength))
		# t_fund =    (t * fund_tran_ratio) 			* (fund_tran_amp 		/ np.abs(fund_tran_amp))		* (np.exp( 2j*pi * ( Lt  - w/2) * n_eff_fund  / wavelength))
		# t_first =   (t * first_order_tran_ratio) 	* (first_order_tran_amp / np.abs(first_order_tran_amp)) * (np.exp( 2j*pi * ( Lt  - w/2) * n_eff_first / wavelength))

		r_fund = 	(r * fund_refl_ratio) 			* np.exp( 1j*np.angle(fund_refl_amp) 		+ 2j*pi * (-Lr1 - w/2) * n_eff_fund  / wavelength)	# Lr1 is negative because it is the (negative) position of the monitor, not the distace from the monitor to the center.
		r_first = 	(r * first_order_refl_ratio) 	* np.exp( 1j*np.angle(first_order_refl_amp) + 2j*pi * (-Lr1 - w/2) * n_eff_first / wavelength)
		t_fund =    (t * fund_tran_ratio) 			* np.exp( 1j*np.angle(fund_tran_amp) 		+ 2j*pi * ( Lt  - w/2) * n_eff_fund  / wavelength)
		t_first =   (t * first_order_tran_ratio) 	* np.exp( 1j*np.angle(first_order_tran_amp) + 2j*pi * ( Lt  - w/2) * n_eff_first / wavelength)

		if mode == 0:
			r00 = r_fund
			r01 = r_first

			t00 = t_fund
			t01 = t_first

			su0 = sqrt(np.abs(Su*S_correction))
			sd0 = sqrt(np.abs(Sd*S_correction))
			# su0 = sqrt(Su)
			# sd0 = sqrt(Sd)

			r00s.append(r00)
			r01s.append(r01)
			t00s.append(t00)
			t01s.append(t01)
			su0s.append(su0)
			sd0s.append(sd0)

			if only_fund:
				r10s.append(0)
				r11s.append(0)
				t10s.append(0)
				t11s.append(1)
				su1s.append(0)
				sd1s.append(0)
		else:
			r10 = r_fund
			r11 = r_first

			t10 = t_fund
			t11 = t_first

			su1 = sqrt(np.abs(Su*S_correction))
			sd1 = sqrt(np.abs(Sd*S_correction))
			# su1 = sqrt(Su)
			# sd1 = sqrt(Sd)

			r10s.append(r10)
			r11s.append(r11)
			t10s.append(t10)
			t11s.append(t11)
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
			f1.write("Reflection Percentage: %s\n" % (R * 100))
			f1.write("Transmission Percentage: %s\n" % (T * 100))
			f1.write("Total Loss Percentage: %s\n" % (S * 100))
			f1.write("Percentage of Light Accounted For: %s\n" % (NET))
			f1.write("Upper Loss Percentage: %s\n" % (Su * 100))
			f1.write("Lower Loss Percentage: %s\n" % (Sd * 100))
			f1.write("Percentage of Total Loss Accounted For: %s\n" % (NET_LOSS))
			f1.write("Normalized Upper Loss Percentage: %s\n" % (norm_Su * 100))
			f1.write("\n \n")

			f1.write("FUNDAMENTAL MODE \n")
			f1.write("n_eff:   %s\n" % (n_eff_fund))
			f1.write("Re(r00): %s\n" % (np.real(r00)))
			f1.write("Im(r00): %s\n" % (np.imag(r00)))
			f1.write("Re(r01): %s\n" % (np.real(r01)))
			f1.write("Im(r01): %s\n" % (np.imag(r01)))
			f1.write("Re(t00): %s\n" % (np.real(t00)))
			f1.write("Im(t00): %s\n" % (np.imag(t00)))
			f1.write("Re(t01): %s\n" % (np.real(t01)))
			f1.write("Im(t01): %s\n" % (np.imag(t01)))
			f1.write("Re(su0): %s\n" % (np.real(su0)))
			f1.write("Im(su0): %s\n" % (np.imag(su0)))
			f1.write("Re(sd0): %s\n" % (np.real(sd0)))
			f1.write("Im(sd0): %s\n" % (np.imag(sd0)))
			f1.write("\n")

		else:

			f1.write("FIRST ORDER MODE \n")
			f1.write("n_eff:   %s\n" % (n_eff_first))
			f1.write("Re(r10): %s\n" % (np.real(r10)))
			f1.write("Im(r10): %s\n" % (np.imag(r10)))
			f1.write("Re(r11): %s\n" % (np.real(r11)))
			f1.write("Im(r11): %s\n" % (np.imag(r11)))
			f1.write("Re(t10): %s\n" % (np.real(t10)))
			f1.write("Im(t10): %s\n" % (np.imag(t10)))
			f1.write("Re(t11): %s\n" % (np.real(t11)))
			f1.write("Im(t11): %s\n" % (np.imag(t11)))
			f1.write("Re(su1): %s\n" % (np.real(su1)))
			f1.write("Im(su1): %s\n" % (np.imag(su1)))
			f1.write("Re(sd1): %s\n" % (np.real(sd1)))
			f1.write("Im(sd1): %s\n" % (np.imag(sd1)))

			f1.write("--------------------------------------------------- \n")

		sim.reset_meep()

	#-------------------------------------------------------------

# notch(0)
# notch(.1)
#
# quit()

notch(0)

# old_widths = range(4, 32, 2)

# widths = [elt/100 for elt in old_widths]
# widths = [0, .1, .2]
# widths = [0, .09, .1, .11, .2]
# widths = [0, .05, .1, .15, .2, .25, .3]
# widths = [0, .02, .04, .06, .08, .1, .12, .14, .16, .18, .2]

# p0 = t00s[0] * np.exp(n_eff_funds[0]  * 2j * pi * (Lt + widths/2) / wavelength)
# p1 = t11s[0] * np.exp(n_eff_firsts[0] * 2j * pi * (Lt + widths/2) / wavelength)

# p0t = t00s[0]
# p1t = t11s[0]
p0t = [ t00s[0] * np.exp( n_eff_funds[0]  * 2j * pi * (width/2) / wavelength) for width in widths ]
p1t = [ t11s[0] * np.exp( n_eff_firsts[0] * 2j * pi * (width/2) / wavelength) for width in widths ]


# p0t = [ t00s[0] * np.exp( n_eff_funds[0]  * 2j * pi * width / wavelength) for width in widths ]
# p1t = [ t11s[0] * np.exp( n_eff_firsts[0] * 2j * pi * width / wavelength) for width in widths ]

# p0r = [ t00s[0] * np.exp( n_eff_funds[0]  * 2j * pi * (Lt + width/2) / wavelength) for width in widths ]
# p1r = [ t11s[0] * np.exp( n_eff_firsts[0] * 2j * pi * (Lt + width/2) / wavelength) for width in widths ]

p0r = [ t00s[0] * np.exp( n_eff_funds[0]  * 2j * pi * (width/2) / wavelength) for width in widths ]
p1r = [ t11s[0] * np.exp( n_eff_firsts[0] * 2j * pi * (width/2) / wavelength) for width in widths ]

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

# notch(.1)

for width in widths:
	notch(width)

def writeme(value, phasecompensation):
	# print("%s\n" % np.real(np.divide(value, phasecompensation)))
	# print("%s\n" % np.imag(np.divide(value, phasecompensation)))
	# print("%s\n" % (np.real(np.divide(value, phasecompensation))))
	# print("%s\n" % (np.imag(np.divide(value, phasecompensation))))
	f2.write("%s\n" % list(np.real(np.divide(value, phasecompensation))))		# 7
	f2.write("%s\n" % list(np.imag(np.divide(value, phasecompensation))))		# 8

f2.write("%s\n" % (ws))			# 0
f2.write("%s\n" % (Rs))			# 1
f2.write("%s\n" % (Ts))			# 2
f2.write("%s\n" % (Ss))			# 3
f2.write("%s\n" % (Sus))			# 4
f2.write("%s\n" % (Sds))			# 5
f2.write("%s\n" % (norm_Sus))		# 6

writeme(r00s, p0r);
writeme(r01s, p0r);
writeme(t00s, p0t);
print("%s\n" % list(np.real(t00s)))
print("%s\n" % list(np.imag(t00s)))
print("%s\n" % (np.real(p0t)))
print("%s\n" % (np.imag(p0t)))
print("%s\n" % list(np.real(np.divide(t00s, p0t))))
print("%s\n" % list(np.imag(np.divide(t00s, p0t))))
writeme(t01s, p0t);
writeme(sd0s, p0r);
writeme(su0s, p0r);

writeme(r10s, p1r);
writeme(r11s, p1r);
writeme(t10s, p1t);
writeme(t11s, p1t);
writeme(sd1s, p1r);
writeme(su1s, p1r);

f2.write("%s\n" % (n_eff_funds))	# 31
f2.write("%s\n" % (n_eff_firsts))	# 32

f1.close()
f2.close()

plt.subplot(2, 1, 1)
plt.plot(ws, Rs,'bo-',label='reflectance')
plt.plot(ws, Ts,'ro-',label='transmittance')
plt.plot(ws, Ss,'go-',label='net loss')
plt.plot(ws, Sus, 'co-', label='upper loss')
plt.plot(ws, Sds, 'mo-', label='lower loss')
plt.plot(ws, norm_Sus, 'ko-', label='normalized upper loss')
plt.axis([0, 300.0, 0.0, 100.0])
plt.xlabel("wavelength (nm)")
plt.ylabel("percentage (%)")
plt.legend(loc="center right")

plt.subplot(2, 1, 2)
plt.plot(ws, list(np.angle(np.divide(r00s, p0r))),'bo-',label='reflectance')
plt.plot(ws, list(np.angle(np.divide(t00s, p0t))),'ro-',label='transmittance')
plt.axis([0, 300.0, -pi, pi])
plt.xlabel("wavelength (nm)")
plt.ylabel("Phase (rad)")
plt.legend(loc="center right")

plt.show()
