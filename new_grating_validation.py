import os
from meep import mpb
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt, cos, sin
import sys

def validation(dw, dl):

	name = "w=%snm-dl=%snm" % (int(dw), int(dl))
	print(name);

	outputdir = "grating_validation-out-" + name;
	print(outputdir);

	os.system("rm -r " + outputdir + "/")
	f1 = open('grating_validation-' + name + '.out', 'w')

	# print("#----------------------------------------")
	# print("NOTCH WIDTH: %s nanometers" % (w * 1000))
	# print("#----------------------------------------")

	#T = 5;
	#T = 10;
	#T = 20;
	T = 40;
	a = 16 						# length of cell
	# a = 20 						# length of cell
	h = 0.4 					# height of Waveguide
	dpml = 1;
	# monitorheight = .8
	monitorheight = .6
	H = monitorheight + 2*dpml 	#height of cell
	#resolution = 10
	resolution = 10
	H = a
	sxy = H - 2*dpml
	sxy2 = H - 3*dpml
	sxy2 = H

	Lr2 = 	2.00-a/2	# Position of reflection monitor2
	Ls = 	1.75-a/2	# Position of source
	Lr1 = 	1.50-a/2	# Position of reflection monitor1
	Lt = 	a/2-1.50	# Position of transmisison monitor
	# monitorheight = 4*h

	only_fund = True

	ALN420_bool = False

	if ALN420_bool:
		case = "AlN420"
		wavelength = 0.4203
		#widths = [0, .09, .1, .11, .2]
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
		#widths = [0, .09, .1, .11]

		n_e = 2.2012
		n_o = 2.3048
		n_x = n_o
		n_y = n_o
		n_z = n_e

		n_u = 1.4569
		n_l = 1.4569

		e = .4

		h = .2
		hu = 6
		hl = 6

		ang = 80

	# f1 = open('notch-' + case + '.out', 'w')
	# f2 = open('notch-' + case + '.txt', 'w')

	hu = min(hu, (sxy2-h)/2)
	hl = min(hl, (sxy2-h)/2)

	n_m = max(n_x, n_y, n_z)
	n_c = max(n_l, n_u)

	fcen = 1 / wavelength
	df = 0.05

	default_material = 	mp.Medium(epsilon 		= 1)
	upper_material = 	mp.Medium(epsilon 		= n_u ** 2)
	core_material = 	mp.Medium(epsilon_diag 	= mp.Vector3(n_x ** 2, n_y ** 2, n_z ** 2))
	lower_material = 	mp.Medium(epsilon 		= n_l ** 2)

	#Waveguide Geometry
	cell = mp.Vector3(a, H)

	geometry = [mp.Block(cell, 						center = mp.Vector3(0,0), 				material = default_material),
				mp.Block(mp.Vector3(sxy2, hu+h/2), 	center = mp.Vector3(0,  (hu+h/2)/2), 	material = upper_material),
				mp.Block(mp.Vector3(sxy2, hl+h/2), 	center = mp.Vector3(0, -(hl+h/2)/2), 	material = lower_material),
				mp.Block(mp.Vector3(sxy2, h), 		center = mp.Vector3(0,0), 				material = core_material)]

	#lengths = [50., 740., 1470., 2170., 2850., 3390., 3780., 4170., 4560., 4950.]
	# lengths = [50.,   810.,  1560.,  1920.,  2630.,  3370.,  4110.,  4690.,  5020.,  5600.]
	# lengths = range(150,12*300,300)
	# lengths = range(150,12*300,350)
	# lengths = [640., 500., 680., 670., 490., 650., 620., 510., 150., 600.]
	# lengths = [650., 695., 675., 630., 405., 710., 350., 625., 490., 600.]
	# lengths = [675., 675., 665., 680., 325., 625., 160., 805., 310., 600.]
	# lengths = [335.0, 375.0, 330.0, 130.0, 325.0, 355.0, 350.0, 340.0, 330.0, 125.0, 325.0, 360.0, 325.0, 390.0, 320.0, 485.0, 385.0, 250.0, 110.0, 300.0]
	# lengths = [490.0, 315.0, 315.0, 340.0, 305.0, 460.0, 385.0, 290.0, 335.0, 480.0, 290.0, 350.0, 475.0, 480.0, 170.0, 305.0, 330.0, 155.0, 325.0, 300.0]
	# lengths = [290, 495, 285, 495, 295, 485, 300, 485, 285, 500, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# lengths = [290, 495, 285, 495, 295, 485, 300, 485, 285, 500, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# lengths = [395, 395, 395, 395, 395, 395, 395, 395, 395, 395, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# lengths = [250, 220, 255, 210, 255, 190, 270, 200, 275, 185, 275, 175, 285, 175, 280, 205, 270, 190, 280, 300]
	# lengths = [260, 220, 240, 220, 270, 205, 275, 190, 245, 200, 290, 170, 290, 185, 280, 190, 255, 185, 300, 300]

	# lengths = [885, 900, 920, 945, 965, 995, 335, 590, 585, 100]
	# lengths = [915, 915, 990, 940, 1000, 815, 605, 405, 850, 100]
	# lengths = [165, 355, 165, 370, 155, 170, 380, 155, 170, 390, 180, 140, 185, 400, 400, 405, 195, 405, 190, 100]
	# lengths = [355, 555, 365, 160, 365, 175, 155, 165, 170, 165, 390, 160, 165, 195, 380, 415, 395, 405, 485, 100];
	# lengths = [160, 155, 160, 165, 160, 175, 150, 170, 395, 400, 400, 410, 400, 380, 380];
	# lengths = [115, 125, 130, 130, 160, 140, 185, 410, 495, 400];
	# lengths = [115, 125, 130, 130, 160, 140, 185];
	# lengths = [110, 340, 130, 135, 160, 165, 375];
	# lengths = [320, 345, 125, 140, 145, 145, 120];
	# lengths = [115, 340, 130, 135, 160, 165, 495];
	# lengths = [115, 340, 130, 135, 160, 165, 500];
	# lengths = [120, 345, 135, 160, 155, 170, 250]
	# lengths = [110, 125, 135, 120, 185, 100]
	# lengths = [110, 125, 130, 125, 405, 165]
	# lengths = [120, 135, 135, 165, 185, 165];
	# lengths = [120, 135, 135, 165, 175, 200, 165];
	# lengths = [110, 285, 265, 220, 100, 195, 165];
	# lengths = [130, 125, 130, 130, 150, 150, 185];
	# lengths = [130, 125, 130, 130, 150, 150, 190];
	# widths =  [50,  100, 100, 100, 100, 100, 100, 100];
	# lengths = [121, 320, 138, 132, 172, 173];
	# widths =  [50,  100, 100, 100, 100, 100, 100];
	# lengths = [122, 320, 136, 130]
	# widths =  [50,  100, 100, 100, 100];
	# lengths = [146, 309, 336, 130, 140, 169, 137];
	# lengths = [143, 313, 328, 135, 132, 167, 165];
	# lengths = [178, 150, 309, 100, 258, 100, 208];
	# lengths = [181, 152, 307, 100, 259, 100, 199];
	# widths =  [50,  50,  100, 100, 100, 100, 100, 100];
	#
	# lengths = [152, 307, 100, 259];
	# widths =  [50,  100, 100, 100, 100];
	#
	# # lengths = [175, 148, 113, 100, 258, 100]
	# lengths = [174, 145, 316, 99,  271, 92,  88]
	# widths =  [50,  50,  100, 100, 100, 100, 100, 100];
	#
	# # lengths = [157, 336, 104, 276, 88, 89]
	# lengths = [158, 335, 105, 275, 90, 87]
	# widths =  [50,  75,  100, 100, 100, 100, 100];
	#
	# lengths = [156, 336, 104, 274, 88]
	# widths =  [50,  75,  100, 100, 100, 100];
	#
	# lengths = [160, 340, 299, 95, 281, 299, 280];
	# widths =  [50,  75,  100, 100, 100, 75, 50, 50];

	input_lengths =  [30, 25, 30, 30, 22, 28, 18, 32, 16, 25]
	lengths = [10 * length for length in input_lengths]
	widths = 11 * [100]

	assert(len(lengths) + 1 == len(widths));

	for i in range(0, len(lengths)):
		lengths[i] += dl;

	for i in range(0, len(widths)):
		widths[i] += dw;

	x = -widths[0]/2000.;

	for i in range(0, len(lengths)):
		x -= (widths[i+1] + lengths[i]) / 1000.

	x /= 2;

	# grate_positions = []
	# w = .050
	# x = 3*dpml - a/2;
	# h2 =
	for i in range(0, len(lengths)):
		geometry.append(mp.Block(mp.Vector3(widths[i]/1000., e*h),
			center = mp.Vector3(x, h * (1 - e)/2),
			material = default_material))
		x += (widths[i+1] + lengths[i]) / 1000.
		# grate_positions.append(x)

	geometry.append(mp.Block(mp.Vector3(widths[len(widths)-1]/1000., e*h),
		center = mp.Vector3(x, h * (1 - e)/2),
		material = default_material))

	# sources = [mp.EigenModeSource(mp.ContinuousSource(frequency = fcen),
	# 							  size = mp.Vector3(0,H),
	# 							  center = mp.Vector3(2*dpml-a/2, 0))]
	sources = [ mp.EigenModeSource(mp.ContinuousSource(frequency = fcen),
				size = mp.Vector3(0,H),
				center = mp.Vector3(-2*x, 0))]

	# pml_layers = [mp.PML(0.2)]
	pml_layers = [mp.Absorber(thickness=dpml)]

	# resolution = 50

	sim = mp.Simulation(cell_size = cell,
						boundary_layers = pml_layers,
						geometry = geometry,
						sources = sources,
						resolution = resolution)

	refl_fr1 = 	mp.FluxRegion(center=mp.Vector3(Lr1,0), 	size=mp.Vector3(0,monitorheight))	# Reflected flux 1
	refl_fr2 = 	mp.FluxRegion(center=mp.Vector3(Lr2,0), 	size=mp.Vector3(0,monitorheight))	# Reflected flux 2
	tran_fr = 	mp.FluxRegion(center=mp.Vector3(Lt,0), 		size=mp.Vector3(0,monitorheight))	# Transmitted flux
	su_fr = 	mp.FluxRegion(center=mp.Vector3(0, monitorheight/2),	size=mp.Vector3(a,0))	# Flux loss above the waveguide
	sd_fr = 	mp.FluxRegion(center=mp.Vector3(0,-monitorheight/2), 	size=mp.Vector3(a,0))	# Flux loss below the waveguide

	if resolution <= 50:
		epsform = "eps-000000.00"
	else:
		epsform = "eps-000000000"

	# T = 40;

	sim.use_output_directory(outputdir)
	# sim.use_output_directory()
	sim.run(mp.at_beginning(mp.output_epsilon),
			mp.at_every(1, mp.output_png(mp.Ez, "-RZc bluered -A " + outputdir + "/grating_validation-" + epsform + ".h5 -a gray:.2")),
			# mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
			until = T)

	nfreq = 1

	refl1 = sim.add_flux(fcen, df, nfreq, refl_fr1)
	refl2 = sim.add_flux(fcen, df, nfreq, refl_fr2)
	tran = 	sim.add_flux(fcen, df, nfreq, tran_fr)
	su = 	sim.add_flux(fcen, df, nfreq, su_fr)
	sd = 	sim.add_flux(fcen, df, nfreq, sd_fr)

	geometry_lattice = mp.Lattice(size=mp.Vector3(0, monitorheight, 0))

	# bandNum = 4

	# geometry = [mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
	# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
	# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material),
	# 			mp.Block(mp.Vector3(mp.inf, h, mp.inf), center=mp.Vector3(0,0,0), material=core_material)]

	geometrympb = [	mp.Block(mp.Vector3(mp.inf, monitorheight, 	mp.inf), center = mp.Vector3(0,0), 				material = default_material),
					mp.Block(mp.Vector3(mp.inf, hu+h/2, 		mp.inf), center = mp.Vector3(0,  (hu+h/2)/2), 	material = upper_material),
					mp.Block(mp.Vector3(mp.inf, hl+h/2, 		mp.inf), center = mp.Vector3(0, -(hl+h/2)/2), 	material = lower_material),
					mp.Block(mp.Vector3(mp.inf, h, 				mp.inf), center = mp.Vector3(0, 0), 			material = core_material)]

	ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
						resolution=resolution,
						# ensure_periodicity=False,
						geometry=geometrympb)

	os.system("convert " + outputdir + "/grating_validation-ez-*.png gratingFull-" + name + ".gif");
	os.system("rm " + outputdir + "/grating_validation-ez-*.png");

	# nearfield = sim.add_near2far(fcen, 0, 1,
	# 	mp.Near2FarRegion(mp.Vector3(0,  0.5*sxy), size=mp.Vector3(sxy)),
	# 	mp.Near2FarRegion(mp.Vector3(0, -0.5*sxy), size=mp.Vector3(sxy), weight=-1.0),
	# 	mp.Near2FarRegion(mp.Vector3( 0.5*sxy), size=mp.Vector3(0,sxy)),
	# 	mp.Near2FarRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(0,sxy), weight=-1.0))

	nearfield = sim.add_near2far(fcen, 0, 1,
		mp.Near2FarRegion(mp.Vector3(0,  	 	0.5*sxy), size=mp.Vector3(sxy)))#,
			# mp.Near2FarRegion(mp.Vector3(-0.5*sxy, 	0.3*sxy), size=mp.Vector3(0, 0.4*sxy), weight=-1.0),
		# mp.Near2FarRegion(mp.Vector3(0.5*sxy, 	0.3*sxy), size=mp.Vector3(0, 0.4*sxy)))

	sim.run(mp.at_every(wavelength/20 , mp.output_png(mp.Ez, "-RZc bluered -A " + outputdir + "/grating_validation-" + epsform + ".h5 -a gray:.2")), until=19*wavelength/20)
	sim.run(until = 19*wavelength)

	band = 1

	k = ms.find_k(
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

	print(k)

	neff = k[0]/fcen;

	print(neff)

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

	refl_val = sim.get_array(center=mp.Vector3(Lr1,0), size=mp.Vector3(0,monitorheight-2/resolution), component=mp.Ez, cmplx=True)
	tran_val = sim.get_array(center=mp.Vector3(Lt, 0), size=mp.Vector3(0,monitorheight-2/resolution), component=mp.Ez, cmplx=True)

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

	refl1_flux = 	mp.get_fluxes(refl1)
	refl2_flux = 	mp.get_fluxes(refl2)
	tran_flux = 	mp.get_fluxes(tran)
	su_flux = 		mp.get_fluxes(su)
	sd_flux = 		mp.get_fluxes(sd)

	index = 0

	R = 	-refl1_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
	T = 	 tran_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
	S = 	(refl2_flux[index] - tran_flux[index]) 	/ (refl2_flux[index] - refl1_flux[index])
	Su = 	 su_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])
	Sd = 	-sd_flux[index] 						/ (refl2_flux[index] - refl1_flux[index])

	print("Reflected Power: ", -R)
	print("Transmitted Power: ", -T)
	print("Scattered Power: ", +S)
	print("UpScattered Power: ", -Su)
	print("DownScattered Power: ", -Sd)

	print("Total Power: ", -R-T-S)

	pt = mp.Vector3(0,0)
	# sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

	os.system("convert " + outputdir + "/grating_validation-ez-*.png grating-" + name + ".gif");

	r = 800
	npts = 1000;

	Su = 0

	for n in range(npts):
		# print(n);
		ff = sim.get_farfield(nearfield, mp.Vector3(r*cos(2*pi*(n/npts)), r*sin(2*pi*(n/npts))))

		Ex = ff[0]
		Ey = ff[1]
		Ez = ff[2]

		Hx = ff[3]
		Hy = ff[4]
		Hz = ff[5]

		Py=((Ez * Hx)-(Ex * Hz)).real;
		Pz=((Ex * Hy)-(Ey * Hx)).real;
		Pr=sqrt((Py ** 2)+(Pz ** 2));

		Su += Pr

		#print("Pr: ", Pr)

		f1.write("{}, {}, ".format(n,2*pi*n/npts))
		f1.write(", ".join([str(f).strip('()').replace('j', 'i') for f in ff]))
		f1.write("\n")
		# print(n)

	print("Power Scattered Upward (Farfield): ", Su)
	print("Power Scattered Upward (Farfield, Normalized)", Su / (-R-T+S))
	#print("Power Reflected: ", fund_refl_power.real)
	#3print("Power Transmitted: ", fund_tran_power.real)

	f1.close();

print(sys.argv)
print(len(sys.argv))

if len(sys.argv) > 2:
	print(int(sys.argv[1]))
	print(int(sys.argv[2]))
	validation(int(sys.argv[1]), int(sys.argv[2]));
elif len(sys.argv) > 1:
	validation(int(sys.argv[1]), 0);
else:
	validation(0, 0);
