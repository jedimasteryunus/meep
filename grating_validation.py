import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt, cos, sin

os.system("rm -r grating_validation-out/")

f1 = open('grating_validation.out', 'w')

# print(f1)

def validation(w):

	print("#----------------------------------------")
	print("NOTCH WIDTH: %s nanometers" % (w * 1000))
	print("#----------------------------------------")

	# w is the width of the notch in the waveguide

	#Waveguide Math
	# a = 12 # length of cell
	# h = 0.2 # height of Waveguide
	# e = 0.4 # etch fraction [0 -> 1]
	# H = 3*a
	# # H = 6 #height of cell

	# h2 = h * (1 - e) #height of the etched region
	# d = 1 #spacing of the waveguides
	# posd = d / 2 #pos half d
	# negd = -d / 2 #neg half d
	#
	# ang = 80 #angle of the sidewall
	# angrad = ang * pi / 180
	# bottomoffset = h / tan(angrad)
	# c = bottomoffset / 2
	# r = sqrt(bottomoffset ** 2 + h ** 2)
	#
	# wavelength = 0.6372
	# # wavelength = 0.4203
	# fcen = 1 / wavelength
	# df = 0.05
	#
	# n_e = 2.2012 # refractive index inside waveguide (at wavelength 0.6372)
	# n_c = 1.4569 # refractive index outside waveguide (at wavelength 0.6372)
	#
	# #Waveguide Geometry
	# cell = mp.Vector3(H, H)
	#
	# default_material = mp.Medium(epsilon = n_c ** 2)
	# core_material = 	mp.Medium(epsilon = n_e ** 2)


	a = 16 						# length of cell
	h = 0.4 					# height of Waveguide
	dpml = 1;
	# monitorheight = .8
	monitorheight = .6
	H = monitorheight + 2*dpml 	#height of cell
	resolution = 50
	H = a
	sxy = H - 2*dpml
	# sxy2 = H - 3*dpml
	sxy2 = H

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
		hu = 6
		hl = 6

		ang = 80

	# f1 = open('notch-' + case + '.out', 'w')
	# f2 = open('notch-' + case + '.txt', 'w')

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


	# angrad = ang * pi / 180
	# bottomoffset = e * h / tan(angrad)
	#
	# vertices = [mp.Vector3( w/2 + bottomoffset, (.5 + e)*h),
	# 			mp.Vector3(-w/2 - bottomoffset, (.5 + e)*h),
	# 			mp.Vector3(-w/2 + bottomoffset, (.5 - e)*h),
	# 			mp.Vector3( w/2 - bottomoffset, (.5 - e)*h)]
	#
	# if bottomoffset > w/2:
	# 	ratio = (w/2)/bottomoffset
	# 	vertices = [mp.Vector3( w/2, h/2),
	# 				mp.Vector3(-w/2, h/2),
	# 				mp.Vector3(0, (.5 - e*ratio)*h)]
	#
	# print(vertices)

	#Waveguide Geometry
	cell = mp.Vector3(a, H)

	geometry = [mp.Block(cell, 						center = mp.Vector3(0,0), 				material = default_material),
				mp.Block(mp.Vector3(sxy2, hu+h/2), 	center = mp.Vector3(0,  (hu+h/2)/2), 	material = upper_material),
				mp.Block(mp.Vector3(sxy2, hl+h/2), 	center = mp.Vector3(0, -(hl+h/2)/2), 	material = lower_material),
				mp.Block(mp.Vector3(sxy2, h), 		center = mp.Vector3(0,0), 				material = core_material)]

	#old_grate_positions = [50., 740., 1470., 2170., 2850., 3390., 3780., 4170., 4560., 4950.]
	# old_grate_positions = [50.,   810.,  1560.,  1920.,  2630.,  3370.,  4110.,  4690.,  5020.,  5600.]
	# old_grate_positions = range(150,12*300,300)
	# old_grate_positions = range(150,12*300,350)
	# old_grate_positions = [640., 500., 680., 670., 490., 650., 620., 510., 150., 600.]
	# old_grate_positions = [650., 695., 675., 630., 405., 710., 350., 625., 490., 600.]
	# old_grate_positions = [675., 675., 665., 680., 325., 625., 160., 805., 310., 600.]
	# old_grate_positions = [335.0, 375.0, 330.0, 130.0, 325.0, 355.0, 350.0, 340.0, 330.0, 125.0, 325.0, 360.0, 325.0, 390.0, 320.0, 485.0, 385.0, 250.0, 110.0, 300.0]
	# old_grate_positions = [490.0, 315.0, 315.0, 340.0, 305.0, 460.0, 385.0, 290.0, 335.0, 480.0, 290.0, 350.0, 475.0, 480.0, 170.0, 305.0, 330.0, 155.0, 325.0, 300.0]
	# old_grate_positions = [290, 495, 285, 495, 295, 485, 300, 485, 285, 500, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# old_grate_positions = [290, 495, 285, 495, 295, 485, 300, 485, 285, 500, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# old_grate_positions = [395, 395, 395, 395, 395, 395, 395, 395, 395, 395, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# old_grate_positions = [250, 220, 255, 210, 255, 190, 270, 200, 275, 185, 275, 175, 285, 175, 280, 205, 270, 190, 280, 300]
	# old_grate_positions = [260, 220, 240, 220, 270, 205, 275, 190, 245, 200, 290, 170, 290, 185, 280, 190, 255, 185, 300, 300]

	# old_grate_positions = [885, 900, 920, 945, 965, 995, 335, 590, 585, 100]
	# old_grate_positions = [915, 915, 990, 940, 1000, 815, 605, 405, 850, 100]
	# old_grate_positions = [165, 355, 165, 370, 155, 170, 380, 155, 170, 390, 180, 140, 185, 400, 400, 405, 195, 405, 190, 100]
	# old_grate_positions = [355, 555, 365, 160, 365, 175, 155, 165, 170, 165, 390, 160, 165, 195, 380, 415, 395, 405, 485, 100];
	# old_grate_positions = [160, 155, 160, 165, 160, 175, 150, 170, 395, 400, 400, 410, 400, 380, 380];
	# old_grate_positions = [115, 125, 130, 130, 160, 140, 185, 410, 495, 400];
	# old_grate_positions = [115, 125, 130, 130, 160, 140, 185];
	old_grate_positions = [110, 340, 130, 135, 160, 165, 375];
	# geometry = [mp.Block(cell,
	# 				center = mp.Vector3(0,0),
	# 				material = default_material),
	# 			mp.Block(mp.Vector3(a -2.1*dpml, h),
	# 				center = mp.Vector3(0,0),
	# 				material = core_material)]#,
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[0], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[1], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[2], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[3], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[4], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[5], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[6], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[7], h2/2),
				#     material = default_material),
				# mp.Block(mp.Vector3(w, e*h),
				#     center = mp.Vector3(grate_positions[8], h2/2),
				#     material = default_material),
				# ]

	x = 0;

	for elt in old_grate_positions:
		x -= w + elt / 1000

	x /= 2;

	# grate_positions = []
	w = .050
	# x = 3*dpml - a/2;
	# h2 =
	for elt in old_grate_positions:
		geometry.append(mp.Block(mp.Vector3(w, e*h),
			center = mp.Vector3(x, h * (1 - e)/2),
			material = default_material))
		x += w + elt / 1000
		# grate_positions.append(x)

		# geometry.append(mp.Block(mp.Vector3(w, e*h),
		# 	center = mp.Vector3(x, h * (1 - e)/2),
		# 	material = default_material))

	sources = [mp.EigenModeSource(mp.ContinuousSource(frequency = fcen),
								  size = mp.Vector3(0,H),
								  center = mp.Vector3(2*dpml-a/2, 0))]

	# pml_layers = [mp.PML(0.2)]
	pml_layers = [mp.Absorber(thickness=dpml)]

	# resolution = 50

	sim = mp.Simulation(cell_size = cell,
						boundary_layers = pml_layers,
						geometry = geometry,
						sources = sources,
						resolution = resolution)


	#------------------------------------------------------
	#FOR GENERATING THE ELECTRIC FIELD GIF
	#Note: After running this program, write the following commands in Terminal:
		# $ source deactivate mp
		# $ cd grating_validation-out/
		# $ python ../GratingValidationIP.py

	sim.use_output_directory()
	sim.run(mp.at_beginning(mp.output_epsilon),
			mp.at_every(1, mp.output_png(mp.Ez, "-RZc bluered -A grating_validation-out/grating_validation-eps-000000.00.h5 -a gray:.2")),
			# mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z)),
			until = 40)

	os.system("convert grating_validation-out/grating_validation-ez-*.png gratingFull.gif");
	os.system("rm grating_validation-out/grating_validation-ez-*.png");

	# nearfield = sim.add_near2far(fcen, 0, 1,
	# 	mp.Near2FarRegion(mp.Vector3(0,  0.5*sxy), size=mp.Vector3(sxy)),
	# 	mp.Near2FarRegion(mp.Vector3(0, -0.5*sxy), size=mp.Vector3(sxy), weight=-1.0),
	# 	mp.Near2FarRegion(mp.Vector3( 0.5*sxy), size=mp.Vector3(0,sxy)),
	# 	mp.Near2FarRegion(mp.Vector3(-0.5*sxy), size=mp.Vector3(0,sxy), weight=-1.0))

	nearfield = sim.add_near2far(fcen, 0, 1, mp.Near2FarRegion(mp.Vector3(0,  0.5*sxy), size=mp.Vector3(sxy)))

	sim.run(mp.at_every(wavelength/20 , mp.output_png(mp.Ez, "-RZc bluered -A grating_validation-out/grating_validation-eps-000000.00.h5 -a gray:.2")), until=19*wavelength/20)
	sim.run(until = 19*wavelength)

	pt = mp.Vector3(0,0)
	# sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

	os.system("convert grating_validation-out/grating_validation-ez-*.png grating.gif");

	r = 800
	npts = 1000;

	for n in range(npts):
		print(n);
		ff = sim.get_farfield(nearfield, mp.Vector3(r*cos(2*pi*(n/npts)), r*sin(2*pi*(n/npts))))
		f1.write("{}, {}, ".format(n,2*pi*n/npts))
		f1.write(", ".join([str(f).strip('()').replace('j', 'i') for f in ff]))
		f1.write("\n")
		print(n)
	f1.close();

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

	pt = mp.Vector3(9.75,0)

	sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

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
	'''
# validation(0.05)
validation(0.1)
