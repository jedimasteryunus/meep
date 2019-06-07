import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt, cos, sin
import sys
import operator
import time

def validation(dw, dl, farfield_bool):

    name = "dw=%snm-dl=%snm" % (int(dw), int(dl))
    print(name)

    outputdir = "grating_validation-out-" + name
    print(outputdir)

    #os.system("rm -r " + outputdir + "/")

    T = 40
    a = 24							# length of cell
    h = 0.4 					# height of Waveguide
    monitorheight = 3 * h

    ALN420_bool = False

    if farfield_bool:
        f1 = open('grating_validation-Farfield.out', 'w')
        dpml = 1
        H = a  # height of cell
    else:
        f1 = open('grating_validation-Nearfield.out', 'w')
        dpml = 3 * h
        H = monitorheight + 2 * dpml  # height of cell

    resolution = 50

    sxy = a - 2*dpml
    sxy2 = a

    if ALN420_bool:
        case = "AlN420"
        wavelength = 0.4203

        n_e = 2.2415
        n_o = 2.1888
        n_x = n_o
        n_y = n_e
        n_z = n_o

        n_u = 1.4681
        n_l = 1.7800

        e = .5

        h = .2
        hu = .2     #Height of upper cladding
        hl = 100	#Height of lower cladding

        ang = 90
    else:
        case = "LN637"
        wavelength = 0.6372

        n_e = 2.2012
        n_o = 2.3048
        n_x = n_o
        n_y = n_o
        n_z = n_e

        n_u = 1.4569
        n_l = 1.4569

        e = .4

        h = .2
        hu = 6		#Height of upper cladding
        hl = 6		#Height of lower cladding

        ang = 80

    hu = min(hu, (sxy2-h)/2)    #Height of upper cladding
    hl = min(hl, (sxy2-h)/2)	#Height of lower cladding

    n_m = max(n_x, n_y, n_z)
    n_c = max(n_l, n_u)

    fcen = 1 / wavelength
    df = 0.05

    default_material = mp.Medium(epsilon=1)
    upper_material = mp.Medium(epsilon=n_u ** 2)
    core_material = mp.Medium(epsilon_diag=mp.Vector3(n_x ** 2, n_y ** 2, n_z ** 2))
    lower_material = mp.Medium(epsilon=n_l ** 2)

    # Waveguide Geometry
    cell = mp.Vector3(a, H)

    geometry = [mp.Block(cell, 						center=mp.Vector3(0, 0), 				material=default_material),
                mp.Block(mp.Vector3(sxy2, hu+h/2), 	center=mp.Vector3(0,
                                                                      (hu+h/2)/2), 	material=upper_material),
                mp.Block(mp.Vector3(sxy2, hl+h/2), 	center=mp.Vector3(0, -
                                                                      (hl+h/2)/2), 	material=lower_material),
                mp.Block(mp.Vector3(sxy2, h), 		center=mp.Vector3(0, 0), 				material=core_material)]

    input_lengths = [31, 20, 36, 32, 15, 4, 31, 16, 32, 16, 32, 16, 16, 16, 32, 16, 32, 32, 32, 9]
    num_notches = len(input_lengths)
    #lengths = input_lengths
    lengths = [10 * length for length in input_lengths]
    widths = (num_notches + 1) * [100]

    try:
        assert(len(lengths) + 1 == len(widths))
    except:
        print("WARNING: There is a mismatch between the number of notches and the number of inter-notch lengths.")
        quit()

    for i in range(0, len(widths)):
        widths[i] += dw

    for i in range(0, len(lengths)):
        lengths[i] += dl

    x = -widths[0]/2000.

    for i in range(0, len(lengths)):
        x -= (widths[i+1] + lengths[i]) / 1000.

    x /= 2

    x_list = []
    for i in range(0, len(lengths)):
        x_list.append(x)
        geometry.append(mp.Block(mp.Vector3(widths[i]/1000., e*h),
                                 center=mp.Vector3(x, h * (1 - e)/2),
                                 material=default_material))
        x += (widths[i+1] + lengths[i]) / 1000.

    geometry.append(mp.Block(mp.Vector3(widths[len(widths)-1]/1000., e*h),
                             center=mp.Vector3(x, h * (1 - e)/2),
                             material=default_material))

    Ls = -2*x  # Position of source
    Lr1 = -2*x - 1.25  # Position of reflection monitor 1
    Lr2 = -2*x + 1.25  # Position of reflection monitor 2
    Lt = -Lr2  # Position of transmission monitor
    scatter_monitor_size = 2 * Lt

    assert(Lr1 > - a / 2 + dpml)  # Make sure that nothing we care about is sitting inside the PML

    sources = [mp.EigenModeSource(mp.ContinuousSource(frequency=fcen),
                                  size=mp.Vector3(0, H),
                                  center=mp.Vector3(Ls, 0))]

    pml_layers = [mp.Absorber(thickness=dpml)]

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        resolution=resolution)

    sim.reset_meep()

    nfreq = 1

    refl_fr1 = mp.FluxRegion(center=mp.Vector3(Lr1, 0),  size=mp.Vector3(0, monitorheight))
    refl_fr2 = mp.FluxRegion(center=mp.Vector3(Lr2, 0),  size=mp.Vector3(0, monitorheight))
    tran_fr = mp.FluxRegion(center=mp.Vector3(Lt, 0),      size=mp.Vector3(0, monitorheight))
    su_fr = mp.FluxRegion(center=mp.Vector3(0, monitorheight/2),
                          size=mp.Vector3(scatter_monitor_size, 0))
    sd_fr = mp.FluxRegion(center=mp.Vector3(0, -monitorheight/2),
                          size=mp.Vector3(scatter_monitor_size, 0))

    if resolution <= 50:
        epsform = "eps-000000.00"
    else:
        epsform = "eps-000000000"

    sim.use_output_directory(outputdir)

    # Simulation run for TRANSIENT state
    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.at_every(1, mp.output_png(mp.Ez, "-RZc bluered -A " + outputdir +
                                         "/grating_validation-" + epsform + ".h5 -a gray:.2")),
            until=T)

    if farfield_bool:
        os.system("convert " + outputdir + "/grating_validation-ez-*.png " +
                  outputdir + "/N=" + str(num_notches) + "-FarfieldTransient.gif")
    else:
        os.system("convert " + outputdir + "/grating_validation-ez-*.png " +
                  outputdir + "/N=" + str(num_notches) + "-NearfieldTransient.gif")

    os.system("rm " + outputdir + "/grating_validation-ez-*.png")

    # Delete after printing output
    refl1 = sim.add_flux(fcen, df, nfreq, refl_fr1)
    refl2 = sim.add_flux(fcen, df, nfreq, refl_fr2)
    tran = sim.add_flux(fcen, df, nfreq, tran_fr)
    su = sim.add_flux(fcen, df, nfreq, su_fr)
    sd = sim.add_flux(fcen, df, nfreq, sd_fr)

    nearfield = sim.add_near2far(fcen, 0, 1,
                                 mp.Near2FarRegion(mp.Vector3(0, 0.5*sxy), size=mp.Vector3(sxy)))  #

    old_refl1_flux = mp.get_fluxes(refl1)
    old_refl2_flux = mp.get_fluxes(refl2)
    old_tran_flux = mp.get_fluxes(tran)
    old_su_flux = mp.get_fluxes(su)
    old_sd_flux = mp.get_fluxes(sd)

    # Simulation run for STEADY state
    sim.run(mp.at_every(1, mp.output_png(mp.Ez, "-RZc bluered -A " + outputdir +
                                         "/grating_validation-" + epsform + ".h5 -a gray:.2")), until=19*wavelength/20)
    sim.run(until=19*wavelength)

    # Delete after printing output
    refl1_flux = mp.get_fluxes(refl1)
    refl2_flux = mp.get_fluxes(refl2)
    tran_flux = mp.get_fluxes(tran)
    su_flux = mp.get_fluxes(su)
    sd_flux = mp.get_fluxes(sd)

    incident_flux = refl2_flux[0] - refl1_flux[0]

    pt = mp.Vector3(0, 0)

    if farfield_bool:
        os.system("convert " + outputdir + "/grating_validation-ez-*.png " +
                  outputdir + "/N=" + str(num_notches) + "-FarfieldSteadyState.gif")
    else:
        os.system("convert " + outputdir + "/grating_validation-ez-*.png " +
                  outputdir + "/N=" + str(num_notches) + "-NearfieldSteadyState.gif")

    os.system("rm " + outputdir + "/grating_validation-ez-*.png")

    r = 800
    npts = 1000

    if farfield_bool:
        Su = 0

        for n in range(npts):

            ff = sim.get_farfield(nearfield, mp.Vector3(r*cos(2*pi*(n/npts)), r*sin(2*pi*(n/npts))))

            Ex = ff[0]
            Ey = ff[1]
            Ez = ff[2]

            Hx = ff[3]
            Hy = ff[4]
            Hz = ff[5]

            Py = ((Ez * Hx)-(Ex * Hz)).real
            Pz = ((Ex * Hy)-(Ey * Hx)).real
            Pr = sqrt((Py ** 2)+(Pz ** 2))

            Su += Pr

            f1.write("{}, {}, ".format(n, 2*pi*n/npts))
            f1.write(", ".join([str(f).strip('()').replace('j', 'i') for f in ff]))
            f1.write("\n")

    sim.reset_meep()

    print("Number of Notches: ", num_notches)

    print("X List: ", x_list)
    print("Position Comparison:", "Lr1 = %s," % (Lr1), "-2x = %s," %
          (-2*x), "Lr2 = %s," % (Lr2), "Lt = %s," % (Lt))

    left_notch_boundary = x_list[0] - widths[0] / 2000.
    right_notch_boundary = x_list[-1] + widths[-1] / 2000.
    print("Left Notch Boundary: ", left_notch_boundary)
    print("Right Notch Boundary: ", right_notch_boundary)

    assert(Lr2 < left_notch_boundary)
    assert(Lt > right_notch_boundary)

    print("Refl1 Flux: ", refl1_flux[0])
    print("Old Refl1 Flux: ", old_refl1_flux[0])
    print("Refl2 Flux: ", refl2_flux[0])
    print("Old Refl2 Flux: ", old_refl2_flux[0])
    print("Tran Flux: ", tran_flux[0])
    print("Old Tran Flux: ", old_tran_flux[0])
    print("Upward-Scattered Flux: ", su_flux[0])
    print("Old Upward-Scattered Flux: ", old_su_flux[0])
    print("Downward-Scattered Flux: ", sd_flux[0])
    print("Old Downward-Scattered Flux: ", old_sd_flux[0])

    print("Gauss' Law for Magnetism Sanity Check: ",
          su_flux[0] - sd_flux[0] - refl2_flux[0] + tran_flux[0])

    print("Incident Flux: ", incident_flux)

    print("Power Scattered Upward: ", su_flux[0])
    print("Power Scattered Upward (Normalized): ", su_flux[0] / incident_flux)

    print("Grating Efficiency: %s Percent" % (100 * su_flux[0] / incident_flux))

    '''
	if farfield_bool:
		print("Power Scattered Upward (Farfield): ", Su)
		print("Power Scattered Upward (Farfield, Normalized): ", Su / incident_flux)

		print("Grating Efficiency (Farfield): %s Percent" % (100 * Su / incident_flux))

		f1.close()

		result = 100 * Su / incident_flux

		return result
	'''

    f1.close()

    result = 100 * su_flux[0] / incident_flux

    '''
	# Memory cleanup for future validation runs
	del sim
	del sources
	del input_lengths
	del num_notches
	del widths
	del lengths
	del refl_fr1
	del refl1
	del refl1_flux
	del refl_fr2
	del refl2
	del refl2_flux
	del tran_fr
	del tran
	del tran_flux
	del su_fr
	del su
	del su_flux
	del sd_fr
	del sd
	del sd_flux
	del incident_flux
	'''

    return result


def sweep(dw_lower_bound, dw_upper_bound, dl_lower_bound, dl_upper_bound, farfield_bool):
    start = time.time()
    efficiency_dict = dict()
    for dw in range(dw_lower_bound, dw_upper_bound, 10):
        for dl in range(dl_lower_bound, dl_upper_bound, 10):
            efficiency_dict[(dw, dl)] = validation(dw, dl, farfield_bool)
    sorted_efficiency_dict = sorted(efficiency_dict.items(),
                                    key=operator.itemgetter(1), reverse=True)
    print("Sorted Efficiency Dictionary: ", sorted_efficiency_dict)
    end = time.time()
    runtime = end - start
    print("Total Runtime: %s seconds" % (runtime))


dw_lower_bound = -10  # Note: This bound is inclusive
dw_upper_bound = 30  # Note: This bound is exclusive
dl_lower_bound = dw_lower_bound  # Note: This bound is inclusive
dl_upper_bound = dw_upper_bound  # Note: This bound is exclusive
# dl_lower_bound = 0 #Note: This bound is inclusive
# dl_upper_bound = 30 #Note: This bound is exclusive

try:
    assert(sys.argv[-1] == "nearfield" or sys.argv[-1] == "farfield")
except:
    print("WARNING: You did not specify whether to run a nearfield simulation or a farfield simulation.")
    quit()


if len(sys.argv) > 3:
    dw = int(sys.argv[1])
    dl = int(sys.argv[2])
    farfield = sys.argv[3]

    if farfield == "farfield":
        farfield_bool = True
    else:
        farfield_bool = False

    print(dw)
    print(dl)
    print(farfield_bool)

    validation(dw, dl, farfield_bool)

elif len(sys.argv) > 2:

    farfield = sys.argv[2]

    if farfield == "farfield":
        farfield_bool = True
    else:
        farfield_bool = False

    if sys.argv[1] == "sweep":
        sweep(dw_lower_bound, dw_upper_bound, dl_lower_bound, dl_upper_bound, farfield_bool)
    else:
        validation(int(sys.argv[1]), 0, farfield_bool)

else:

    farfield = sys.argv[1]

    if farfield == "farfield":
        farfield_bool = True
    else:
        farfield_bool = False

    validation(0, 0, farfield_bool)
