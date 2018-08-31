import os
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from math import pi, tan, sqrt

#os.system("rm -r notch-out/")

f1 = open('grating_validation.out', 'w')

def validation(w):

    print("#----------------------------------------")
    print("NOTCH WIDTH: %s nanometers" % (w * 1000))
    print("#----------------------------------------")

    # w is the width of the notch in the waveguide

    #Waveguide Math
    a = 8 # length of cell
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

    #old_grate_positions = [50., 740., 1470., 2170., 2850., 3390., 3780., 4170., 4560., 4950.]
    old_grate_positions = [50.,   800.,  1570.,  1860.,  2580.,  3310.,  4060.,  4640.,  5370.,  5950.]

    grate_positions = []
    for elt in old_grate_positions:
        grate_positions.append(elt / 1000 - 2)

    geometry = [mp.Block(mp.Vector3(a + 2 * h, h),
                center = mp.Vector3(0,0),
                material = mp.Medium(epsilon = n_e ** 2)),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[0], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[1], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[2], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[3], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[4], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[5], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[6], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[7], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[8], h2 / 2),
                material = mp.Medium()),
                mp.Block(mp.Vector3(w, e * h),
                center = mp.Vector3(grate_positions[9], h2 / 2),
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

validation(0.1)
