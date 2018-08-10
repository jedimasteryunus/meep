from math import pi, sqrt, exp
import ast
import numpy as np
import matplotlib.pyplot as plt

wavelength = 0.637
neff = 2.2012 # refractive index inside waveguide (at wavelength 0.6372)

def widthToIndex(width): 
	return int(width / 20 - 2)

def importData(fname): 
	f = open(fname, "r")
	data = []
	for line in f.readlines():
    		line = ast.literal_eval(line)
    		data.append(line)
	f.close()
	return data

fname = "notch.txt"
data = importData(fname)

def transmission(width): 
	index = widthToIndex(width)
	t = sqrt(data[2][index])
	return t

def reflection(width): 
	index = widthToIndex(width)
	r = sqrt(data[1][index])
	return r

def scatter(width): 
    index = widthToIndex(width)
    s = sqrt(data[-1][index])
    return s

def scatter_matrix(width): 
	t = transmission(width)
	r = reflection(width)

	return np.matrix([  [ 1/t,   r/t             ],
                        [ r/t,   (r*r + t*t)/t   ]  ], 
                        dtype=complex)

def propagate_matrix(length): 
    return np.matrix([  [ np.exp(-pi*2j*neff*length/wavelength),              0                     ],
                        [ 0,                                 np.exp(pi*2j*neff*length/wavelength)   ]  ], 
                        dtype=complex)


def getAmplitudes(widths, lengths):
    N = widths.shape[0]                                 # N grates.

    a = np.zeros((2, 2*N), dtype=complex)                              # N times [[ a  b  ]; [ b' a' ]].

    a[0, 2*N-1] = 1                                     # Set the initial vector (b = 1, a' = 0).
    a[:, 2*N-2] = a[:, 2*N-1] * scatter_matrix(widths[N-1])   # And find a, b' for the first grate.

    for ii in range(N-1):                               # Now do this for the rest of the grates.
        a[:, 2*ii+1] = a[:, 2*ii+2] * propagate_matrix(lengths[ii])
        a[:, 2*ii] = a[:, 2*ii+1]  * scatter_matrix(widths[ii])

    return a 


def E(x,W):
    return np.exp(-x**2 / W**2)

def getOverlap(widths, lengths, amplitudes, W, grating_length, num_notches):
    
    s = np.zeros(num_notches, dtype = complex)

    for i in range(0, num_notches): 
        print("Scatter:", scatter(widths[i]))
        print("Amplitudes:", amplitudes[2*i,0] + amplitudes[2*i+1,1])

    for i in range(0, num_notches): 
        s[i] = scatter(widths[i])*(amplitudes[2*i,0] + amplitudes[2*i+1,1])

    final_gamma = 0
    for X in range(0, grating_length + 100, 100):
        gamma = 0
        for i in range(0, num_notches): 
            gamma += s[i] * E(lengths[i] - X, W)
        gamma = np.abs(gamma)**2
        final_gamma = max(gamma, final_gamma)
    return final_gamma


def main():
    N = 10
    NA = 0.2 #Numerical Aperture
    W = wavelength / (pi * NA) #mode field diameter
    grating_length = 500


    widths =    .1 * np.ones(N)
    lengths =   .3 * np.ones(N)

    gammaprev = 0

    optimize = 1

    while optimize:
        a =  getAmplitudes(widths, lengths)
        gamma = getOverlap(widths, lengths, a, W, grating_length, N)

        if gamma < gammaprev:               # If our figure of merit got worse...
            widths =        widthsprev      # ...then revert the widths and lengths.
            lengths =       lengthsprev 
        else:                               # If our figure of merit got better...
            widthsprev =    widths          # ...then proceed.
            lengthsprev =   lengths 

        widths =    evolve(widths)      # Do something to change the widths (perhaps annealing process)
        lengths =   evolve(lengths)

main() 
