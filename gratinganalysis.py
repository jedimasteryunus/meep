from math import pi, sqrt
import ast
import numpy as np
import matplotlib.pyplot as plt

wavelength = 0.637

def widthToIndex(width): 
	return width / 20 - 2

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
	t = transmission(width)
	r = reflection(width)

	return np.matrix([  [ 1/t,   r/t             ],
                            [ r/t,   (r*r + t*t)/t   ]  ], dtype=complex)

def propagate(length): 
	return np.matrix([  [ exp(-pi*2j*neff*length/wavelength),    0                               ],
                            [ 0,                                 exp(pi*2j*neff*length/wavelength)  ]  ], dtype=complex)


def getAmplitudes(widths, lengths):
    N = widths.shape[0]                                 # N grates.

    a = np.zeros((2, 2*N))                              # N times [[ a  b  ]; [ b' a' ]].

    a[0, 2*N-1] = 1                                     # Set the initial vector (b = 1, a' = 0).
    a[:, 2*N-2] = scatter(widths[N-1]) * a[:, 2*N-1]    # And find a, b' for the first grate.

    for ii in range(N-1):                               # Now do this for the rest of the grates.
        a[:, 2*ii+1] = propogate(lengths[ii]) * a[:, 2*ii+2] 
        a[:, 2*ii] = scatter(   widths[ii]) * a[:, 2*ii+1] 

    return a 


def getOverlap(widths, lengths, amplitudes, NA):
    gamma = 0  # Finish me...

    return gamma 

def main():
    N = 10

    widths =    .1 * np.ones(N)
    lengths =   .3 * np.ones(N)

    gammaprev = 0

    optimize = 1

    while optimize:
        a =  getAmplitudes(widths, lengths)
        gamma = getOverlap(widths, lengths, a, NA)

        if gamma < gammaprev:               # If our figure of merit got worse...
            widths =        widthsprev      # ...then revert the widths and lengths.
            lengths =       lengthsprev 
        else:                               # If our figure of merit got better...
            widthsprev =    widths          # ...then proceed.
            lengthsprev =   lengths 

        widths =    evolve(widths)      # Do something to change the widths (perhaps annealing process)
        lengths =   evolve(lengths)

main() 
