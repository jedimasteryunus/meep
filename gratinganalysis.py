from math import pi, sqrt
from random import random, randint
import time
import ast
import numpy as np
import matplotlib.pyplot as plt

wavelength = 637
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
print(data)

def transmission(width): 
	index = widthToIndex(width)
	t = sqrt(data[2][index] / 100)
	return t

def reflection(width): 
	index = widthToIndex(width)
	r = sqrt(data[1][index] / 100)
	return r

def scatter(width): 
    index = widthToIndex(width)
    s = sqrt(data[-1][index] / 100)
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

    #for i in range(0, num_notches): 
        #print("Scatter:", scatter(widths[i]))
        #print("Amplitudes:", amplitudes[0,2*i] + amplitudes[1,2*i+1])

    for i in range(0, num_notches): 
        s[i] = scatter(widths[i])*(amplitudes[0,2*i] + amplitudes[1,2*i+1])

    final_gamma = 0
    for X in range(0, grating_length + 100, 100):
        gamma = 0
        for i in range(0, num_notches): 
            gamma += s[i] * E(lengths[i] - X, W)
        gamma = np.abs(gamma)**2
        final_gamma = max(gamma, final_gamma)
    return final_gamma

def neighbor(lengths):
    N = len(lengths)
    sensitivity = 20
    rand_index = randint(0, N - 1)
    rand_val = randint(5, 25) * sensitivity
    new_lengths = lengths.copy()
    new_lengths[rand_index] = rand_val
    return new_lengths

def acceptance_probability(c_old, c_new, T):
    frac = (c_old - c_new) / T
    prob = np.exp(frac)
    return prob

def anneal(widths, lengths, amplitudes, W, grating_length, num_notches):
    gamma = getOverlap(widths, lengths, amplitudes, W, grating_length, num_notches)

    T = 1.0
    T_min = 0.00001
    alpha = 0.9

    while T > T_min:
        i = 1

        while i <= 100:
            new_lengths = neighbor(lengths)
            gamma_new = getOverlap(widths, lengths, amplitudes, W, grating_length, num_notches)
            #print(gamma_new)
            ap = acceptance_probability(gamma, gamma_new, T)

            if ap > random():
                lengths = new_lengths
                gamma = gamma_new

            i += 1

        T = T*alpha

    return lengths, gamma

def main():
    N = 10
    NA = 200 #Numerical Aperture
    W = wavelength / (pi * NA) #mode field diameter
    grating_length = 500

    widths =    100 * np.ones(N)
    lengths =   100 * np.ones(N)

    amplitudes = getAmplitudes(widths, lengths)

    #gamma = getOverlap(widths, lengths, amplitudes, W, grating_length, N) #Test for getOverlap function

    start = time.time()

    #print(gamma) #Test for getOverlap function

    print(anneal(widths, lengths, amplitudes, W, grating_length, N))

    end = time.time()

    print("Run Time: ", end - start)

'''
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
'''
main() 
