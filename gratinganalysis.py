from math import pi, sqrt
from random import random, randint, choice
import time
import ast
import numpy as np
import matplotlib.pyplot as plt

wavelength = 637
neff = 2.2012 # refractive index inside waveguide (at wavelength 0.6372)

debug = 0

# DATA #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
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
    t = sqrt(data[2][index] / 100)
    return t

def reflection(width):
    index = widthToIndex(width)
    r = sqrt(data[1][index] / 100)
    return r

def scatter(width):
    index = widthToIndex(width)
    s = sqrt(data[-1][index])
    return s

# MATRICES #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
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

# AMPLITUDES #### #### #### #### #### #### #### #### #### ####
def getAmplitudes(widths, lengths):
    N = widths.shape[0]                                 # N grates.

    a = np.zeros((2, 2*N), dtype=complex)                              # N times [[ a  b  ]; [ b' a' ]].

    a[0, 2*N-1] = 1                                     # Set the initial vector (b = 1, a' = 0).
    a[:, 2*N-2] = a[:, 2*N-1] * scatter_matrix(widths[N-1])   # And find a, b' for the first grate.

    for ii in range(N-2, -1, -1):                               # Now do this for the rest of the grates.
        a[:, 2*ii+1] = a[:, 2*ii+2] * propagate_matrix(lengths[ii])
        a[:, 2*ii] =   a[:, 2*ii+1]  * scatter_matrix(widths[ii])

    return a


# OVERLAP #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
def E(x,W):
    return np.exp(-x**2 / W**2)

def getOverlap(widths, lengths, amplitudes, W, grating_length, num_notches):

    s = np.zeros(num_notches, dtype = complex)
    x = np.zeros(num_notches)

    #for i in range(0, num_notches):
        #print("Scatter:", scatter(widths[i]))
        #print("Amplitudes:", amplitudes[0,2*i] + amplitudes[1,2*i+1])

    currentx = 0;

    for i in range(0, num_notches):
        s[i] = scatter(widths[i])*(amplitudes[0,2*i] + amplitudes[1,2*i+1])
        x[i] = currentx + widths[i]/2;
        currentx += lengths[i] + widths[i];

    s /= amplitudes[0,0]

    if debug:
        S = np.zeros(num_notches)

        for i in range(0, num_notches):
            S[i] = 100 * np.abs(s[i])**2

        print(S)

    final_reflection =      np.abs(amplitudes[1,0]/amplitudes[0,0])**2
    final_transmission =    np.abs(1/amplitudes[0,0])**2

    final_gamma = 0
    final_X = 0;
    sensitivity = 20;
    Wround = int(W/100)*100;
    for X in range(Wround, int(np.max(x)) + sensitivity-Wround, sensitivity):
        gamma = 0
        integral = 0
        for i in range(0, num_notches):
            gamma       += E(x[i] - X, W) * s[i]
            integral    += E(x[i] - X, W) * E(x[i] - X, W)
        gamma = (np.abs(gamma/integral)**2)
        if gamma > final_gamma:
            final_gamma = gamma
            final_X = X;
        # print gamma
    # print final_X;

    print("Lengths:           ", lengths)
    print("Grate positions:   ", x)

    return [final_gamma, final_transmission, final_reflection, final_X]

# ANNEALING #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
def neighbor(lengths):
    N = len(lengths)
    sensitivity = 10
    rand_index = randint(0, N - 2)
    new_lengths = lengths.copy()

    lowerbound = 100
    upperbound = 700

    # rand_val = randint(5, 25) * sensitivity
    rand_val = choice(range(lowerbound, upperbound+sensitivity, sensitivity))
    new_lengths[rand_index] = rand_val

    # if random() > .5:
    #     new_lengths[rand_index] += sensitivity;
    # else:
    #     new_lengths[rand_index] -= sensitivity;
    #
    # if new_lengths[rand_index] < lowerbound:
    #     new_lengths[rand_index] += 2*sensitivity;
    #
    # if new_lengths[rand_index] > upperbound:
    #     new_lengths[rand_index] -= 2*sensitivity;

    return new_lengths

def acceptance_probability(c_old, c_new, T):
    return np.exp((c_new - c_old) / T)

def anneal(widths, lengths, W, grating_length, num_notches):
    gtrx = getOverlap(widths, lengths, getAmplitudes(widths, lengths), W, grating_length, num_notches)

    T = 0.01
    T_min = 0.001
    alpha = 0.99

    while T > T_min:
        i = 1

        while i <= 100:
            print("log(T): ", np.log10(T))

            new_lengths = neighbor(lengths)
            new_gtrx = getOverlap(widths, new_lengths, getAmplitudes(widths, new_lengths), W, grating_length, num_notches)

            ap = acceptance_probability(gtrx[0], new_gtrx[0], T)
            r = random();

            if ap > r:
                if gtrx[0] < new_gtrx[0]:
                    print(" +\t",)
                else:
                    print(" -\t",)

                lengths = new_lengths
                gtrx = new_gtrx
            else:
                print("  \t",)

            print(    "S = {:.2f}%,\tT = {:.2f}%,\tR = {:.2f}%,\tX = {};".format(gtrx[0]*100, gtrx[1]*100, gtrx[2]*100, gtrx[3]))
            print(    "S'= {:.2f}%,\tT'= {:.2f}%,\tR'= {:.2f}%,\tX'= {}.".format(new_gtrx[0]*100, new_gtrx[1]*100, new_gtrx[2]*100, new_gtrx[3]))

            i += 1

        T = T*alpha

    return lengths, gtrx

# MAIN #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# [590. 630. 600. 580. 290. 290. 290. 290. 670. 300.]
# [590. 640. 570. 580. 290. 270. 310. 280. 570. 300.]
# [570. 160. 600. 570. 300. 290. 580. 580. 580. 300.]

def main():
    N = 10
    NA = .200 #Numerical Aperture
    W = wavelength / (pi * NA) #mode field diameter
    grating_length = 2000

    widths =    100 * np.ones(N)
    # lengths =   300 * np.ones(N)
    lengths = np.array([590., 630., 600., 580., 290., 290., 290., 290., 670., 300.]);

    # print scatter(100), scatter(100)**2
    #
    # quit();

    # for l in range(0, 520, 20):
    #     print propagate_matrix(l)
    # quit()

    amplitudes = getAmplitudes(widths, lengths)

    #gamma = getOverlap(widths, lengths, amplitudes, W, grating_length, N) #Test for getOverlap function

    start = time.time()

    #print(gamma) #Test for getOverlap function

    print(anneal(widths, lengths, W, grating_length, N))

    end = time.time()

    print("Run Time: ", end - start)
    
main()