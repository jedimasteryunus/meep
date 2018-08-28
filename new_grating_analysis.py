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

    real_t00 = data[11][index]
    imag_t00 = data[12][index]
    real_t01 = data[13][index]
    imag_t01 = data[14][index]

    real_t10 = data[23][index]
    imag_t10 = data[24][index]
    real_t11 = data[25][index]
    imag_t11 = data[26]index]

    t00 = complex(real_t00, imag_t00)
    t01 = complex(real_t01, imag_t01)

    t10 = complex(real_t10, imag_t10)
    t11 = complex(real_t11, imag_t11)

    return {"t00": t00, "t01": t01, "t10": t10, "t11": t11}

def reflection(width):
    index = widthToIndex(width)

    real_r00 = data[7][index]
    imag_r00 = data[8][index]
    real_r01 = data[9][index]
    imag_r01 = data[10][index]

    real_r10 = data[19][index]
    imag_r10 = data[20][index]
    real_r11 = data[21][index]
    imag_r11 = data[22]index]

    r00 = complex(real_r00, imag_r00)
    r01 = complex(real_r01, imag_r01)

    r10 = complex(real_r10, imag_r10)
    r11 = complex(real_r11, imag_r11)

    return {"r00": r00, "r01": r01, "r10": r10, "r11": r11}

def scatter(width):
    index = widthToIndex(width)

    real_su0 = data[15][index]
    imag_su0 = data[16][index]

    real_su1 = data[27][index]
    imag_su1 = data[28][index]

    su0 = complex(real_su0, imag_su0)

    su1 = complex(real_su1, imag_su1)

    return {"su0": su0. "su1": su1}

# MATRICES #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
def scatter_matrix(width):
    t = transmission(width)
    r = reflection(width)

    t00 = t["t00"]
    t01 = t["t01"]
    t10 = t["t10"]
    t11 = t["t11"]

    r00 = r["r00"]
    r01 = r["r01"]
    r10 = r["r10"]
    r11 = r["r11"]

    denom = t00 * t11 - t01 * t10

    return np.matrix([  [t11 / denom,
                        (t10 * np.conj(r01) - t11 * np.conj(r00)) / denom,
                        -t10 / denom,
                        (t10 * np.conj(r11) - t11 * np.conj(r10)) / denom],

                        [(t11 * r00 - t01 * r10) / denom,
                        np.conj(t00) + r00 * (t10 * np.conj(r01) - t11 * np.conj(r00)) / denom + r10 * (t01 + np.conj(r00) - t00 * np.conj(r01)) / denom,
                        (t00 * r01 - t10 * r00) / denom,
                        np.conj(t10) + r00 * (t10 * np.conj(r11) - t11 * np.conj(r10)) / denom + r10 * (t01 * np.conj(r10) - t00 * np.conj(r11)) / denom],

                        [-t01 / denom,
                        (t01 * np.conj(r00) - t00 * np.conj(r01)) / denom,
                        t00 / denom,
                        (t01 * np.conj(r10) - t00 * np.conj(r11)) / denom],

                        [(t11 * r01 - t01 * r11) / denom,
                        np.cong(t01) + r01 * (t10 * np.conj(r01) - t11 * np.conj(r00)) / denom + r11 * (t01 * np.conj(r00) - t00 * np.conj(r01)) / denom,
                        (t00 * r11 - t10 * r01) / denom,
                        np.conj(t11) + r01 * (t10 * np.conj(r11) - t11 * np.conj(r10)) / denom + r11 * (t01 * np.conj(r10) - t00 * np.conj(r11)) / denom] ],

                        dtype=complex)

def propagate_matrix(length):
    n_eff_fund = data[31][0]
    n_eff_first = data[32][0]
    return np.matrix([  [np.exp(-pi*2j*n_eff_fund*length/wavelength), 0, 0, 0],
                        [0, np.exp(pi*2j*n_eff_fund*length/wavelength), 0, 0],
                        [0, 0, np.exp(-pi*2j*n_eff_first*length/wavelength), 0],
                        [0, 0, 0, np.exp(pi*2j*n_eff_first*length/wavelength)]  ],
                        dtype=complex)

# AMPLITUDES #### #### #### #### #### #### #### #### #### ####
def getAmplitudes(widths, lengths):
    N = widths.shape[0]                                 # N grates.

    a = np.zeros((4, 2*N), dtype=complex)                              # N times [[ a  b  ]; [ b' a' ]].

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
            if np.log10(T) == -2.0:
                print("Progress:  0.0%")
            else:
                print("Progress: ", -(np.log10(T) + 2) * 100, "%")

            new_lengths = neighbor(lengths)
            new_gtrx = getOverlap(widths, new_lengths, getAmplitudes(widths, new_lengths), W, grating_length, num_notches)

            ap = acceptance_probability(gtrx[0], new_gtrx[0], T)
            r = random();

            if ap > r:
                if gtrx[0] < new_gtrx[0]:
                    print("Gamma INCREASED this step",)
                else:
                    print("Gamma DECREASED this step",)

                lengths = new_lengths
                gtrx = new_gtrx
            else:
                print("Gamma DID NOT CHANGE this step",)

            print(    "S = {:.2f}%,\tT = {:.2f}%,\tR = {:.2f}%,\tX = {};".format(gtrx[0]*100, gtrx[1]*100, gtrx[2]*100, gtrx[3]))
            print(    "S'= {:.2f}%,\tT'= {:.2f}%,\tR'= {:.2f}%,\tX'= {}.".format(new_gtrx[0]*100, new_gtrx[1]*100, new_gtrx[2]*100, new_gtrx[3]), "\n")

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

    widths = 100 * np.ones(N)
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

    anneal_tup = anneal(widths, lengths, W, grating_length, N)

    print("Lengths: ", anneal_tup[0])
    print("GTRX List: ", anneal_tup[1])

    end = time.time()

    print("Run Time: ", end - start)

main()
