from math import pi, sqrt, atan
from random import random, randint, choice
import time
import ast
import csv
import numpy as np
import matplotlib.pyplot as plt
import cProfile

wavelength = 637
#wavelength = 420
# neff = 2.2012 # refractive index inside waveguide (at wavelength 0.6372)

# N = 20
# z = 20e3

# N = 15
# z = 10e3
# NA = .150 #Numerical Aperture

N = 8
z = 2.5e3
NA = .300 #Numerical Aperture

w = wavelength / (pi * NA) #mode field diameter

k = 2*pi/wavelength

zr = pi*w*w/wavelength

W = w*sqrt(1 + (z/zr)**2)
R = z*(1 + (zr/z)**2)
phi = atan(z/zr)

debug = 0

input_file = 'transpose-notch-LN637.csv'
notch_file = open(input_file, 'r')

width_dict = dict()
data_list = list(csv.reader(notch_file, delimiter=','))

n_eff = float(data_list[-1][0])
print(n_eff, type(n_eff))

for row in data_list[:-1]:
	width_dict[float(row[0])] = dict()
	width_dict[float(row[0])]['r00'] = complex(float(row[7]), float(row[8]))
	width_dict[float(row[0])]['t00'] = complex(float(row[9]), float(row[10]))
	width_dict[float(row[0])]['sd0'] = complex(float(row[11]), float(row[12]))
	width_dict[float(row[0])]['su0'] = complex(float(row[13]), float(row[14]))

output_file = 'new_new_grating_analysis.out'
f = open(output_file, 'w')

# DATA #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# def widthToIndex(width):
#     # return int(width / 20 - 2)
#     # return 1
#     return 2

# def importData(fname):
#     f = open(fname, "r")
#     data = []
#     for line in f.readlines():
#             line = ast.literal_eval(line)
#             data.append(line)
#     f.close()
#     return data

# fname = "notch2.txt"
# fname = "AlN-a=16/notch-AlN420-test.txt"
# data = importData(fname)

def widthCheck(width):
	if width in width_dict.keys():
		return True
	else:
		return False

def transmission(width):
	if not widthCheck(width):
		print("Alert! The width you have entered may not have any MEEP data.")

	# index = widthToIndex(width)
	#
	# real_t00 = data[11][index]
	# imag_t00 = data[12][index]
	# real_t01 = data[13][index]
	# imag_t01 = data[14][index]
	#
	# real_t10 = data[23][index]
	# imag_t10 = data[24][index]
	# real_t11 = data[25][index]
	# imag_t11 = data[26][index]
	#
	# t00 = complex(real_t00, imag_t00)
	# t01 = complex(real_t01, imag_t01)
	#
	# t10 = complex(real_t10, imag_t10)
	# t11 = complex(real_t11, imag_t11)
	#
	# return {"t00": t00, "t01": t01, "t10": t10, "t11": t11}

	# t00 = complex(0.927068099, -1.334914783)
	# t00 = 0.927068099 * np.exp(-1.334914783j);
	t00 = width_dict[width]["t00"]

	return {"t00": t00}

def reflection(width):

	if not widthCheck(width):
		print("Alert! The width you have entered may not have any MEEP data.")

	# index = widthToIndex(width)
	#
	# real_r00 = data[7][index]
	# imag_r00 = data[8][index]
	# real_r01 = data[9][index]
	# imag_r01 = data[10][index]
	#
	# real_r10 = data[19][index]
	# imag_r10 = data[20][index]
	# real_r11 = data[21][index]
	# imag_r11 = data[22][index]
	#
	# r00 = complex(real_r00, imag_r00)
	# r01 = complex(real_r01, imag_r01)
	#
	# r10 = complex(real_r10, imag_r10)
	# r11 = complex(real_r11, imag_r11)
	#
	# return {"r00": r00, "r01": r01, "r10": r10, "r11": r11}

	# r00 = complex(0.083756091, 0.311361372)
	# r00 = 0.083756091 * np.exp(0.311361372j);
	r00 = width_dict[width]["r00"]

	return {"r00": r00}

def scatter(width):

	if not widthCheck(width):
		print("Alert! The width you have entered may not have any MEEP data.")

	# index = widthToIndex(width)
	#
	# real_su0 = data[15][index]
	# imag_su0 = data[16][index]
	# real_su1 = data[27][index]
	# imag_su1 = data[28][index]
	#
	# real_sd0 = data[17][index]
	# imag_sd0 = data[18][index]
	# real_sd1 = data[29][index]
	# imag_sd1 = data[30][index]
	#
	# su0 = complex(real_su0, imag_su0)
	# su1 = complex(real_su1, imag_su1)
	# sd0 = complex(real_sd0, imag_sd0)
	# sd1 = complex(real_sd1, imag_sd1)
	#
	# return {"su0": su0, "su1": su1, "sd0": sd0, "sd1": sd1}

	# su0 = complex(0.098703181, 1.141884602)
	# sd0 = complex(0.243604087, 1.141884602)
	# su0 = 0.098703181 * np.exp(1.141884602j);
	# sd0 = 0.243604087 * np.exp(1.141884602j);

	su0 = width_dict[width]["su0"]
	sd0 = width_dict[width]["sd0"]

	return {"su0": su0, "sd0": sd0}

# MATRICES #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

def scatter_matrix(width):

	if not widthCheck(width):
		print("Alert! The width you have entered may not have any MEEP data.")

	t00 = width_dict[width]["t00"]
	r00 = width_dict[width]["r00"]

	return np.matrix([  [ 1/t00,     np.conj(r00)/t00             ],
						[ r00/t00,   (r00*np.conj(r00) + t00*np.conj(t00))/t00   ]  ],
						dtype=complex)

def propagate_matrix(length):
	return np.matrix([  [ np.exp(pi*2j*n_eff*length/wavelength),              0                     ],
						[ 0,                                 np.exp(-pi*2j*n_eff*length/wavelength)   ]  ],
						dtype=complex)

	# return np.matrix([  [np.exp(          pi * 2j * n_eff_fund  * length / wavelength), 0, 0, 0],
	#                     [0, np.exp(      -pi * 2j * n_eff_fund  * length / wavelength),    0, 0],
	#                     [0, 0, np.exp(    pi * 2j * n_eff_first * length / wavelength),       0],
	#                     [0, 0, 0, np.exp(-pi * 2j * n_eff_first * length / wavelength)         ]  ],
	#                     dtype=complex)

# AMPLITUDES #### #### #### #### #### #### #### #### #### ####
def getAmplitudes(widths, lengths):
	N = widths.shape[0]                                 # N grates.

	a_fund = np.zeros((2, 2*N), dtype=complex)                              # N times [[ a  b  ]; [ b' a' ]].

	a_fund[0, 2*N-1] = 1                                     # Set the initial vector (b = 1, a' = 0).
	a_fund[:, 2*N-2] = a_fund[:, 2*N-1] * scatter_matrix(widths[N-1])   # And find a, b' for the first grate.

	# print(scatter_matrix(widths[0]))

	for ii in range(N-2, -1, -1):                               # Now do this for the rest of the grates.
		a_fund[:, 2*ii+1] = propagate_matrix(lengths[ii]).dot(a_fund[:, 2*ii+2])
		a_fund[:, 2*ii] =   scatter_matrix(widths[ii])   .dot(a_fund[:, 2*ii+1])

	# a_first = np.zeros((4, 2*N), dtype=complex)                              # N times [[ a  b  ]; [ b' a' ]].
	#
	# a_first[2, 2*N-1] = 1                                     # Set the initial vector (b = 1, a' = 0).
	# a_first[:, 2*N-2] = a_first[:, 2*N-1] * scatter_matrix(widths[N-1])   # And find a, b' for the first grate.
	#
	# for ii in range(N-2, -1, -1):                               # Now do this for the rest of the grates.
	#     # a_first[:, 2*ii+1] = a_first[:, 2*ii+2] * propagate_matrix(lengths[ii])
	#     # a_first[:, 2*ii] =   a_first[:, 2*ii+1]  * scatter_matrix(widths[ii])
	#     a_first[:, 2*ii+1] = propagate_matrix(lengths[ii]).dot(a_first[:, 2*ii+2])
	#     a_first[:, 2*ii] =   scatter_matrix(widths[ii])   .dot(a_first[:, 2*ii+1])
	#
	# a = a_fund - (a_fund[2,0]/a_first[2,0]) * a_first

	# print(np.abs(a_fund[2,0]/a_first[2,0]))

	return a_fund


# OVERLAP #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# def E(x,w,z):
#     k = 2*pi/wavelength
#
#     zr = pi*w*w/wavelength
#
#     W = w*sqrt(1 + (z/zr)**2)
#     R = z*(1 + (zr/z)**2)
#     phi = atan(z/zr)
#
#     return np.exp(-x**2 / W**2) * np.exp(-1j *(k*z + k*x*x/(2*R) - phi) )

def E(x,w,z):
	# return np.exp(-x**2 / W**2) * np.exp(-1j * (k*z + k*x*x/(2*R) - phi) )
	return np.exp(-x**2 / W**2) * np.exp(1j * k*x*x/(2*R) )
	# d = (2 * pi / wavelength) * (400/800);
	# return d if x == 0 else ( d * np.sin(x * d) / (pi * x * d)) * np.exp(-1j *(k*z + k*x*x/(2*R) - phi) )
	# return ( 50000 * np.sin(x / d) / (pi * x / d)) * np.exp(-1j *(k*z + k*x*x/(2*R) - phi) )

def getOverlap(widths, lengths, amplitudes, W, Z, num_notches, outputScatter):
	other_notches = 5;

	tot_notches = num_notches + 2*other_notches

	other_notch_space = 200

	s =     np.zeros(tot_notches, dtype = complex)
	sd =    np.zeros(tot_notches, dtype = complex)
	x =     np.zeros(tot_notches)
	dx =    np.ones(tot_notches)*other_notch_space

	#for i in range(0, num_notches):
		#print("Scatter:", scatter(widths[i]))
		#print("Amplitudes:", amplitudes[0,2*i] + amplitudes[1,2*i+1])

	currentx = -widths[0]/2;

	for i in range(0, num_notches):
		scatter_dict = scatter(widths[i])
		s[i+other_notches] =  scatter_dict["su0"]*(amplitudes[0,2*i] + amplitudes[1,2*i+1]) #+ scatter_dict["su1"]*(amplitudes[2,2*i] + amplitudes[3,2*i+1])
		sd[i+other_notches] = scatter_dict["sd0"]*(amplitudes[0,2*i] + amplitudes[1,2*i+1]) #+ scatter_dict["sd1"]*(amplitudes[2,2*i] + amplitudes[3,2*i+1])
		x[i+other_notches] = currentx + widths[i]/2;

		if 		i == 0:
			dx[i+other_notches] = widths[i] + lengths[i]
		elif 	i == num_notches-1:
			dx[i+other_notches] = widths[i] + lengths[i-1]
		else:
			dx[i+other_notches] = widths[i] + (lengths[i-1] + lengths[i])/2

		if 	i != num_notches-1:
			currentx += lengths[i] + widths[i];

	for i in range(0, other_notches):
		x[other_notches-i-1] = x[other_notches-i] - other_notch_space - 100;
		x[other_notches+num_notches+i] = x[other_notches+num_notches-1+i] + other_notch_space + 100;

	s  /= amplitudes[0,0]
	sd /= amplitudes[0,0]

	# if debug:
	#     S = np.zeros(num_notches)
	#
	#     for i in range(0, num_notches):
	#         S[i] = 100 * np.abs(s[i])**2
	#
	#     print(S)


	scatter_up = 0
	scatter_down = 0

	for i in range(0, tot_notches):
		scatter_up      += np.abs( s[i])**2
		scatter_down    += np.abs(sd[i])**2

	# np.set_printoptions(suppress=True)
	# np.set_printoptions(precision=3)
	# print(np.abs(amplitudes/amplitudes[0,0])**2)
	# print(np.angle(amplitudes/amplitudes[0,0]))

	final_reflection =          np.abs(amplitudes[1,0]                  /amplitudes[0,0])**2
	final_transmission =        np.abs(amplitudes[0,2*num_notches-1]    /amplitudes[0,0])**2
	# final_reflection_first =    np.abs(amplitudes[3,0]                  /amplitudes[0,0])**2
	# final_transmission_first =  np.abs(amplitudes[2,2*num_notches-1]    /amplitudes[0,0])**2
	final_reflection_first =    0
	final_transmission_first =  0

	final_gamma = 0
	final_X = 0;
	final_scatter_up = 0
	final_scatter_down = 0
	sensitivity = 100;
	Wround = int(min(currentx, W)/200)*100;
	# print(range(Wround, int(np.max(x)) + sensitivity-Wround, sensitivity));

	gammav2 = True

	for X in range(Wround, int(currentx) + sensitivity-Wround, sensitivity):
		gamma = 0
		integral = 0

		for i in range(0, tot_notches):
			if gammav2:
				# gamma           += s[i] * np.conj(E(x[i] - X, W, Z)) * dx[i]
				# integral        += (np.abs(E(x[i] - X, W, Z))**2) * (dx[i]**2)
				gamma           += s[i] * np.conj(E(x[i] - X, W, Z)) * sqrt(dx[i])
				integral        += (np.abs(E(x[i] - X, W, Z))**2) * dx[i]
			else:
				gamma           += s[i] * np.conj(E(x[i] - X, W, Z))
				integral        += (np.abs(E(x[i] - X, W, Z))**2)

		gamma = (np.abs(gamma/integral)**2)
		# gamma = (np.abs(gamma)**2)

		# if outputScatter:
		# 	print("{:.2f}\t{:.2f}\t{:.2f}".format(X, 100*gamma, 100*integral))

		if gamma > final_gamma:
			final_gamma = gamma
			final_X = X;

	if outputScatter:
		print(range(Wround, int(currentx) + sensitivity-Wround, sensitivity))

		gamma = 0
		integral = 0

		for i in range(0, tot_notches):
			if gammav2:
				# gamma           += s[i] * np.conj(E(x[i] - X, W, Z)) * dx[i]
				# integral        += (np.abs(E(x[i] - X, W, Z))**2) * (dx[i]**2)
				gamma           += s[i] * np.conj(E(x[i] - X, W, Z)) * sqrt(dx[i])
				integral        += (np.abs(E(x[i] - X, W, Z))**2) * dx[i]
			else:
				gamma           += s[i] * np.conj(E(x[i] - X, W, Z))
				integral        += (np.abs(E(x[i] - X, W, Z))**2)

		# for i in range(0, tot_notches):
		# 	# print("{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(x[i], 100*(np.abs(s[i])**2)/dx[i], np.angle(s[i]), 100*(np.abs(E(x[i] - final_X, W))**2)/dx[i]/integral))
		# 	# print("{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(x[i], 100*(np.abs(s[i])**2) / dx[i], np.angle(s[i]), 100*(np.abs(E(x[i] - final_X, W))**2)*dx[i]/integral))
		# 	print("{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}".format(x[i], 100*(np.abs(s[i])**2), np.angle(s[i]), 100*(np.abs(E(x[i] - final_X, W, Z))**2)/integral))

		# print(np.angle(np.transpose(amplitudes)/amplitudes[0,0]));

		# plt.plot(x, 100*(np.abs(s)**2)/dx,'bo-',label='Scatter')
		# plt.plot(x, 1*(np.angle(s) + pi),'ro-',label='Phase')
		# plt.plot(x, 100*(np.abs(E(x - final_X, W))**2)/dx/integral,'go-',label='Match')
		plt.gcf().clear()
		plt.plot(x, 500*(np.abs(s)**2),'bo-',label='Scatter')
		plt.plot(x, 5000*(np.abs(s)**2)/sqrt(dx[i]),'co-',label='Scatter')
		plt.plot(x, 1*(np.angle(s) + pi),'ro-',label='Phase')
		plt.plot(x, 10*(np.abs(E(x - final_X, W, Z))**2),'go-',label='Match')
		plt.plot(x, 1*(np.angle(E(x - final_X, W, Z)) + pi),'ko-',label='Match')
		plt.plot(x, dx/100,'yo-',label='Match')
		plt.draw()
		plt.pause(.00001)

		# print [X, gamma]
	# print final_X;

	# total_scatter_up =      np.sum(np.abs(s)**2);
	# total_scatter_down =    np.sum(np.abs(sd)**2);

	return [final_gamma, final_transmission, final_reflection, final_transmission_first, final_reflection_first, final_X, scatter_up, scatter_down] #, total_scatter_up, total_scatter_down]

# ANNEALING #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
def neighbor(lengths):
	N = len(lengths)
	sensitivity = 1
	rand_index = randint(0, N - 1)
	new_lengths = lengths.copy()

	lowerbound = 50
	upperbound = 400
	# if rand_index < 5:
	#     upperbound = 900

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

def anneal(widths, lengths, W, Z, num_notches):
	gtrx = getOverlap(widths, lengths, getAmplitudes(widths, lengths), W, Z, num_notches, False)

	# T = 0.001
	# T_min = 0.0001
	# T = 0.001
	# T_min = 0.0001
	T =     0.000001
	T_min = 0.000000001
	alpha = 0.975

	# T = T_min/alpha;

	output = False;

	while T > T_min:
		i = 1

		while i <= 200:
			new_lengths = neighbor(lengths)
			new_gtrx = getOverlap(widths, new_lengths, getAmplitudes(widths, new_lengths), W, Z, num_notches, output or i == 1)

			ap = acceptance_probability(gtrx[0], new_gtrx[0], T)
			r = random();

			if ap > r:
				lengths = new_lengths
				gtrx = new_gtrx

			if output or i == 1:

				print("Lengths:           ", list(lengths))
				#f.write("Lengths:           %s" % (lengths))

				# print("Grate positions:   ", x)
				#f.write("Grate positions:   %s" % (x))

				if ap > r:
					if gtrx[0] < new_gtrx[0]:
						print("Gamma INCREASED this step",)
					else:
						print("Gamma DECREASED this step",)
				else:
					print("Gamma DID NOT CHANGE this step",)

				print("γ = {:.2f}%,\tT0 = {:.2f}%,\tR0 = {:.2f}%,\tT1 = {:.2f}%,\tR1 = {:.2f}%,\tS  = {:.2f}%,\tP =  {:.2f}%,\tX  = {};".format(gtrx[0]*100, gtrx[1]*100, gtrx[2]*100, gtrx[3]*100, gtrx[4]*100, gtrx[6]*100, (gtrx[1] + gtrx[2] + gtrx[3] + gtrx[4] + gtrx[6] + gtrx[7])*100, gtrx[5]))
				# print("S = {:.2f}%,\tT0 = {:.2f}%,\tR0 = {:.2f}%,\tT1 = {:.2f}%,\tR1 = {:.2f}%,\tX  = {},\tSu = {:.2f}%,\tSd = {:.2f}%;".format(gtrx[0]*100, gtrx[1]*100, gtrx[2]*100, gtrx[3]*100, gtrx[4]*100, gtrx[5], gtrx[6]*100, gtrx[7]*100))
				#f.write("S = {:.2f}%,\tT = {:.2f}%,\tR = {:.2f}%,\tX = {};".format(gtrx[0]*100, gtrx[1]*100, gtrx[2]*100, gtrx[3]))

				print("γ'= {:.2f}%,\tT0'= {:.2f}%,\tR0'= {:.2f}%,\tT1'= {:.2f}%,\tR1'= {:.2f}%,\tS' = {:.2f}%,\tP'=  {:.2f}%,\tX' = {}.".format(new_gtrx[0]*100, new_gtrx[1]*100, new_gtrx[2]*100, new_gtrx[3]*100, new_gtrx[4]*100, new_gtrx[6]*100, (new_gtrx[1] + new_gtrx[2] + new_gtrx[3] + new_gtrx[4] + new_gtrx[6] + new_gtrx[7])*100, new_gtrx[5]))
				# print("S'= {:.2f}%,\tT0'= {:.2f}%,\tR0'= {:.2f}%,\tT1'= {:.2f}%,\tR1'= {:.2f}%,\tX' = {},\tSu'= {:.2f}%,\tSd'= {:.2f}%.\n".format(new_gtrx[0]*100, new_gtrx[1]*100, new_gtrx[2]*100, new_gtrx[3]*100, new_gtrx[4]*100, new_gtrx[5], new_gtrx[6]*100, new_gtrx[7]*100))
				#f.write("S'= {:.2f}%,\tT'= {:.2f}%,\tR'= {:.2f}%,\tX'= {}.".format(new_gtrx[0]*100, new_gtrx[1]*100, new_gtrx[2]*100, new_gtrx[3]))

'''
				if np.log10(T) == -2.0:
					print("Progress:  0.0%")
				else:
					print("Progress: ", -(np.log10(T) + 2) * 100, "%")
'''
			i += 1

		T = T*alpha

	return lengths, gtrx

# MAIN #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

# [590. 630. 600. 580. 290. 290. 290. 290. 670. 300.]
# [590. 640. 570. 580. 290. 270. 310. 280. 570. 300.]
# [570. 160. 600. 570. 300. 290. 580. 580. 580. 300.]
# [640. 540. 330. 170. 600.]

def main():
	# print(transmission(100))
	# print(reflection(100))
	# print(scatter(100))
	# print(np.angle(scatter_matrix(100)))
	# for l in range(0, 350, 50):
	#     print(l)
	#     print(np.angle(propagate_matrix(l)))
	# #
	# print(np.angle(propagate_matrix(300)))
	# return

	# t = transmission(100);
	# r = reflection(100);
	# s = scatter(100);
	#
	# T0 = np.abs(t["t00"])**2 + np.abs(t["t01"])**2;
	# R0 = np.abs(r["r00"])**2 + np.abs(r["r01"])**2;
	# Su0 = np.abs(s["su0"])**2
	# Sd0 = np.abs(s["sd0"])**2
	#
	# T1 = np.abs(t["t10"])**2 + np.abs(t["t11"])**2;
	# R1 = np.abs(r["r10"])**2 + np.abs(r["r11"])**2;
	# Su1 = np.abs(s["su1"])**2
	# Sd1 = np.abs(s["sd1"])**2
	#
	# print(T0, R0, Su0, Sd0, T0+R0+Su0+Sd0)
	# print(T1, R1, Su1, Sd1, T1+R1+Su1+Sd1)
	#
	# return

	plt.axis([0, 7000, 0.0, 20.0])
	plt.xlabel("Position (nm)")
	plt.ylabel("scatter (%)")
	# plt.legend(loc="center right")
	plt.ion()
	plt.show()

	# grating_length = 2000

	# widths =    90 * np.ones(N)
	widths =    100 * np.ones(N)
	# widths =    50 * np.ones(N)
	# lengths =   300 * np.ones(N)
	# lengths = [490.0, 315.0, 315.0, 340.0, 305.0, 460.0, 385.0, 290.0, 335.0, 480.0, 290.0, 350.0, 475.0, 480.0, 170.0, 305.0, 330.0, 155.0, 325.0, 300.0]

	# lengths = np.array([590., 630., 600., 580., 290., 290., 290., 290., 670., 300.]);
	# lengths = [395, 395, 395, 395, 395, 395, 395, 395, 395, 395, 300, 500, 455, 325, 475, 335, 470, 325, 280, 300]
	# lengths = [250, 220, 255, 210, 255, 190, 270, 200, 275, 185, 275, 175, 285, 175, 280, 205, 270, 190, 280, 300]
	# lengths = 100 * np.ones(N)
	# lengths = [500, 205, 310, 405, 500, 400, 320, 400, 400, 400, 310, 410, 405, 405, 405, 405, 405, 405, 445, 100];
	# lengths = [500, 205, 310, 405, 500, 400, 320, 400, 400, 400, 310, 410, 185, 185, 185, 185, 185, 185, 225, 100];
	# lengths = [165, 355, 165, 370, 155, 170, 380, 155, 170, 390, 180, 140, 185, 400, 400, 405, 195, 405, 190, 100];
	# lengths = [355, 555, 365, 160, 365, 175, 155, 165, 170, 165, 390, 160, 165, 195, 380, 415, 395, 405, 485, 100];
	# lengths = [355, 555, 365, 160, 365, 175, 155, 165, 170, 165, 390, 160, 165, 195, 380];
	# lengths = [155, 155, 160, 155, 165, 165, 155, 175, 180, 150, 175, 400, 395, 180, 380];
	# lengths = [160, 155, 160, 165, 160, 175, 150, 170, 395, 400, 400, 410, 400, 380, 380];
	# lengths = [160, 155, 160, 165, 160, 175, 150];
	# lengths = [110, 120, 135, 125, 175, 210, 255];
	# lengths = [110, 340, 130, 135, 160, 165, 375];
	# lengths = [320, 345, 125, 140, 145, 145, 120];
	# lengths = [310, 340, 120, 125, 145, 135, 105];
	# lengths = [115, 340, 130, 135, 160, 165];
	# lengths = [120, 135, 135, 165, 185, 165, 165];
	# lengths = [120, 110, 275, 100, 275, 220, 165];
	# lengths = [125, 500, 145, 135, 165, 155, 170];
	# lengths = [120, 500, 145, 130, 165, 155];
	# lengths = [120, 300, 145, 130];
	lengths = [143, 313, 328, 135, 132, 167, 165];
	lengths = [181, 152, 307, 100, 259, 100, 199];
	lengths = [181, 152, 307, 100, 259, 100];
	lengths = [177, 145, 113, 100, 259, 100];
	lengths = [177, 145, 113, 100, 259, 100, 100];
	lengths = [174, 145, 316, 100, 267, 100, 389]
	lengths = [173, 146, 500, 93,  87,  85];
	lengths = [158, 335, 105, 275, 90];
	lengths = [160, 340, 299, 95, 281, 299, 280];
	# widths =  [50,  75,  100, 100, 100, 75, 50, 50];
	# lengths = [600, 500, 400, 400, 300, 300, 200, 200, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
	# lengths = [600, 500, 400, 400, 300, 300, 200, 200, 100, 100]

	# lengths = [205, 230, 360, 390, 385, 375, 285, 150, 125, 155, 310, 255, 130, 170, 345, 245, 395, 320, 145, 100]


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

	anneal_tup = anneal(widths, lengths, w, z, N)

	print("Lengths: ", anneal_tup[0])
	f.write("Lengths: %s \n" % (anneal_tup[0]))

	print("GTRX List: ", anneal_tup[1])
	f.write("GTRX List: %s \n" % (anneal_tup[1]))

	end = time.time()

	print("Run Time: ", end - start)
	f.write("Run Time: %s \n" % (end-start))

# cProfile.run('main()')
main()
f.close()
