import ast
import matplotlib.pyplot as plt

output_file = "notch.txt"

f = open(output_file, "r")
list_of_lines = []
for line in f.readlines():
    line = ast.literal_eval(line)
    list_of_lines.append(line)
f.close()

ws = list_of_lines[0]
Rs = list_of_lines[1]
Ts = list_of_lines[2]
Ss = list_of_lines[3]
Sus = list_of_lines[4]
Sds = list_of_lines[5]
norm_Sus = list_of_lines[6]

plt.plot(ws, Rs,'bo-',label='reflectance')
plt.plot(ws, Ts,'ro-',label='transmittance')
plt.plot(ws, Ss,'go-',label='net loss')
plt.plot(ws, Sus, 'co-', label='upper loss')
plt.plot(ws, Sds, 'mo-', label='lower loss')
plt.plot(ws, norm_Sus, 'ko-', label='normalized upper loss')
plt.axis([40.0, 300.0, 0.0, 100.0])
plt.xlabel("wavelength (nm)")
plt.ylabel("percentage (%)")
plt.legend(loc="center right")
plt.show()
