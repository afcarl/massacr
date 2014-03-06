#
# secondary.py
#
# plot logK(T) for secondary minerals
#
# log K = a + bT + c/T + dlog10(T) + e/T2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,42,2) + 273.15

albite = [-1.1694e+001, 1.4429e-002, 1.3784e+004, -7.2866e+000, -1.6136e+006]
kaolinite =[1.6835e+001, -7.8939e-003, 7.7636e+003, -1.2190e+001, -3.2354e+005]
saponite = [9.8888e+000, 1.4320e-002, 1.9418e+004, -1.5259e+001, -1.3716e+006]
stilbite = [-2.4483e+001, 3.0987e-002, 2.8013e+004, -1.5802e+001, -3.4491e+006]
celadonite = [-3.3097e+001, 1.7989e-002, 1.8919e+004, -2.1219e+000, -2.0588e+006]

albite_k = np.zeros((len(temps)))
kaolinite_k = np.zeros((len(temps)))
saponite_k = np.zeros((len(temps)))
stilbite_k = np.zeros((len(temps)))
celadonite_k = np.zeros((len(temps)))

for i in range(len(temps)):
    albite_k[i] = albite[0] + albite[1]*temps[i] + albite[2]/temps[i] + \
               albite[3]*np.log10(temps[i]) + albite[4]/(temps[i]**2)
    kaolinite_k[i] = kaolinite[0] + kaolinite[1]*temps[i] + kaolinite[2]/temps[i] + \
               kaolinite[3]*np.log10(temps[i]) + kaolinite[4]/(temps[i]**2)
    saponite_k[i] = saponite[0] + saponite[1]*temps[i] + saponite[2]/temps[i] + \
               saponite[3]*np.log10(temps[i]) + saponite[4]/(temps[i]**2)
    stilbite_k[i] = stilbite[0] + stilbite[1]*temps[i] + stilbite[2]/temps[i] + \
               stilbite[3]*np.log10(temps[i]) + stilbite[4]/(temps[i]**2)
    celadonite_k[i] = celadonite[0] + celadonite[1]*temps[i] + celadonite[2]/temps[i] + \
               celadonite[3]*np.log10(temps[i]) + celadonite[4]/(temps[i]**2)

print temps.shape
print albite_k.shape


fig=plt.figure()

ax = plt.subplot(1,1,1)
p = plt.plot(temps-273.15,albite_k, label='albite')
p = plt.plot(temps-273.15,kaolinite_k, label='kaolinite')
p = plt.plot(temps-273.15,saponite_k, label='saponite-mg')
p = plt.plot(temps-273.15,stilbite_k, label='stilbite')
p = plt.plot(temps-273.15,celadonite_k, label='celadonite')

plt.title('SECONDARY MINERAL STABILITY',fontsize=10)
plt.ylabel('log$_{10}$(K)',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':10}, ncol=1)


plt.savefig('secondary.png')
