# benchRes.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

res100 = [21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0]
ra100 = [3.642, 3.556, 3.475, 3.409, 3.341, 3.333, 3.306]

res1000 = [41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0, 110.0]
ra1000 = [13.176197503855343,
13.800794074121205,
14.215223171531104,
14.415501741798980,
14.572451065421175,
14.672856820478509,
14.734969790697797,
14.770670602292542]


plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)

################
#   ra = 100   #
################

fig=plt.figure()
ax1=fig.add_subplot(1,2,1)
r1 = plt.scatter(res100, ra100, s=80, color='r', marker='^')
plt.plot(res100,ra100,color='r')
plt.xlabel('LATERAL/VERTICAL GRID POINTS')
plt.xticks(res100)
plt.ylabel('NUSSELT NUMBER')
plt.grid(b=True)

################
#   ra = 1000  #
################

ax1=fig.add_subplot(1,2,2)
r1 = plt.scatter(res1000, ra1000, s=80, color='b', marker='^')
plt.plot(res1000,ra1000,color='b')
plt.xlabel('LATERAL/VERTICAL GRID POINTS')
plt.xticks(res1000)
plt.ylabel('NUSSELT NUMBER')
plt.grid(b=True)

plt.savefig('benchRes.png')
