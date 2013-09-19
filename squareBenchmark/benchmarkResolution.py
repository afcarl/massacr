# benchRes.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

res100 = [21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0]
ra100 = [3.642, 3.556, 3.475, 3.409, 3.341, 3.333, 3.306, 3.284, 3.266]
res1000 = [21.0, 31.0, 41.0, 51.0, 61.0]
ra1000 = [13.642, 13.556, 13.475, 13.409, 13.0]

fig=plt.figure()
ax1=fig.add_subplot(1,1,1)
r1 = plt.scatter(res, ra100, s=80, color='r', marker='^')
plt.plot(res,ra100,color='r')
#r2 = plt.scatter(res, ra1000, s=80, color='b', marker='^')
#plt.plot(res,ra1000,color='b')
plt.xlabel('LATERAL/VERTICAL GRID POINTS')
plt.xticks(res)
plt.ylabel('NUSSELT NUMBER')
plt.grid(b=True)
#plt.legend([r1, r2], {'Ra = 1000.0', 'Ra = 100.0'})

plt.savefig('benchmarkResolution.png')
