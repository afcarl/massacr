#plotExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

# dic-limit, T=10
k200 = np.loadtxt('k200.txt')
k250 = np.loadtxt('k250.txt')
k400 = np.loadtxt('k400.txt')
k600 = np.loadtxt('k600.txt')

# alkalinity also increases in these experiments
alk_k = [k200[-1,2],k250[-1,2],k400[-1,2], k600[-1,2]]

# c-limit, T=10
c08t10 = np.loadtxt('c08t10.txt')
c10t10 = np.loadtxt('c10t10.txt')
c12t10 = np.loadtxt('c12t10.txt')
c14t10 = np.loadtxt('c14t10.txt')
c16t10 = np.loadtxt('c16t10.txt')
c18t10 = np.loadtxt('c18t10.txt')
c20t10 = np.loadtxt('c20t10.txt')
c22t10 = np.loadtxt('c22t10.txt')
c24t10 = np.loadtxt('c24t10.txt')
c26t10 = np.loadtxt('c26t10.txt')
c28t10 = np.loadtxt('c28t10.txt')
c30t10 = np.loadtxt('c30t10.txt')
c32t10 = np.loadtxt('c32t10.txt')

# carbon limit
c_t10 = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
         22.0, 24.0, 26.0, 28.0, 30.0, 32.0]

# final carbonate amount
caco3_t10 = [c08t10[-1,33], c10t10[-1,33], c12t10[-1,33],
             c14t10[-1,33],c16t10[-1,33],
             c18t10[-1,33], c20t10[-1,33], c22t10[-1,33],
             c24t10[-1,33], c26t10[-1,33], c28t10[-1,33],
             c30t10[-1,33], c32t10[-1,33]]

# final basalt almount
glass_t10 = [c08t10[-1,47], c10t10[-1,47], c12t10[-1,47],
             c14t10[-1,47],c16t10[-1,47],
             c18t10[-1,47], c20t10[-1,47], c22t10[-1,47],
             c24t10[-1,47], c26t10[-1,47], c28t10[-1,47],
             c30t10[-1,47], c32t10[-1,47]]

# final alk
alk_t10 = [c08t10[-1,2], c10t10[-1,2], c12t10[-1,2],
           c14t10[-1,2],c16t10[-1,2],
             c18t10[-1,2], c20t10[-1,2], c22t10[-1,2],
             c24t10[-1,2], c26t10[-1,2], c28t10[-1,2],
           c30t10[-1,2], c32t10[-1,2]]

# final ph
ph_t10 = [c08t10[-1,1], c10t10[-1,1], c12t10[-1,1],
          c14t10[-1,1],c16t10[-1,1],
             c18t10[-1,1], c20t10[-1,1], c22t10[-1,1],
             c24t10[-1,1], c26t10[-1,1], c28t10[-1,1],
          c30t10[-1,1], c32t10[-1,1]]


print alk_t10



fig=plt.figure()

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6)


##############
# PLOT CACO3 #
##############

ax = plt.subplot(2,2,1)
p = plt.plot(c_t10[2:],caco3_t10[2:],'o--')


plt.ylabel('CaCO$_3$ [mol]',fontsize=6)
plt.xlabel('MAX CO$_2$',fontsize=6)

###################
# PLOT ALKALINITY #
###################

ax = plt.subplot(2,2,2)
p = plt.plot(c_t10[2:],alk_t10[2:],'o--')


plt.ylabel('ALKALINITY [mol / kgw]',fontsize=6)
plt.xlabel('MAX CO$_2$',fontsize=6)

####################
# REMAINING BASALT #
####################

ax = plt.subplot(2,2,3)
p = plt.plot(c_t10[2:],glass_t10[2:],'o--')

plt.ylabel('REMAINING BASALTIC GLASS [mol]',fontsize=6)
plt.xlabel('MAX CO$_2$',fontsize=6)


######
# pH #
######

ax = plt.subplot(2,2,4)
p = plt.plot(c_t10[2:],ph_t10[2:],'o--')


plt.ylabel('pH ',fontsize=6)
plt.xlabel('MAX CO$_2$',fontsize=6)



plt.subplots_adjust(hspace=.4, wspace=.4)
plt.savefig('t10.png')
