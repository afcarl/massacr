# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,42,2)



#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
#############################

# load temps
t02 = np.loadtxt('r90t02.txt')
t04 = np.loadtxt('r90t04.txt')
t06 = np.loadtxt('r90t06.txt')
t08 = np.loadtxt('r90t08.txt')
t10 = np.loadtxt('r90t10.txt')
t12 = np.loadtxt('r90t12.txt')
t14 = np.loadtxt('r90t14.txt')
t16 = np.loadtxt('r90t16.txt')
t18 = np.loadtxt('r90t18.txt')
t20 = np.loadtxt('r90t20.txt')
t22 = np.loadtxt('r90t22.txt')
t24 = np.loadtxt('r90t24.txt')
t26 = np.loadtxt('r90t26.txt')
t28 = np.loadtxt('r90t28.txt')
t30 = np.loadtxt('r90t30.txt')
t32 = np.loadtxt('r90t32.txt')
t34 = np.loadtxt('r90t34.txt')
t36 = np.loadtxt('r90t36.txt')
t38 = np.loadtxt('r90t38.txt')
t40 = np.loadtxt('r90t40.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90_alk = np.zeros((len(temps)))
t90_alkflux = np.zeros((len(temps)))
t90_glass = np.zeros((len(temps)))
t90_water = np.zeros((len(temps)))
t90_HCO3 = np.zeros((len(temps)))
t90_CO3 = np.zeros((len(temps)))


t90_kaolinite = np.zeros((len(temps)))
t90_stilbite = np.zeros((len(temps)))
t90_saponite = np.zeros((len(temps)))
t90_albite = np.zeros((len(temps)))
t90_celadonite = np.zeros((len(temps)))
t90_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
    t90_alkflux[i] = gx[-1,3] # [mol / kgw / kyr]
    t90_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90_HCO3[i] = gx[np.argmax(abs(gx[:,25])),13]
    t90_CO3[i] = gx[np.argmax(abs(gx[:,25])),14]

    t90_kaolinite[i] = np.max(abs(gx[:,19]))
    t90_stilbite[i] = np.max(abs(gx[:,15]))
    t90_saponite[i] = np.max(abs(gx[:,23]))
    t90_albite[i] = np.max(abs(gx[:,21]))
    t90_celadonite[i] = np.max(abs(gx[:,25]))
    t90_quartz[i] = np.max(abs(gx[:,47]))

t90_alkflux0 = 2*t90_CO3 + t90_HCO3

# H+ concentration change from minerals
t90_dH_clay = ((t90_kaolinite*6.0 + t90_stilbite*8.72 + t90_saponite*7.32 + \
          t90_albite*4.0 + t90_celadonite*6.0) / t90_water) 
print t90_dH_clay
# H+ consumption by basalt dissolution
t90_dH_diss = -t90_glass * .5 / t90_water
print t90_dH_diss




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 1e-4.45            #
#############################

# load temps
t02 = np.loadtxt('r90t02p445.txt')
t04 = np.loadtxt('r90t04p445.txt')
t06 = np.loadtxt('r90t06p445.txt')
t08 = np.loadtxt('r90t08p445.txt')
t10 = np.loadtxt('r90t10p445.txt')
t12 = np.loadtxt('r90t12p445.txt')
t14 = np.loadtxt('r90t14p445.txt')
t16 = np.loadtxt('r90t16p445.txt')
t18 = np.loadtxt('r90t18p445.txt')
t20 = np.loadtxt('r90t20p445.txt')
t22 = np.loadtxt('r90t22p445.txt')
t24 = np.loadtxt('r90t24p445.txt')
t26 = np.loadtxt('r90t26p445.txt')
t28 = np.loadtxt('r90t28p445.txt')
t30 = np.loadtxt('r90t30p445.txt')
t32 = np.loadtxt('r90t32p445.txt')
t34 = np.loadtxt('r90t34p445.txt')
t36 = np.loadtxt('r90t36p445.txt')
t38 = np.loadtxt('r90t38p445.txt')
t40 = np.loadtxt('r90t40p445.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90p4_alk = np.zeros((len(temps)))
t90p4_alkflux = np.zeros((len(temps)))
t90p4_glass = np.zeros((len(temps)))
t90p4_water = np.zeros((len(temps)))
t90p4_HCO3 = np.zeros((len(temps)))
t90p4_CO3 = np.zeros((len(temps)))


t90p4_kaolinite = np.zeros((len(temps)))
t90p4_stilbite = np.zeros((len(temps)))
t90p4_saponite = np.zeros((len(temps)))
t90p4_albite = np.zeros((len(temps)))
t90p4_celadonite = np.zeros((len(temps)))
t90p4_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90p4_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
    t90p4_alkflux[i] = gx[-1,3] # [mol / kgw / kyr]
    t90p4_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90p4_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90p4_HCO3[i] = gx[np.argmax(abs(gx[:,25])),13]
    t90p4_CO3[i] = gx[np.argmax(abs(gx[:,25])),14]

    t90p4_kaolinite[i] = np.max(abs(gx[:,19]))
    t90p4_stilbite[i] = np.max(abs(gx[:,15]))
    t90p4_saponite[i] = np.max(abs(gx[:,23]))
    t90p4_albite[i] = np.max(abs(gx[:,21]))
    t90p4_celadonite[i] = np.max(abs(gx[:,25]))
    t90p4_quartz[i] = np.max(abs(gx[:,47]))

t90p4_alkflux0 = 2*t90p4_CO3 + t90p4_HCO3

# H+ concentration change from minerals
t90p4_dH_clay = ((t90p4_kaolinite*6.0 + t90p4_stilbite*8.72 + t90p4_saponite*7.32 + \
          t90p4_albite*4.0 + t90p4_celadonite*6.0) / t90p4_water) 
print t90_dH_clay
# H+ consumption by basalt dissolution
t90p4_dH_diss = -t90p4_glass * .5 / t90p4_water
print t90p4_dH_diss











fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)



##################
# PLOT ALK FINAL #
##################

ax = plt.subplot(2,2,1)


#p = plt.plot(temps,t90_dH_diss,'r-',linewidth=1,label='ALKdiss')
#p = plt.plot(temps,-t90_dH_clay,'b-',linewidth=1,label='ALKclay')

p = plt.plot(temps,t90_dH_clay+t90_dH_diss,'k-',linewidth=1,label='ALKclay')
p = plt.plot(temps,t90p4_dH_clay+t90p4_dH_diss,'r-',linewidth=1,label='ALKclay')

plt.title('ALKALINITY',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



########################
# PLOT ALK FLUX (DIFF) #
########################

ax = plt.subplot(2,2,2)


p = plt.plot(temps,t90_alkflux,'k^-',linewidth=2,label='p3')
p = plt.plot(temps,t90p4_alkflux,'r^-',linewidth=2,label='p4')

plt.title('ALKALINITY FLUX',fontsize=8)
plt.ylabel('',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT GLASS FOR MINERALS #
###########################

ax = plt.subplot(2,2,4)

p = plt.plot(temps,t90_glass,'k^-',linewidth=2,label='p3')
p = plt.plot(temps,t90p4_glass,'r^-',linewidth=2,label='p4')


plt.title('BASALT DISSOLUTION RATE',fontsize=8)
plt.ylabel('[mol kyr$^{-1}$]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=6)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





##############
# PLOT CLAYS #
##############

ax = plt.subplot(2,2,3)

p = plt.plot(temps,t90_celadonite,'k-',linewidth=1,label='celadonite')
p = plt.plot(temps,t90_stilbite,'k--',linewidth=1,label='stilbite')
p = plt.plot(temps,t90_kaolinite,'k:',linewidth=2,label='kaolinite')
p = plt.plot(temps,t90_albite,'k-.',linewidth=1,label='albite')
p = plt.plot(temps,t90_saponite,'k-',linewidth=2,label='saponite')


plt.title('mineral production rate', fontsize=8)
plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
