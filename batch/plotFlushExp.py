# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,42,2)


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 99/01        #
# A.K.A.                    #
# t_RES = 3.14e12 (100kyr)  #
#############################

# load temps
t02 = np.loadtxt('r99t02.txt')
t04 = np.loadtxt('r99t04.txt')
t06 = np.loadtxt('r99t06.txt')
t08 = np.loadtxt('r99t08.txt')
t10 = np.loadtxt('r99t10.txt')
t12 = np.loadtxt('r99t12.txt')
t14 = np.loadtxt('r99t14.txt')
t16 = np.loadtxt('r99t16.txt')
t18 = np.loadtxt('r99t18.txt')
t20 = np.loadtxt('r99t20.txt')
t22 = np.loadtxt('r99t22.txt')
t24 = np.loadtxt('r99t24.txt')
t26 = np.loadtxt('r99t26.txt')
t28 = np.loadtxt('r99t28.txt')
t30 = np.loadtxt('r99t30.txt')
t32 = np.loadtxt('r99t32.txt')
t34 = np.loadtxt('r99t34.txt')
t36 = np.loadtxt('r99t36.txt')
t38 = np.loadtxt('r99t38.txt')
t40 = np.loadtxt('r99t40.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t99_alk = np.zeros((len(temps)))
t99_alkflux = np.zeros((len(temps)))
t99_glass = np.zeros((len(temps)))
t99_kaolinite = np.zeros((len(temps)))
t99_stilbite = np.zeros((len(temps)))
t99_saponite = np.zeros((len(temps)))
t99_albite = np.zeros((len(temps)))
t99_celadonite = np.zeros((len(temps)))
t99_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t99_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
##    t90_glass[i] = bit[0,-2,57]
    t99_alkflux[i] = np.max(abs(gx[:,3]))
    t99_glass[i] = np.max(abs(gx[:,59]))
    t99_kaolinite[i] = np.max(abs(gx[:,19]))
    t99_stilbite[i] = np.max(abs(gx[:,15]))
    t99_saponite[i] = np.max(abs(gx[:,23]))
    t99_albite[i] = np.max(abs(gx[:,21]))
    t99_celadonite[i] = np.max(abs(gx[:,25]))
    t99_quartz[i] = np.max(abs(gx[:,47]))




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
    t90_alkflux[i] = np.max(abs(gx[:,3])) # [mol / kgw / kyr]
    t90_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90_water[i] = bit[0,np.max(abs(gx[:,59])),63]
    t90_HCO3[i] = bit[0,np.max(abs(gx[:,59])),13]
    t90_CO3[i] = bit[0,np.max(abs(gx[:,59])),14]

    t90_kaolinite[i] = np.max(abs(gx[:,19]))
    t90_stilbite[i] = np.max(abs(gx[:,15]))
    t90_saponite[i] = np.max(abs(gx[:,23]))
    t90_albite[i] = np.max(abs(gx[:,21]))
    t90_celadonite[i] = np.max(abs(gx[:,25]))
    t90_quartz[i] = np.max(abs(gx[:,47]))

# units: mol ca2+ kyr^-1 kgw^-1
##t90_ALKdiss = (t90_glass * .0069) / t90_water
##t90_dalkflux = np.gradient(t90_alkflux)
##t90_dALKdiss = 2.0 * np.gradient(t90_ALKdiss)
##t90_dALKclay = (t90_dalkflux) - t90_dALKdiss


# TRYING FROM THE OTHER ANGLE

# change in concentration [mol/kgw] of H+ per timestep
t90_dH = (t90_kaolinite*6.0 + t90_stilbite*8.72 + t90_saponite*7.32 + \
          t90_albite*4.0 + t90_celadonite*6.0) / t90_water
print "t90_dH"
print t90_dH
print " "
# change in concentration [mol/kgw] of HCO3- per timestep
t90_dHCO3clay = (t90_CO3 / 4.69e-11) * t90_dH
t90_dALKclay = t90_dHCO3clay
t90_dALKdiss = t90_alkflux - t90_dALKclay
print "t90_alkflux (total?)"
print t90_alkflux
print " "
print "t90_dALKclay"
print t90_dALKclay
print " "
print "t90_dALKdiss"
print t90_dALKdiss
print " "






fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)



##################
# PLOT ALK FINAL #
##################

ax = plt.subplot(2,2,1)


p = plt.plot(temps,t90_alk,'k^-',linewidth=2,label='mixing 90/10')

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


p = plt.plot(temps,t90_alkflux,'k^-',linewidth=2,label='mixing 90/10')

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

p = plt.plot(temps,t90_glass,'k^-',linewidth=2,label='mixing 90/10')


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
