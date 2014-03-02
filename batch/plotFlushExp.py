# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,42,2)
m=-8 # alk/other step
n=-8 # precip step




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 350 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c475p.txt')
t04 = np.loadtxt('r90t04c475p.txt')
t06 = np.loadtxt('r90t06c475p.txt')
t08 = np.loadtxt('r90t08c475p.txt')
t10 = np.loadtxt('r90t10c475p.txt')
t12 = np.loadtxt('r90t12c475p.txt')
t14 = np.loadtxt('r90t14c475p.txt')
t16 = np.loadtxt('r90t16c475p.txt')
t18 = np.loadtxt('r90t18c475p.txt')
t20 = np.loadtxt('r90t20c475p.txt')
t22 = np.loadtxt('r90t22c475p.txt')
t24 = np.loadtxt('r90t24c475p.txt')
t26 = np.loadtxt('r90t26c475p.txt')
t28 = np.loadtxt('r90t28c475p.txt')
t30 = np.loadtxt('r90t30c475p.txt')
t32 = np.loadtxt('r90t32c475p.txt')
t34 = np.loadtxt('r90t34c475p.txt')
t36 = np.loadtxt('r90t36c475p.txt')
t38 = np.loadtxt('r90t38c475p.txt')
t40 = np.loadtxt('r90t40c475p.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90c475_alk = np.zeros((len(temps)))
t90c475_alkflux = np.zeros((len(temps)))
t90c475_alkflux0 = np.zeros((len(temps)))
t90c475_glass = np.zeros((len(temps)))
t90c475_water = np.zeros((len(temps)))
t90c475_HCO3 = np.zeros((len(temps)))
t90c475_CO3 = np.zeros((len(temps)))


t90c475_kaolinite = np.zeros((len(temps)))
t90c475_stilbite = np.zeros((len(temps)))
t90c475_saponite = np.zeros((len(temps)))
t90c475_albite = np.zeros((len(temps)))
t90c475_celadonite = np.zeros((len(temps)))
t90c475_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90c475_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
    t90c475_alkflux[i] = gx[m,3]*1000.0
    t90c475_alkflux0[i] = bit[0,m,3]*1000.0 
    t90c475_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90c475_glass[i] = gx[n,59]

    t90c475_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90c475_HCO3[i] = bit[0,m,13]
    t90c475_CO3[i] = bit[0,m,14]

    t90c475_kaolinite[i] = gx[n,19]
    t90c475_stilbite[i] = gx[n,15]
    t90c475_saponite[i] = gx[n,23]
    t90c475_albite[i] = gx[n,21]
    t90c475_celadonite[i] = gx[n,25]
    t90c475_quartz[i] = gx[n,47]

##    t90c475_kaolinite[i] = np.max(abs(gx[:,19]))
##    t90c475_stilbite[i] = np.max(abs(gx[:,15]))
##    t90c475_saponite[i] = np.max(abs(gx[:,23]))
##    t90c475_albite[i] = np.max(abs(gx[:,21]))
##    t90c475_celadonite[i] = np.max(abs(gx[:,25]))
##    t90c475_quartz[i] = np.max(abs(gx[:,47]))


# H+ concentration change from minerals
t90c475_dH_clay = ((t90c475_kaolinite*6.0 + t90c475_stilbite*8.72 + t90c475_saponite*7.32 + \
          t90c475_albite*4.0 + t90c475_celadonite*6.0) / t90c475_water) 
# H+ consumption by basalt dissolution
t90c475_dH_diss = -t90c475_glass * .5 / t90c475_water







#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 350 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c350p.txt')
t04 = np.loadtxt('r90t04c350p.txt')
t06 = np.loadtxt('r90t06c350p.txt')
t08 = np.loadtxt('r90t08c350p.txt')
t10 = np.loadtxt('r90t10c350p.txt')
t12 = np.loadtxt('r90t12c350p.txt')
t14 = np.loadtxt('r90t14c350p.txt')
t16 = np.loadtxt('r90t16c350p.txt')
t18 = np.loadtxt('r90t18c350p.txt')
t20 = np.loadtxt('r90t20c350p.txt')
t22 = np.loadtxt('r90t22c350p.txt')
t24 = np.loadtxt('r90t24c350p.txt')
t26 = np.loadtxt('r90t26c350p.txt')
t28 = np.loadtxt('r90t28c350p.txt')
t30 = np.loadtxt('r90t30c350p.txt')
t32 = np.loadtxt('r90t32c350p.txt')
t34 = np.loadtxt('r90t34c350p.txt')
t36 = np.loadtxt('r90t36c350p.txt')
t38 = np.loadtxt('r90t38c350p.txt')
t40 = np.loadtxt('r90t40c350p.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90_alk = np.zeros((len(temps)))
t90_alkflux = np.zeros((len(temps)))
t90_alkflux0 = np.zeros((len(temps)))
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
    t90_alkflux[i] = gx[m,3]*1000.0 
    t90_alkflux0[i] = bit[0,m,3]*1000.0
    t90_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90_glass[i] = gx[n,59]
    
    t90_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90_HCO3[i] = bit[0,m,13]
    t90_CO3[i] = bit[0,m,14]

    t90_kaolinite[i] = gx[n,19]
    t90_stilbite[i] = gx[n,15]
    t90_saponite[i] = gx[n,23]
    t90_albite[i] = gx[n,21]
    t90_celadonite[i] = gx[n,25]
    t90_quartz[i] = gx[n,47]

##    t90_kaolinite[i] = np.max(abs(gx[:,19]))
##    t90_stilbite[i] = np.max(abs(gx[:,15]))
##    t90_saponite[i] = np.max(abs(gx[:,23]))
##    t90_albite[i] = np.max(abs(gx[:,21]))
##    t90_celadonite[i] = np.max(abs(gx[:,25]))
##    t90_quartz[i] = np.max(abs(gx[:,47]))


# H+ concentration change from minerals
t90_dH_clay = ((t90_kaolinite*6.0 + t90_stilbite*8.72 + t90_saponite*7.32 + \
          t90_albite*4.0 + t90_celadonite*6.0) / t90_water) 
# H+ consumption by basalt dissolution
t90_dH_diss = -t90_glass * .5 / t90_water





#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 225 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c225p.txt')
t04 = np.loadtxt('r90t04c225p.txt')
t06 = np.loadtxt('r90t06c225p.txt')
t08 = np.loadtxt('r90t08c225p.txt')
t10 = np.loadtxt('r90t10c225p.txt')
t12 = np.loadtxt('r90t12c225p.txt')
t14 = np.loadtxt('r90t14c225p.txt')
t16 = np.loadtxt('r90t16c225p.txt')
t18 = np.loadtxt('r90t18c225p.txt')
t20 = np.loadtxt('r90t20c225p.txt')
t22 = np.loadtxt('r90t22c225p.txt')
t24 = np.loadtxt('r90t24c225p.txt')
t26 = np.loadtxt('r90t26c225p.txt')
t28 = np.loadtxt('r90t28c225p.txt')
t30 = np.loadtxt('r90t30c225p.txt')
t32 = np.loadtxt('r90t32c225p.txt')
t34 = np.loadtxt('r90t34c225p.txt')
t36 = np.loadtxt('r90t36c225p.txt')
t38 = np.loadtxt('r90t38c225p.txt')
t40 = np.loadtxt('r90t40c225p.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90c225_alk = np.zeros((len(temps)))
t90c225_alkflux = np.zeros((len(temps)))
t90c225_alkflux0 = np.zeros((len(temps)))
t90c225_glass = np.zeros((len(temps)))
t90c225_water = np.zeros((len(temps)))
t90c225_HCO3 = np.zeros((len(temps)))
t90c225_CO3 = np.zeros((len(temps)))


t90c225_kaolinite = np.zeros((len(temps)))
t90c225_stilbite = np.zeros((len(temps)))
t90c225_saponite = np.zeros((len(temps)))
t90c225_albite = np.zeros((len(temps)))
t90c225_celadonite = np.zeros((len(temps)))
t90c225_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90c225_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
    t90c225_alkflux[i] = gx[m,3]*1000.0
    t90c225_alkflux0[i] = bit[0,m,3]*1000.0
    t90c225_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90c225_glass[i] = gx[n,59]
    
    t90c225_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90c225_HCO3[i] = bit[0,m,13]
    t90c225_CO3[i] = bit[0,m,14]

    t90c225_kaolinite[i] = gx[n,19]
    t90c225_stilbite[i] = gx[n,15]
    t90c225_saponite[i] = gx[n,23]
    t90c225_albite[i] = gx[n,21]
    t90c225_celadonite[i] = gx[n,25]
    t90c225_quartz[i] = gx[n,47]

##    t90c225_kaolinite[i] = np.max(abs(gx[:,19]))
##    t90c225_stilbite[i] = np.max(abs(gx[:,15]))
##    t90c225_saponite[i] = np.max(abs(gx[:,23]))
##    t90c225_albite[i] = np.max(abs(gx[:,21]))
##    t90c225_celadonite[i] = np.max(abs(gx[:,25]))
##    t90c225_quartz[i] = np.max(abs(gx[:,47]))


# H+ concentration change from minerals
t90c225_dH_clay = ((t90c225_kaolinite*6.0 + t90c225_stilbite*8.72 + t90c225_saponite*7.32 + \
          t90c225_albite*4.0 + t90c225_celadonite*6.0) / t90c225_water) 
# H+ consumption by basalt dissolution
t90c225_dH_diss = -t90c225_glass * .5 / t90c225_water







#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 100 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c100p.txt')
t04 = np.loadtxt('r90t04c100p.txt')
t06 = np.loadtxt('r90t06c100p.txt')
t08 = np.loadtxt('r90t08c100p.txt')
t10 = np.loadtxt('r90t10c100p.txt')
t12 = np.loadtxt('r90t12c100p.txt')
t14 = np.loadtxt('r90t14c100p.txt')
t16 = np.loadtxt('r90t16c100p.txt')
t18 = np.loadtxt('r90t18c100p.txt')
t20 = np.loadtxt('r90t20c100p.txt')
t22 = np.loadtxt('r90t22c100p.txt')
t24 = np.loadtxt('r90t24c100p.txt')
t26 = np.loadtxt('r90t26c100p.txt')
t28 = np.loadtxt('r90t28c100p.txt')
t30 = np.loadtxt('r90t30c100p.txt')
t32 = np.loadtxt('r90t32c100p.txt')
t34 = np.loadtxt('r90t34c100p.txt')
t36 = np.loadtxt('r90t36c100p.txt')
t38 = np.loadtxt('r90t38c100p.txt')
t40 = np.loadtxt('r90t40c100p.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90c100_alk = np.zeros((len(temps)))
t90c100_alkflux = np.zeros((len(temps)))
t90c100_alkflux0 = np.zeros((len(temps)))
t90c100_glass = np.zeros((len(temps)))
t90c100_water = np.zeros((len(temps)))
t90c100_HCO3 = np.zeros((len(temps)))
t90c100_CO3 = np.zeros((len(temps)))


t90c100_kaolinite = np.zeros((len(temps)))
t90c100_stilbite = np.zeros((len(temps)))
t90c100_saponite = np.zeros((len(temps)))
t90c100_albite = np.zeros((len(temps)))
t90c100_celadonite = np.zeros((len(temps)))
t90c100_quartz = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90c100_alk[i] = np.max(bit[0,:,3])
    gx, gy = np.gradient(bit[0,:,:])
    t90c100_alkflux[i] = gx[m,3]*1000.0
    t90c100_alkflux0[i] = bit[0,m,3]*1000.0
    t90c100_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
    t90c100_glass[i] = gx[n,59]
    
    t90c100_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
    t90c100_HCO3[i] = bit[0,m,13]
    t90c100_CO3[i] = bit[0,m,14]

    t90c100_kaolinite[i] = gx[n,19]
    t90c100_stilbite[i] = gx[n,15]
    t90c100_saponite[i] = gx[n,23]
    t90c100_albite[i] = gx[n,21]
    t90c100_celadonite[i] = gx[n,25]
    t90c100_quartz[i] = gx[n,47]

##    t90c100_kaolinite[i] = np.max(abs(gx[:,19]))
##    t90c100_stilbite[i] = np.max(abs(gx[:,15]))
##    t90c100_saponite[i] = np.max(abs(gx[:,23]))
##    t90c100_albite[i] = np.max(abs(gx[:,21]))
##    t90c100_celadonite[i] = np.max(abs(gx[:,25]))
##    t90c100_quartz[i] = np.max(abs(gx[:,47]))


# H+ concentration change from minerals
t90c100_dH_clay = ((t90c100_kaolinite*6.0 + t90c100_stilbite*8.72 + t90c100_saponite*7.32 + \
          t90c100_albite*4.0 + t90c100_celadonite*6.0) / t90c100_water) 
# H+ consumption by basalt dissolution
t90c100_dH_diss = -t90c100_glass * .5 / t90c100_water


















fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

########################
# PLOT ALK FLUX (DIFF) #
########################

ax = plt.subplot(1,1,1)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,t90c475_alkflux,'c^-',linewidth=2,label='475ppm CO$_2$')
p = plt.plot(temps,t90_alkflux,'k^-',linewidth=2,label='350ppm CO$_2$')
p = plt.plot(temps,t90c225_alkflux,'b^-',linewidth=2,label='225ppm CO$_2$')
p = plt.plot(temps,t90c100_alkflux,'r^-',linewidth=2,label='100ppm CO$_2$')


#plt.text(6, .023, r"+ ALK TO OCEAN", horizontalalignment='center', fontsize=12)
#plt.text(23, -.02, r"- ALK TO OCEAN", horizontalalignment='center', fontsize=12)

plt.text(7, .03, r"+ SYSTEM PRODUCES ALKALINITY", horizontalalignment='center', fontsize=8)
plt.text(23, -.02, r"- SYSTEM CONSUMES ALKALINITY", horizontalalignment='center', fontsize=8)

plt.title('ALKALINITY FLUX',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':10}, ncol=1)


plt.savefig('flushAlk.png')









fig=plt.figure()

plt.rc('xtick', labelsize=7) 
plt.rc('ytick', labelsize=7)



##################
# PLOT dH COMBOS #
##################

ax = plt.subplot(2,2,1)

print t90_dH_clay
print t90_dH_diss
# clay >> diss

p = plt.plot(temps[0:20],t90c475_dH_clay[0:20]-t90c475_dH_diss[0:20],
             'c-',linewidth=1,label='p475 both')
p = plt.plot(temps[0:20],t90_dH_clay[0:20]-t90_dH_diss[0:20],
             'k-',linewidth=1,label='p350')
p = plt.plot(temps[0:20],t90c225_dH_clay[0:20]-t90c225_dH_diss[0:20],
             'b-',linewidth=1,label='p225')
p = plt.plot(temps[0:20],t90c100_dH_clay[0:20]-t90c100_dH_diss[0:20],
             'r-',linewidth=1,label='p100')

##p = plt.plot(temps,-t90c475_dH_diss,'c-',linewidth=1,label='p475 diss')
##p = plt.plot(temps,-t90_dH_diss,'k-',linewidth=1,label='p350 diss')
##p = plt.plot(temps,-t90c225_dH_diss,'b-',linewidth=1,label='p225')
##p = plt.plot(temps,-t90c100_dH_diss,'r-',linewidth=1,label='p100')


plt.title('PRODUCTION - H+ CONSUMPTION',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)




##################
# PLOT dH COMBOS #
##################


ax = plt.subplot(2,2,2)

p = plt.plot(temps,t90c475_dH_clay,'c--',linewidth=1,label='p475 clay')
p = plt.plot(temps,t90_dH_clay,'k--',linewidth=1,label='p350 clay')
p = plt.plot(temps,t90c225_dH_clay,'b--',linewidth=1,label='p225')
p = plt.plot(temps,t90c100_dH_clay,'r--',linewidth=1,label='p100')

##print t90c475_HCO3
##p = plt.plot(temps,t90c475_HCO3+2.0*t90c475_CO3,'c^-',linewidth=1,label='p475 clay')
##p = plt.plot(temps,t90_HCO3+2.0*t90_CO3,'k^-',linewidth=1,label='p350 clay')
##p = plt.plot(temps,t90c225_HCO3+2.0*t90c225_CO3,'b^-',linewidth=1,label='p225')
##p = plt.plot(temps,t90c100_HCO3+2.0*t90c100_CO3,'r^-',linewidth=1,label='p100')

plt.title('H+ PRODUCTION',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT GLASS FOR MINERALS #
###########################

ax = plt.subplot(2,2,4)

p = plt.plot(temps,-t90c475_glass,'c^-',linewidth=2,label='p350')
p = plt.plot(temps,-t90_glass,'k^-',linewidth=2,label='p350')
p = plt.plot(temps,-t90c225_glass,'b^-',linewidth=2,label='p225')
p = plt.plot(temps,-t90c100_glass,'r^-',linewidth=2,label='p100')


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

p = plt.plot(temps,t90c475_albite,'c-.',linewidth=1,label='albite')
p = plt.plot(temps,t90c225_albite,'b-.',linewidth=1,label='albite')
p = plt.plot(temps,t90c100_albite,'r-.',linewidth=1,label='albite')

p = plt.plot(temps,t90c475_saponite,'c-',linewidth=2)
p = plt.plot(temps,t90c225_saponite,'b-',linewidth=2)
p = plt.plot(temps,t90c100_saponite,'r-',linewidth=2)

p = plt.plot(temps,t90c475_kaolinite,'c:',linewidth=2)
p = plt.plot(temps,t90c225_kaolinite,'b:',linewidth=2)
p = plt.plot(temps,t90c100_kaolinite,'r:',linewidth=2)

##p = plt.plot(temps,t90c100_celadonite,'r-',linewidth=1,label='100ppm')
##p = plt.plot(temps,t90c100_stilbite,'r--',linewidth=1)
##p = plt.plot(temps,t90c100_kaolinite,'r:',linewidth=2)
##p = plt.plot(temps,t90c100_albite,'r-.',linewidth=1)
##p = plt.plot(temps,t90c100_saponite,'r-',linewidth=2)


plt.title('mineral production rate', fontsize=8)
plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
