# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,12,2)


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 99/01        #
#############################

# load temps
t02 = np.loadtxt('r99t02.txt')
t04 = np.loadtxt('r99t04.txt')
t06 = np.loadtxt('r99t06.txt')
t08 = np.loadtxt('r99t08.txt')
t10 = np.loadtxt('r99t10.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]]]

##t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
##     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
##     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
##     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90_alk = np.zeros((len(temps)))
t90_glass = np.zeros((len(temps)))
t90_kaolinite = np.zeros((len(temps)))
t90_stilbite = np.zeros((len(temps)))
t90_saponite = np.zeros((len(temps)))
t90_albite = np.zeros((len(temps)))
t90_celadonite = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90_alk[i] = bit[0,-2,3] - .002#*.1
    t90_glass[i] = bit[0,-2,57]
    t90_kaolinite[i] = bit[0,-3,18]
    t90_stilbite[i] = bit[0,-3,14]
    t90_saponite[i] = bit[0,-3,22]
    t90_albite[i] = bit[0,-3,20]
    t90_celadonite[i] = bit[0,-3,24]




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 99/01        #
#############################

# load temps
t02 = np.loadtxt('r99t02a4.txt')
t04 = np.loadtxt('r99t04a4.txt')
t06 = np.loadtxt('r99t06a4.txt')
t08 = np.loadtxt('r99t08a4.txt')
t10 = np.loadtxt('r99t10a4.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]]]

##t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
##     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
##     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
##     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

# grab for each experiment
t90a4_alk = np.zeros((len(temps)))
t90a4_glass = np.zeros((len(temps)))
t90a4_kaolinite = np.zeros((len(temps)))
t90a4_stilbite = np.zeros((len(temps)))
t90a4_saponite = np.zeros((len(temps)))
t90a4_albite = np.zeros((len(temps)))
t90a4_celadonite = np.zeros((len(temps)))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90a4_alk[i] = bit[0,-2,3] - .002#*.1
    t90a4_glass[i] = bit[0,-2,57]
    t90a4_kaolinite[i] = bit[0,-3,18]
    t90a4_stilbite[i] = bit[0,-3,14]
    t90a4_saponite[i] = bit[0,-3,22]
    t90a4_albite[i] = bit[0,-3,20]
    t90a4_celadonite[i] = bit[0,-3,24]







fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)



#########################
# PLOT ALK FOR MINERALS #
#########################

ax = plt.subplot(2,2,1)


p = plt.plot(temps,t90_alk,'k^-',linewidth=2,label='mixing 99/01')

plt.title('ALKALINITY',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$ per 1000 mol basalt]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT GLASS FOR MINERALS #
###########################

ax = plt.subplot(2,2,2)

p = plt.plot(temps,t90_glass,'k^-',linewidth=2,label='mixing 99/01')


plt.title('REMAINING BASALT',fontsize=8)
plt.ylabel('ALKALINITY',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=6)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)




########################
# PLOT ALK NO MINERALS #
########################

ax = plt.subplot(2,2,4)


plt.title('',fontsize=10)
plt.ylabel('',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

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


##p = plt.plot(temps,t60c16_celadonite,'g-',linewidth=1,label='mixing 60/40, alk 1.6e-3')
##p = plt.plot(temps,t60c16_stilbite,'g--',linewidth=1)
##p = plt.plot(temps,t60c16_kaolinite,'g:',linewidth=2)
##p = plt.plot(temps,t60c16_albite,'g-.',linewidth=1)
##p = plt.plot(temps,t60c16_saponite,'g-',linewidth=2)


##p = plt.plot(temps,t60_celadonite,'c-',linewidth=1,label='mixing 60/40, alk 2.0e-3')
##p = plt.plot(temps,t60_stilbite,'c--',linewidth=1)
##p = plt.plot(temps,t60_kaolinite,'c:',linewidth=2)
##p = plt.plot(temps,t60_albite,'c-.',linewidth=1)
##p = plt.plot(temps,t60_saponite,'c-',linewidth=2)


#p = plt.plot(temps,t80_celadonite,'r-',linewidth=3,label='mixing 80/20')
#p = plt.plot(temps,t80_stilbite,'r--',linewidth=3)


plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT 2,2 #
###########################

ax = plt.subplot(2,2,4)



plt.title('',fontsize=8)
plt.ylabel('',fontsize=6)
plt.xlabel('',fontsize=6)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':8}, ncol=2)


plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
