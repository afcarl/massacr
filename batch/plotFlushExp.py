# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 6.28e10     #
#############################
substep = 6.28e9

# load temps
t02 = np.loadtxt('t02s.txt')
t04 = np.loadtxt('t04s.txt')
t06 = np.loadtxt('t06s.txt')
t08 = np.loadtxt('t08s.txt')
t10 = np.loadtxt('t10s.txt')
t12 = np.loadtxt('t12s.txt')
t14 = np.loadtxt('t14s.txt')
t16 = np.loadtxt('t16s.txt')
t18 = np.loadtxt('t18s.txt')
t20 = np.loadtxt('t20s.txt')
temps = np.arange(2,22,2)

subs = range(9,100,1)

t = [[t02[subs,:]], [t04[subs,:]], [t06[subs,:]], [t08[subs,:]], [t10[subs,:]],
     [t12[subs,:]], [t14[subs,:]], [t16[subs,:]], [t18[subs,:]], [t20[subs,:]]]

# grab alkalinity for each experiment
t_alk = np.zeros((10))
t_glass = np.zeros((10))
t_kaolinite = np.zeros((10))
t_stilbite = np.zeros((10))
t_saponite = np.zeros((10))
t_albite = np.zeros((10))
t_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t_alk[i] = bit[0,-2,3]-.002407 # steady state alk
    #t_alk[i] = sum(bit[0,:,3]) # total alk produced for ocean
    t_glass[i] = bit[0,-2,56]
    t_kaolinite[i] = bit[0,-3,17]
    t_stilbite[i] = bit[0,-3,13]
    t_saponite[i] = bit[0,-3,21]
    t_albite[i] = bit[0,-3,19]
    t_celadonite[i] = bit[0,-3,23]




    



fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)



#########################
# PLOT ALK FOR MINERALS #
#########################

ax = plt.subplot(2,2,1)

##p = plt.plot(temps,t_alk/max(t_alk),'k^-',linewidth=2,label='flush 2kyr')
##p = plt.plot(temps,ts1_alk/max(ts1_alk),'r^-',linewidth=2,label='flush 1kyr')
##p = plt.plot(temps,ts250_alk/max(ts250_alk),'g^-',
##             linewidth=2,label='flush 250yr')

#p = plt.plot(temps,ts3_alk/2000.0,'b^-',linewidth=2,label='flush 3kyr')
p = plt.plot(temps,t_alk/2000.0,'k^-',linewidth=2,label='flush 2kyr')
p = plt.plot(temps,ts1_alk/2000.0,'m^-',linewidth=2,label='flush 1kyr')
p = plt.plot(temps,ts500_alk/2000.0,'y^-',
             linewidth=2,label='flush 500yr')
p = plt.plot(temps,ts250_alk/2000.0,'c^-',
             linewidth=2,label='flush 250yr')
p = plt.plot(temps,ts100_alk/2000.0,'g^-',
             linewidth=2,label='flush 100yr')

plt.title('ALK FLUX TO OCEAN FROM SILICATE DISSOLUTION + CLAY FORMATION',
          fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$ per 1000 mol basalt]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT GLASS FOR MINERALS #
###########################

##ax = plt.subplot(2,2,2)
##
##p = plt.plot(temps,t_glass,'k^-',linewidth=2,label='flush 2kyr')
##p = plt.plot(temps,ts1_glass,'m^-',linewidth=2,label='flush 1kyr')
##p = plt.plot(temps,ts500_glass,'y^-',linewidth=2,label='flush 500yr')
##p = plt.plot(temps,ts250_glass,'c^-',linewidth=2,label='flush 250yr')
##p = plt.plot(temps,ts100_glass,'g^-',linewidth=2,label='flush 100yr')
##
##
##plt.title('REMAINING BASALT',fontsize=8)
##plt.ylabel('ALKALINITY',fontsize=6)
##plt.xlabel('T [$^{\circ}$C]',fontsize=6)
##
##handles, labels = ax.get_legend_handles_labels()
##plt.legend(handles[::-1], labels[::-1])
##plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)




########################
# PLOT ALK NO MINERALS #
########################

ax = plt.subplot(2,2,4)

p = plt.plot(temps,tn_alk,'k^-',linewidth=2,label='flush 2kyr')
p = plt.plot(temps,tn1_alk,'m^-',linewidth=2,label='flush 1kyr')
p = plt.plot(temps,tn500_alk,'y^-',linewidth=2,label='flush 500yr')
p = plt.plot(temps,tn250_alk,'c^-',linewidth=2,label='flush 250yr')

plt.title('ALK FLUX ONLY FROM SILICATE DISSOLUTION',fontsize=10)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$ per 1000 mol basalt]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



##############
# PLOT CLAYS #
# 6.26e10    #
##############

ax = plt.subplot(2,2,3)

p = plt.plot(temps,t_celadonite,'k-',linewidth=3,label='flush 2kyr')
p = plt.plot(temps,2.0*ts1_celadonite,'m-',linewidth=3,label='flush 1kyr')
p = plt.plot(temps,8.0*ts500_celadonite,'y-',linewidth=3,label='flush 500yr')
p = plt.plot(temps,8.0*ts250_celadonite,'c-',linewidth=3,label='flush 250yr')
p = plt.plot(temps,20.0*ts100_celadonite,'g-',linewidth=3,label='flush 100yr')

##p = plt.plot(temps,t_stilbite,'k--',linewidth=3)
##p = plt.plot(temps,2.0*ts1_stilbite,'m--',linewidth=3)
##p = plt.plot(temps,8.0*ts500_stilbite,'y--',linewidth=3)
##p = plt.plot(temps,8.0*ts250_stilbite,'c--',linewidth=3)
##p = plt.plot(temps,20.0*ts100_stilbite,'g--',linewidth=3)

plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT ALK FOR ALK CHANGE #
###########################

##ax = plt.subplot(2,2,4)
##
##p = plt.plot(cs,ct10s1_alk/2000.0,'k^-',linewidth=2,label='flush 2kyr')
##
##
##plt.title('ALKALINITY FLUX TO OCEAN',fontsize=8)
##plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$ per 1000 mol basalt]',
##           fontsize=6)
##plt.xlabel('T [$^{\circ}$C]',fontsize=6)
##
##handles, labels = ax.get_legend_handles_labels()
##plt.legend(handles[::-1], labels[::-1])
##plt.legend(handles, labels,loc='best',prop={'size':8}, ncol=2)


plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
