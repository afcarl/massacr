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






#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 3.14e10     #
#############################
substep = 3.14e9

# load temps
t02s1 = np.loadtxt('t02s1kyr.txt')
t04s1 = np.loadtxt('t04s1kyr.txt')
t06s1 = np.loadtxt('t06s1kyr.txt')
t08s1 = np.loadtxt('t08s1kyr.txt')
t10s1 = np.loadtxt('t10s1kyr.txt')
t12s1 = np.loadtxt('t12s1kyr.txt')
t14s1 = np.loadtxt('t14s1kyr.txt')
t16s1 = np.loadtxt('t16s1kyr.txt')
t18s1 = np.loadtxt('t18s1kyr.txt')
t20s1 = np.loadtxt('t20s1kyr.txt')
temps = np.arange(2,22,2)

ts1 = [[t02s1[subs,:]], [t04s1[subs,:]], [t06s1[subs,:]], [t08s1[subs,:]],
       [t10s1[subs,:]], [t12s1[subs,:]], [t14s1[subs,:]], [t16s1[subs,:]],
       [t18s1[subs,:]], [t20s1[subs,:]]]

# grab alkalinity for each experiment
ts1_alk = np.zeros((10))
ts1_glass = np.zeros((10))
ts1_kaolinite = np.zeros((10))
ts1_stilbite = np.zeros((10))
ts1_saponite = np.zeros((10))
ts1_albite = np.zeros((10))
ts1_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(ts1[i])
    ts1_alk[i] = 2.0*bit[0,-2,3]-2.0*.002407 # steady state alk
    #ts1_alk[i] = 2.0*sum(bit[0,:,3]) # total alk produced for ocean
    ts1_glass[i] = bit[0,-2,56]
    ts1_kaolinite[i] = bit[0,-3,17]
    ts1_stilbite[i] = bit[0,-3,13]
    ts1_saponite[i] = bit[0,-3,21]
    ts1_albite[i] = bit[0,-3,19]
    ts1_celadonite[i] = bit[0,-3,23]



#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 9.42e10     #
#############################
substep = 9.42e9

# load temps
t02s3 = np.loadtxt('t02s3kyr.txt')
t04s3 = np.loadtxt('t04s3kyr.txt')
t06s3 = np.loadtxt('t06s3kyr.txt')
t08s3 = np.loadtxt('t08s3kyr.txt')
t10s3 = np.loadtxt('t10s3kyr.txt')
t12s3 = np.loadtxt('t12s3kyr.txt')
t14s3 = np.loadtxt('t14s3kyr.txt')
t16s3 = np.loadtxt('t16s3kyr.txt')
t18s3 = np.loadtxt('t18s3kyr.txt')
t20s3 = np.loadtxt('t20s3kyr.txt')
temps = np.arange(2,22,2)

ts3 = [[t02s3[subs,:]], [t04s3[subs,:]], [t06s3[subs,:]], [t08s3[subs,:]],
       [t10s3[subs,:]], [t12s3[subs,:]], [t14s3[subs,:]], [t16s3[subs,:]],
       [t18s3[subs,:]], [t20s3[subs,:]]]

# grab alkalinity for each experiment
ts3_alk = np.zeros((10))
ts3_glass = np.zeros((10))
ts3_kaolinite = np.zeros((10))
ts3_stilbite = np.zeros((10))
ts3_saponite = np.zeros((10))
ts3_albite = np.zeros((10))
ts3_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(ts3[i])
    ts3_alk[i] = (2.0/3.0)*bit[0,-2,3]-(2.0/3.0)*.002407 # steady state alk
    #ts1_alk[i] = 2.0*sum(bit[0,:,3]) # total alk produced for ocean
    ts3_glass[i] = bit[0,-2,56]
    ts3_kaolinite[i] = bit[0,-3,17]
    ts3_stilbite[i] = bit[0,-3,13]
    ts3_saponite[i] = bit[0,-3,21]
    ts3_albite[i] = bit[0,-3,19]
    ts3_celadonite[i] = bit[0,-3,23]

#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 1.57e10     #
#############################
substep = 1.57e9

# load temps
t02s500 = np.loadtxt('t02s500yr.txt')
t04s500 = np.loadtxt('t04s500yr.txt')
t06s500 = np.loadtxt('t06s500yr.txt')
t08s500 = np.loadtxt('t08s500yr.txt')
t10s500 = np.loadtxt('t10s500yr.txt')
t12s500 = np.loadtxt('t12s500yr.txt')
t14s500 = np.loadtxt('t14s500yr.txt')
t16s500 = np.loadtxt('t16s500yr.txt')
t18s500 = np.loadtxt('t18s500yr.txt')
t20s500 = np.loadtxt('t20s500yr.txt')
temps = np.arange(2,22,2)


ts500 = [[t02s500[subs,:]], [t04s500[subs,:]], [t06s500[subs,:]],
         [t08s500[subs,:]], [t10s500[subs,:]], [t12s500[subs,:]],
         [t14s500[subs,:]], [t16s500[subs,:]],
       [t18s500[subs,:]], [t20s500[subs,:]]]

# grab alkalinity for each experiment
ts500_alk = np.zeros((10))
ts500_glass = np.zeros((10))
ts500_kaolinite = np.zeros((10))
ts500_stilbite = np.zeros((10))
ts500_saponite = np.zeros((10))
ts500_albite = np.zeros((10))
ts500_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(ts500[i])
    ts500_alk[i] = 4.0*bit[0,-2,3]-4.0*.002407 # steady state alk
    #ts250_alk[i] = 8.0*sum(bit[0,:,3]) # total alk produced for ocean
    ts500_glass[i] = bit[0,-2,56]
    ts500_kaolinite[i] = bit[0,-3,17]
    ts500_stilbite[i] = bit[0,-3,13]
    ts500_saponite[i] = bit[0,-3,21]
    ts500_albite[i] = bit[0,-3,19]
    ts500_celadonite[i] = bit[0,-3,23]


#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 0.87e10     #
#############################
substep = .87e9

# load temps
t02s250 = np.loadtxt('t02s250yr.txt')
t04s250 = np.loadtxt('t04s250yr.txt')
t06s250 = np.loadtxt('t06s250yr.txt')
t08s250 = np.loadtxt('t08s250yr.txt')
t10s250 = np.loadtxt('t10s250yr.txt')
t12s250 = np.loadtxt('t12s250yr.txt')
t14s250 = np.loadtxt('t14s250yr.txt')
t16s250 = np.loadtxt('t16s250yr.txt')
t18s250 = np.loadtxt('t18s250yr.txt')
t20s250 = np.loadtxt('t20s250yr.txt')
temps = np.arange(2,22,2)


ts250 = [[t02s250[subs,:]], [t04s250[subs,:]], [t06s250[subs,:]],
         [t08s250[subs,:]], [t10s250[subs,:]], [t12s250[subs,:]],
         [t14s250[subs,:]], [t16s250[subs,:]],
       [t18s250[subs,:]], [t20s250[subs,:]]]

# grab alkalinity for each experiment
ts250_alk = np.zeros((10))
ts250_glass = np.zeros((10))
ts250_kaolinite = np.zeros((10))
ts250_stilbite = np.zeros((10))
ts250_saponite = np.zeros((10))
ts250_albite = np.zeros((10))
ts250_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(ts250[i])
    ts250_alk[i] = 8.0*bit[0,-2,3]-8.0*.002407 # steady state alk
    #ts250_alk[i] = 8.0*sum(bit[0,:,3]) # total alk produced for ocean
    ts250_glass[i] = bit[0,-2,56]
    ts250_kaolinite[i] = bit[0,-3,17]
    ts250_stilbite[i] = bit[0,-3,13]
    ts250_saponite[i] = bit[0,-3,21]
    ts250_albite[i] = bit[0,-3,19]
    ts250_celadonite[i] = bit[0,-3,23]










#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 3.14e09     #
#############################
substep = 3.14e8

# load temps
t02s100 = np.loadtxt('t02s100yr.txt')
t04s100 = np.loadtxt('t04s100yr.txt')
t06s100 = np.loadtxt('t06s100yr.txt')
t08s100 = np.loadtxt('t08s100yr.txt')
t10s100 = np.loadtxt('t10s100yr.txt')
t12s100 = np.loadtxt('t12s100yr.txt')
t14s100 = np.loadtxt('t14s100yr.txt')
t16s100 = np.loadtxt('t16s100yr.txt')
t18s100 = np.loadtxt('t18s100yr.txt')
t20s100 = np.loadtxt('t20s100yr.txt')
temps = np.arange(2,22,2)

ts100 = [[t02s100[subs,:]], [t04s100[subs,:]], [t06s100[subs,:]],
         [t08s100[subs,:]], [t10s100[subs,:]], [t12s100[subs,:]],
         [t14s100[subs,:]], [t16s100[subs,:]],
       [t18s100[subs,:]], [t20s100[subs,:]]]

# grab alkalinity for each experiment
ts100_alk = np.zeros((10))
ts100_glass = np.zeros((10))
ts100_kaolinite = np.zeros((10))
ts100_stilbite = np.zeros((10))
ts100_saponite = np.zeros((10))
ts100_albite = np.zeros((10))
ts100_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(ts100[i])
    ts100_alk[i] = 20.0*bit[0,-2,3]-20.0*.002407 # steady state alk
    #ts250_alk[i] = 8.0*sum(bit[0,:,3]) # total alk produced for ocean
    ts100_glass[i] = bit[0,-2,56]
    ts100_kaolinite[i] = bit[0,-3,17]
    ts100_stilbite[i] = bit[0,-3,13]
    ts100_saponite[i] = bit[0,-3,21]
    ts100_albite[i] = bit[0,-3,19]
    ts100_celadonite[i] = bit[0,-3,23]










#--------------------------------------NOPE-----------------------------------#


#############################
# SIMULATIONS WITH MINERALS #
# FLUSHING STEP 3.14e10     #
#############################

# load temps
c080t10s1 = np.loadtxt('c080t10s1kyr.txt')
c100t10s1 = np.loadtxt('c100t10s1kyr.txt')
c120t10s1 = np.loadtxt('c120t10s1kyr.txt')
c140t10s1 = np.loadtxt('c140t10s1kyr.txt')
c160t10s1 = np.loadtxt('c160t10s1kyr.txt')
c180t10s1 = np.loadtxt('c180t10s1kyr.txt')
c200t10s1 = np.loadtxt('c200t10s1kyr.txt')

cs = np.arange(80,220,20)

ct10s1 = [[c080t10s1[subs,:]], [c100t10s1[subs,:]], [c120t10s1[subs,:]],
         [c140t10s1[subs,:]], [c160t10s1[subs,:]], [c180t10s1[subs,:]],
         [c200t10s1[subs,:]]]

# grab alkalinity for each experiment
ct10s1_alk = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_alk[i] = 20.0*bit[0,-2,3]-20.0*bit[0,-2,3] # steady state alk
    #ts250_alk[i] = 8.0*sum(bit[0,:,3]) # total alk produced for ocean

# grab basalt for each experiment
ct10s1_glass = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_glass[i] = bit[0,-2,56]

# grab kaolinite for each experiment
ct10s1_kaolinite = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_kaolinite[i] = bit[0,-2,16]

# grab stilbite for each experiment
ct10s1_stilbite = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_stilbite[i] = bit[0,-2,12]

# grab saponite for each experiment
ct10s1_saponite = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_saponite[i] = bit[0,-2,20]

# grab albite for each experiment
ct10s1_albite = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_albite[i] = bit[0,-2,18]

# grab celadonite for each experiment
ct10s1_celadonite = np.zeros((7))
for i in range(len(cs)):
    bit = np.asarray(ct10s1[i])
    ct10s1_celadonite[i] = bit[0,-2,22]

    
################################
# SIMULATIONS WITH NO MINERALS #
# FLUSHING STEP 6.28e10        #
################################

# load temps
t02n = np.loadtxt('t02n.txt')
t04n = np.loadtxt('t04n.txt')
t06n = np.loadtxt('t06n.txt')
t08n = np.loadtxt('t08n.txt')
t10n = np.loadtxt('t10n.txt')
t12n = np.loadtxt('t12n.txt')
t14n = np.loadtxt('t14n.txt')
t16n = np.loadtxt('t16n.txt')
t18n = np.loadtxt('t18n.txt')
t20n = np.loadtxt('t20n.txt')

tn = [[t02n[subs,:]], [t04n[subs,:]], [t06n[subs,:]], [t08n[subs,:]],
     [t10n[subs,:]], [t12n[subs,:]], [t14n[subs,:]], [t16n[subs,:]],
     [t18n[subs,:]], [t20n[subs,:]]]

# grab alkalinity for each experiment
tn_alk = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(tn[i])
    #temps_alk[i] = bit[0,-2,3]
    tn_alk[i] = bit[0,-50,3]-.002407 # total alk produced for ocean

################################
# SIMULATIONS WITH NO MINERALS #
# FLUSHING STEP 3.14e10        #
################################

# load temps
t02n1 = np.loadtxt('t02n1kyr.txt')
t04n1 = np.loadtxt('t04n1kyr.txt')
t06n1 = np.loadtxt('t06n1kyr.txt')
t08n1 = np.loadtxt('t08n1kyr.txt')
t10n1 = np.loadtxt('t10n1kyr.txt')
t12n1 = np.loadtxt('t12n1kyr.txt')
t14n1 = np.loadtxt('t14n1kyr.txt')
t16n1 = np.loadtxt('t16n1kyr.txt')
t18n1 = np.loadtxt('t18n1kyr.txt')
t20n1 = np.loadtxt('t20n1kyr.txt')

tn1 = [[t02n1[subs,:]], [t04n1[subs,:]], [t06n1[subs,:]], [t08n1[subs,:]],
     [t10n1[subs,:]], [t12n1[subs,:]], [t14n1[subs,:]], [t16n1[subs,:]],
     [t18n1[subs,:]], [t20n1[subs,:]]]

# grab alkalinity for each experiment
tn1_alk = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(tn1[i])
    #temps_alk[i] = bit[0,-2,3]
    tn1_alk[i] = 2.0*bit[0,-50,3]-2.0*.002407 # total alk produced for ocean
    


################################
# SIMULATIONS WITH NO MINERALS #
# FLUSHING STEP 1.57e10        #
################################

# load temps
t02n500 = np.loadtxt('t02n500yr.txt')
t04n500 = np.loadtxt('t04n500yr.txt')
t06n500 = np.loadtxt('t06n500yr.txt')
t08n500 = np.loadtxt('t08n500yr.txt')
t10n500 = np.loadtxt('t10n500yr.txt')
t12n500 = np.loadtxt('t12n500yr.txt')
t14n500 = np.loadtxt('t14n500yr.txt')
t16n500 = np.loadtxt('t16n500yr.txt')
t18n500 = np.loadtxt('t18n500yr.txt')
t20n500 = np.loadtxt('t20n500yr.txt')

tn500 = [[t02n500[subs,:]], [t04n500[subs,:]], [t06n500[subs,:]],
         [t08n500[subs,:]], [t10n500[subs,:]], [t12n500[subs,:]],
         [t14n500[subs,:]], [t16n500[subs,:]],
     [t18n500[subs,:]], [t20n500[subs,:]]]

# grab alkalinity for each experiment
tn500_alk = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(tn500[i])
    #temps_alk[i] = bit[0,-2,3]
    tn500_alk[i] = 4.0*bit[0,-50,3]-4.0*.002407 # total alk produced for ocean




################################
# SIMULATIONS WITH NO MINERALS #
# FLUSHING STEP .87e10        #
################################

# load temps
t02n250 = np.loadtxt('t02n250yr.txt')
t04n250 = np.loadtxt('t04n250yr.txt')
t06n250 = np.loadtxt('t06n250yr.txt')
t08n250 = np.loadtxt('t08n250yr.txt')
t10n250 = np.loadtxt('t10n250yr.txt')
t12n250 = np.loadtxt('t12n250yr.txt')
t14n250 = np.loadtxt('t14n250yr.txt')
t16n250 = np.loadtxt('t16n250yr.txt')
t18n250 = np.loadtxt('t18n250yr.txt')
t20n250 = np.loadtxt('t20n250yr.txt')

tn250 = [[t02n250[subs,:]], [t04n250[subs,:]], [t06n250[subs,:]],
         [t08n250[subs,:]], [t10n250[subs,:]], [t12n250[subs,:]],
         [t14n250[subs,:]], [t16n250[subs,:]],
     [t18n250[subs,:]], [t20n250[subs,:]]]

# grab alkalinity for each experiment
tn250_alk = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(tn250[i])
    #temps_alk[i] = bit[0,-2,3]
    tn250_alk[i] = 8.0*bit[0,-50,3]-8.0*.002407 # total alk produced for ocean

#--------------------------------------NOPE-----------------------------------#



    



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
