# plotFlushMin.py
#
# same as plotFlushExp.py for batchControl experiments,
# except with no positive saturations (mo minerals mo problems)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,32,2)
m=-8 # alk/other step
n=-8 # precip step

def grab(t):
    # grab for each experiment
    t_alk = np.zeros((len(temps)))
    t_alkflux = np.zeros((len(temps)))
    t_alkflux0 = np.zeros((len(temps)))
    t_glass = np.zeros((len(temps)))
    t_water = np.zeros((len(temps)))
    t_HCO3 = np.zeros((len(temps)))
    t_CO3 = np.zeros((len(temps)))
    t_kaolinite = np.zeros((len(temps)))
    t_stilbite = np.zeros((len(temps)))
    t_saponite = np.zeros((len(temps)))
    t_albite = np.zeros((len(temps)))
    t_celadonite = np.zeros((len(temps)))
    t_quartz = np.zeros((len(temps)))
    t_ph = np.zeros((len(temps)))
    t_al = np.zeros((len(temps)))
    t_na = np.zeros((len(temps)))
    t_si = np.zeros((len(temps)))
    t_kaolinite0 = np.zeros((len(temps)))
    t_albite0 = np.zeros((len(temps)))
    t_mg = np.zeros((len(temps)))
    t_montna = np.zeros((len(temps)))
    t_dolomite = np.zeros((len(temps)))
    t_calcite = np.zeros((len(temps)))
    t_ca = np.zeros((len(temps)))
    for i in range(len(temps)):
        bit = np.asarray(t[i])
        t_alk[i] = np.max(bit[0,:,3])
        gx, gy = np.gradient(bit[0,:,:])
        
        t_alkflux[i] = bit[0,m,3]# + (2.0*bit[0,m,35]/bit[0,n,83])

        t_alkflux0[i] = bit[0,m,3]*1000.0 
        t_glass[i] = np.max(abs(gx[:,79])) # [mol / kyr]
        t_glass[i] = np.abs(gx[n,79])

        t_water[i] = bit[0,n,83]
        t_HCO3[i] = bit[0,m,13]
        t_CO3[i] = bit[0,m,14]

        t_kaolinite[i] = gx[n,19]
        t_stilbite[i] = gx[n,15]
        t_saponite[i] = gx[n,23]
        t_albite[i] = gx[n,21]
        t_celadonite[i] = gx[n,25]
        t_quartz[i] = gx[n,47]

        t_ph[i] = bit[0,m,2]
        t_al[i] = bit[0,m,12]
        t_na[i] = bit[0,m,7]
        t_si[i] = bit[0,m,10]

        t_kaolinite0[i] = bit[0,n,19]
        t_albite0[i] = bit[0,n,21]

        t_mg[i] = bit[0,m,5]
        t_montna[i] = bit[0,m,31]
        
        t_dolomite[i] = bit[0,m,35]
        t_calcite[i] = gx[m,45]

        t_ca[i] = bit[0,m,5]

    # H+ concentration change from minerals
    t_dH_clay = ((t_kaolinite*6.0 + t_stilbite*8.72 + t_saponite*7.32 + \
              t_albite*4.0 + t_celadonite*6.0) / t_water) 
    # H+ consumption by basalt dissolution
    t_dH_diss = -t_glass * .5 / t_water

    out = [t_alk, t_alkflux, t_alkflux0, t_glass, t_water, t_HCO3, t_CO3,
           t_kaolinite, t_stilbite, t_saponite, t_albite, t_celadonite, t_quartz,
           t_dH_clay, t_dH_diss, t_ph, t_al, t_na, t_si, t_kaolinite0, t_albite0,
           t_mg, t_montna, t_dolomite, t_calcite, t_ca]
    
    return out





#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 2.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c02.txt')
t04 = np.loadtxt('r10kt04c02.txt')
t06 = np.loadtxt('r10kt06c02.txt')
t08 = np.loadtxt('r10kt08c02.txt')
t10 = np.loadtxt('r10kt10c02.txt')
t12 = np.loadtxt('r10kt12c02.txt')
t14 = np.loadtxt('r10kt14c02.txt')
t16 = np.loadtxt('r10kt16c02.txt')
t18 = np.loadtxt('r10kt18c02.txt')
t20 = np.loadtxt('r10kt20c02.txt')
t22 = np.loadtxt('r10kt22c02.txt')
t24 = np.loadtxt('r10kt24c02.txt')
t26 = np.loadtxt('r10kt26c02.txt')
t28 = np.loadtxt('r10kt28c02.txt')
t30 = np.loadtxt('r10kt30c02.txt')
##t32 = np.loadtxt('r10kt32c02.txt')
##t34 = np.loadtxt('r10kt34c02.txt')
##t36 = np.loadtxt('r10kt36c02.txt')
##t38 = np.loadtxt('r10kt38c02.txt')
##t40 = np.loadtxt('r10kt40c02.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]]]
#     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c02_alk = out[0]
c02_alkflux = out[1]
c02_alkflux0 = out[2]
c02_glass = out[3]
c02_water = out[4]
c02_HCO3 = out[5]
c02_CO3 = out[6]
c02_kaolinite = out[7]
c02_stilbite = out[8]
c02_saponite = out[9]
c02_albite = out[10]
c02_celadonite = out[11]
c02_quartz = out[12]
c02_dH_clay = out[13]
c02_dH_diss = out[14]
c02_ph = out[15]
c02_al = out[16]
c02_na = out[17]
c02_si = out[18]
c02_kaolinite0 = out[19]
c02_albite0 = out[20]
c02_mg = out[21]
c02_montna = out[22]
c02_dolomite = out[23]
c02_calcite = out[24]
c02_ca = out[25]





#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 3.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c03.txt')
t04 = np.loadtxt('r10kt04c03.txt')
t06 = np.loadtxt('r10kt06c03.txt')
t08 = np.loadtxt('r10kt08c03.txt')
t10 = np.loadtxt('r10kt10c03.txt')
t12 = np.loadtxt('r10kt12c03.txt')
t14 = np.loadtxt('r10kt14c03.txt')
t16 = np.loadtxt('r10kt16c03.txt')
t18 = np.loadtxt('r10kt18c03.txt')
t20 = np.loadtxt('r10kt20c03.txt')
t22 = np.loadtxt('r10kt22c03.txt')
t24 = np.loadtxt('r10kt24c03.txt')
t26 = np.loadtxt('r10kt26c03.txt')
t28 = np.loadtxt('r10kt28c03.txt')
t30 = np.loadtxt('r10kt30c03.txt')
##t32 = np.loadtxt('r10kt32c03.txt')
##t34 = np.loadtxt('r10kt34c03.txt')
##t36 = np.loadtxt('r10kt36c03.txt')
##t38 = np.loadtxt('r10kt38c03.txt')
##t40 = np.loadtxt('r10kt40c03.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]]]
#     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c03_alk = out[0]
c03_alkflux = out[1]
c03_alkflux0 = out[2]
c03_glass = out[3]
c03_water = out[4]
c03_HCO3 = out[5]
c03_CO3 = out[6]
c03_kaolinite = out[7]
c03_stilbite = out[8]
c03_saponite = out[9]
c03_albite = out[10]
c03_celadonite = out[11]
c03_quartz = out[12]
c03_dH_clay = out[13]
c03_dH_diss = out[14]
c03_ph = out[15]
c03_al = out[16]
c03_na = out[17]
c03_si = out[18]
c03_kaolinite0 = out[19]
c03_albite0 = out[20]
c03_mg = out[21]
c03_montna = out[22]
c03_dolomite = out[23]
c03_calcite = out[24]
c03_ca = out[25]



#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 4.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c04.txt')
t04 = np.loadtxt('r10kt04c04.txt')
t06 = np.loadtxt('r10kt06c04.txt')
t08 = np.loadtxt('r10kt08c04.txt')
t10 = np.loadtxt('r10kt10c04.txt')
t12 = np.loadtxt('r10kt12c04.txt')
t14 = np.loadtxt('r10kt14c04.txt')
t16 = np.loadtxt('r10kt16c04.txt')
t18 = np.loadtxt('r10kt18c04.txt')
t20 = np.loadtxt('r10kt20c04.txt')
t22 = np.loadtxt('r10kt22c04.txt')
t24 = np.loadtxt('r10kt24c04.txt')
t26 = np.loadtxt('r10kt26c04.txt')
t28 = np.loadtxt('r10kt28c04.txt')
t30 = np.loadtxt('r10kt30c04.txt')
##t32 = np.loadtxt('r10kt32c04.txt')
##t34 = np.loadtxt('r10kt34c04.txt')
##t36 = np.loadtxt('r10kt36c04.txt')
##t38 = np.loadtxt('r10kt38c04.txt')
##t40 = np.loadtxt('r10kt40c04.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]]]
#     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c04_alk = out[0]
c04_alkflux = out[1]
c04_alkflux0 = out[2]
c04_glass = out[3]
c04_water = out[4]
c04_HCO3 = out[5]
c04_CO3 = out[6]
c04_kaolinite = out[7]
c04_stilbite = out[8]
c04_saponite = out[9]
c04_albite = out[10]
c04_celadonite = out[11]
c04_quartz = out[12]
c04_dH_clay = out[13]
c04_dH_diss = out[14]
c04_ph = out[15]
c04_al = out[16]
c04_na = out[17]
c04_si = out[18]
c04_kaolinite0 = out[19]
c04_albite0 = out[20]
c04_mg = out[21]
c04_montna = out[22]
c04_dolomite = out[23]
c04_calcite = out[24]
c04_ca = out[25]





##############################
# PLOT PCOLOR 2D PARAM SPACE #
##############################

fig=plt.figure()

ax = plt.subplot(1,1,1)

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

dic = np.array([.001])

#grid = np.array([c02_alkflux])

#p = plt.pcolor(dic, temps, np.transpose(grid),
#               cmap=cm.Spectral_r, edgecolors='#444444', linewidth=2)


plt.xlabel('pCO2 [ppm]',fontsize=10)
plt.ylabel('T [$^{\circ}$C]',fontsize=10)

plt.xticks(dic+.0005,dic)
plt.yticks(temps[::-1]+1.0,temps[::-1])

##plt.xticks(pco2,pco2)
##plt.yticks(temps[::-1],temps[::-1])

plt.xlim([.001,.004])
plt.ylim([40.0,2.0])


plt.savefig('pcolor0.png')









############################
# PLOT CALCITE GROWTH RATE #
############################


fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)


ax = plt.subplot(2,2,1)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,c02_calcite*100.0,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03_calcite*100.0,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04_calcite*100.0,linewidth=2,label='DIC = 4.0 mmol/kgw')


plt.title('CALCITE GROWTH RATE',fontsize=10)
plt.ylabel('[mol / yr ]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



ax = plt.subplot(2,2,2)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,c02_alkflux*100.0,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03_alkflux*100.0,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04_alkflux*100.0,linewidth=2,label='DIC = 4.0 mmol/kgw')

plt.title('ALKALINITY FLUX',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
#plt.legend(handles, labels,loc=3,prop={'size':6}, ncol=1)




ax = plt.subplot(2,2,3)

p = plt.plot(temps,c02_glass*100.0,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03_glass*100.0,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04_glass*100.0,linewidth=2,label='DIC = 4.0 mmol/kgw')

plt.title('Ca2+',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

#handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
#plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



plt.savefig('flushAlk0.png')









fig=plt.figure()

plt.rc('xtick', labelsize=7) 
plt.rc('ytick', labelsize=7)



######################
# PLOT PHASE DIAGRAM #
######################

ax = plt.subplot(1,1,1)



plt.title('',fontsize=16)
plt.ylabel('',fontsize=16)
plt.xlabel('',fontsize=16)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':10}, ncol=1)


plt.savefig('flushPhase0.png')






fig=plt.figure()

plt.rc('xtick', labelsize=7) 
plt.rc('ytick', labelsize=7)


###########
# SUBPLOT #
###########

ax = plt.subplot(2,2,2)


plt.title('',fontsize=8)
plt.ylabel('',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)


###########
# SUBPLOT #
###########

ax = plt.subplot(2,2,4)

# [5:11]



plt.title('',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

#(10.0**-t90c475_ph)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)


###########
# SUBPLOT #
###########

ax = plt.subplot(2,2,3)

plt.title('', fontsize=8)
plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp0.png')
