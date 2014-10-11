# plot prelimGrid.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math


plt.rcParams['axes.color_cycle'] = "#FF3300, #FF9900, #FFCC00, \
                                    #00FF00, #339900, #009966, #0000FF, \
                                    #6600CC, #990099"


plt.rcParams['axes.labelsize'] = 7
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7



temps = np.arange(2,36,2)


def grab(t):
    t_alk = np.zeros((len(temps)))
    t_dalk = np.zeros((len(temps)))
    t_calcite = np.zeros((len(temps)))
    t_dcalcite = np.zeros((len(temps)))
    t_ph = np.zeros((len(temps)))
    t_water = np.zeros((len(temps)))
    t_dglass = np.zeros((len(temps)))
    for i in range(len(temps)):
        bit = np.asarray(t[i])
        t_alk[i] = bit[0,3,-1]
        t_dalk[i] = bit[0,3,-1] - bit[0,3,-2]
        t_calcite[i] = bit[0,45,-1]
        t_dcalcite[i] = bit[0,45,-1] - bit[0,45,-2]
        t_ph[i] = bit[0,1,-1]
        t_water[i] = bit[0,83,-1]
        t_dglass[i] = bit[0,80,-1]
        

    output = [t_alk, t_calcite, t_dcalcite, t_ph, t_water, t_dalk, t_dglass]
    return output


    


#############################
# flush 1% / 10 years       #
# t = 3.14e8 (10 years)     #
#############################

t02 = np.loadtxt('f01_t02.txt')
t04 = np.loadtxt('f01_t04.txt')
t06 = np.loadtxt('f01_t06.txt')
t08 = np.loadtxt('f01_t08.txt')
t10 = np.loadtxt('f01_t10.txt')
t12 = np.loadtxt('f01_t12.txt')
t14 = np.loadtxt('f01_t14.txt')
t16 = np.loadtxt('f01_t16.txt')
t18 = np.loadtxt('f01_t18.txt')
t20 = np.loadtxt('f01_t20.txt')
t22 = np.loadtxt('f01_t22.txt')
t24 = np.loadtxt('f01_t24.txt')
t26 = np.loadtxt('f01_t26.txt')
t28 = np.loadtxt('f01_t28.txt')
t30 = np.loadtxt('f01_t30.txt')
t32 = np.loadtxt('f01_t32.txt')
t34 = np.loadtxt('f01_t34.txt')
t36 = np.loadtxt('f01_t36.txt')
t38 = np.loadtxt('f01_t38.txt')
t40 = np.loadtxt('f01_t40.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

f01_dcalcite = out[2] #2.0*out[2] + out[5]*out[4]
f01_dglass = out[6]





#print dcalcite_230


fig=plt.figure()

ax = plt.subplot(2,2,1)

plt.plot(temps, dcalcite_f01)
plt.plot(temps, dcalcite_f02)
plt.plot(temps, dcalcite_f03)
plt.plot(temps, dcalcite_f04)
plt.plot(temps, dcalcite_f05)

plt.title('calcite formation rate')




ax = plt.subplot(2,2,2)

plt.plot(temps, dglass_f01)
plt.plot(temps, dglass_f02)
plt.plot(temps, dglass_f03)
plt.plot(temps, dglass_f04)
plt.plot(temps, dglass_f05)

plt.title('basalt alteration rate')


ax = plt.subplot(2,2,3)


flush = np.array([.01, .02, .03, .04, .05])


grid = np.array([dcalcite_f01, dcalcite_f02, dcalcite_f03,
                 dcalcite_f04, dcalcite_f05])


p = plt.contourf(flush, temps, np.transpose(grid), 20,
               cmap=cm.Spectral_r)


plt.savefig('prelimGrid.png')

