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
# flush 1% / year           #
#############################

t02 = np.loadtxt('c230_t02.txt')
t04 = np.loadtxt('c230_t04.txt')
t06 = np.loadtxt('c230_t06.txt')
t08 = np.loadtxt('c230_t08.txt')
t10 = np.loadtxt('c230_t10.txt')
t12 = np.loadtxt('c230_t12.txt')
t14 = np.loadtxt('c230_t14.txt')
t16 = np.loadtxt('c230_t16.txt')
t18 = np.loadtxt('c230_t18.txt')
t20 = np.loadtxt('c230_t20.txt')
t22 = np.loadtxt('c230_t22.txt')
t24 = np.loadtxt('c230_t24.txt')
t26 = np.loadtxt('c230_t26.txt')
t28 = np.loadtxt('c230_t28.txt')
t30 = np.loadtxt('c230_t30.txt')
t32 = np.loadtxt('c230_t32.txt')
t34 = np.loadtxt('c230_t34.txt')
t36 = np.loadtxt('c230_t36.txt')
t38 = np.loadtxt('c230_t38.txt')
t40 = np.loadtxt('c230_t40.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f01 = out[2] #2.0*out[2] + out[5]*out[4]
dglass_f01 = out[6]


#############################
# flush 2% / year           #
#############################

##t02 = np.loadtxt('c330_t02.txt')
##t04 = np.loadtxt('c330_t04.txt')
##t06 = np.loadtxt('c330_t06.txt')
##t08 = np.loadtxt('c330_t08.txt')
##t10 = np.loadtxt('c330_t10.txt')
##t12 = np.loadtxt('c330_t12.txt')
##t14 = np.loadtxt('c330_t14.txt')
##t16 = np.loadtxt('c330_t16.txt')
##t18 = np.loadtxt('c330_t18.txt')
##t20 = np.loadtxt('c330_t20.txt')
##t22 = np.loadtxt('c330_t22.txt')
##t24 = np.loadtxt('c330_t24.txt')
##t26 = np.loadtxt('c330_t26.txt')
##t28 = np.loadtxt('c330_t28.txt')
##t30 = np.loadtxt('c330_t30.txt')
##t32 = np.loadtxt('c330_t32.txt')
##t34 = np.loadtxt('c330_t34.txt')
##t36 = np.loadtxt('c330_t36.txt')
##t38 = np.loadtxt('c330_t38.txt')
##t40 = np.loadtxt('c330_t40.txt')

t02 = np.loadtxt('c230_t02_f02.txt')
t04 = np.loadtxt('c230_t04_f02.txt')
t06 = np.loadtxt('c230_t06_f02.txt')
t08 = np.loadtxt('c230_t08_f02.txt')
t10 = np.loadtxt('c230_t10_f02.txt')
t12 = np.loadtxt('c230_t12_f02.txt')
t14 = np.loadtxt('c230_t14_f02.txt')
t16 = np.loadtxt('c230_t16_f02.txt')
t18 = np.loadtxt('c230_t18_f02.txt')
t20 = np.loadtxt('c230_t20_f02.txt')
t22 = np.loadtxt('c230_t22_f02.txt')
t24 = np.loadtxt('c230_t24_f02.txt')
t26 = np.loadtxt('c230_t26_f02.txt')
t28 = np.loadtxt('c230_t28_f02.txt')
t30 = np.loadtxt('c230_t30_f02.txt')
t32 = np.loadtxt('c230_t32_f02.txt')
t34 = np.loadtxt('c230_t34_f02.txt')
t36 = np.loadtxt('c230_t36_f02.txt')
t38 = np.loadtxt('c230_t38_f02.txt')
t40 = np.loadtxt('c230_t40_f02.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f02 = out[2] #out[2] + out[5]*out[4]
dglass_f02 = out[6]


#############################
# flush 3% / year           #
#############################

t02 = np.loadtxt('c230_t02_f03.txt')
t04 = np.loadtxt('c230_t04_f03.txt')
t06 = np.loadtxt('c230_t06_f03.txt')
t08 = np.loadtxt('c230_t08_f03.txt')
t10 = np.loadtxt('c230_t10_f03.txt')
t12 = np.loadtxt('c230_t12_f03.txt')
t14 = np.loadtxt('c230_t14_f03.txt')
t16 = np.loadtxt('c230_t16_f03.txt')
t18 = np.loadtxt('c230_t18_f03.txt')
t20 = np.loadtxt('c230_t20_f03.txt')
t22 = np.loadtxt('c230_t22_f03.txt')
t24 = np.loadtxt('c230_t24_f03.txt')
t26 = np.loadtxt('c230_t26_f03.txt')
t28 = np.loadtxt('c230_t28_f03.txt')
t30 = np.loadtxt('c230_t30_f03.txt')
t32 = np.loadtxt('c230_t32_f03.txt')
t34 = np.loadtxt('c230_t34_f03.txt')
t36 = np.loadtxt('c230_t36_f03.txt')
t38 = np.loadtxt('c230_t38_f03.txt')
t40 = np.loadtxt('c230_t40_f03.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f03 = out[2]# + out[5]*out[4]
dglass_f03 = out[6]



#############################
# flush 4% / year           #
#############################

t02 = np.loadtxt('c230_t02_f04.txt')
t04 = np.loadtxt('c230_t04_f04.txt')
t06 = np.loadtxt('c230_t06_f04.txt')
t08 = np.loadtxt('c230_t08_f04.txt')
t10 = np.loadtxt('c230_t10_f04.txt')
t12 = np.loadtxt('c230_t12_f04.txt')
t14 = np.loadtxt('c230_t14_f04.txt')
t16 = np.loadtxt('c230_t16_f04.txt')
t18 = np.loadtxt('c230_t18_f04.txt')
t20 = np.loadtxt('c230_t20_f04.txt')
t22 = np.loadtxt('c230_t22_f04.txt')
t24 = np.loadtxt('c230_t24_f04.txt')
t26 = np.loadtxt('c230_t26_f04.txt')
t28 = np.loadtxt('c230_t28_f04.txt')
t30 = np.loadtxt('c230_t30_f04.txt')
t32 = np.loadtxt('c230_t32_f04.txt')
t34 = np.loadtxt('c230_t34_f04.txt')
t36 = np.loadtxt('c230_t36_f04.txt')
t38 = np.loadtxt('c230_t38_f04.txt')
t40 = np.loadtxt('c230_t40_f04.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f04 = out[2]# + out[5]*out[4]
dglass_f04 = out[6]






#############################
# flush 5% / year           #
#############################

t02 = np.loadtxt('c230_t02_f05.txt')
t04 = np.loadtxt('c230_t04_f05.txt')
t06 = np.loadtxt('c230_t06_f05.txt')
t08 = np.loadtxt('c230_t08_f05.txt')
t10 = np.loadtxt('c230_t10_f05.txt')
t12 = np.loadtxt('c230_t12_f05.txt')
t14 = np.loadtxt('c230_t14_f05.txt')
t16 = np.loadtxt('c230_t16_f05.txt')
t18 = np.loadtxt('c230_t18_f05.txt')
t20 = np.loadtxt('c230_t20_f05.txt')
t22 = np.loadtxt('c230_t22_f05.txt')
t24 = np.loadtxt('c230_t24_f05.txt')
t26 = np.loadtxt('c230_t26_f05.txt')
t28 = np.loadtxt('c230_t28_f05.txt')
t30 = np.loadtxt('c230_t30_f05.txt')
t32 = np.loadtxt('c230_t32_f05.txt')
t34 = np.loadtxt('c230_t34_f05.txt')
t36 = np.loadtxt('c230_t36_f05.txt')
t38 = np.loadtxt('c230_t38_f05.txt')
t40 = np.loadtxt('c230_t40_f05.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f05 = out[2]# + out[5]*out[4]
dglass_f05 = out[6]




#############################
# flush 0.1% / year           #
#############################

t02 = np.loadtxt('c230_t02_f001.txt')
t04 = np.loadtxt('c230_t04_f001.txt')
t06 = np.loadtxt('c230_t06_f001.txt')
t08 = np.loadtxt('c230_t08_f001.txt')
t10 = np.loadtxt('c230_t10_f001.txt')
t12 = np.loadtxt('c230_t12_f001.txt')
t14 = np.loadtxt('c230_t14_f001.txt')
t16 = np.loadtxt('c230_t16_f001.txt')
t18 = np.loadtxt('c230_t18_f001.txt')
t20 = np.loadtxt('c230_t20_f001.txt')
t22 = np.loadtxt('c230_t22_f001.txt')
t24 = np.loadtxt('c230_t24_f001.txt')
t26 = np.loadtxt('c230_t26_f001.txt')
t28 = np.loadtxt('c230_t28_f001.txt')
t30 = np.loadtxt('c230_t30_f001.txt')
t32 = np.loadtxt('c230_t32_f001.txt')
t34 = np.loadtxt('c230_t34_f001.txt')
t36 = np.loadtxt('c230_t36_f001.txt')
t38 = np.loadtxt('c230_t38_f001.txt')
t40 = np.loadtxt('c230_t40_f001.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f0001 = out[2]# + out[5]*out[4]
dcalcite_f0001[-2] = dcalcite_f0001[-3]
dcalcite_f0001[-1] = dcalcite_f0001[-3]
dglass_f0001 = out[6]




#############################
# flush 0.01% / year           #
#############################

t02 = np.loadtxt('c230_t02_f0001.txt')
t04 = np.loadtxt('c230_t04_f0001.txt')
t06 = np.loadtxt('c230_t06_f0001.txt')
t08 = np.loadtxt('c230_t08_f0001.txt')
t10 = np.loadtxt('c230_t10_f0001.txt')
t12 = np.loadtxt('c230_t12_f0001.txt')
t14 = np.loadtxt('c230_t14_f0001.txt')
t16 = np.loadtxt('c230_t16_f0001.txt')
t18 = np.loadtxt('c230_t18_f0001.txt')
t20 = np.loadtxt('c230_t20_f0001.txt')
t22 = np.loadtxt('c230_t22_f0001.txt')
t24 = np.loadtxt('c230_t24_f0001.txt')
t26 = np.loadtxt('c230_t26_f0001.txt')
t28 = np.loadtxt('c230_t28_f0001.txt')
t30 = np.loadtxt('c230_t30_f0001.txt')
t32 = np.loadtxt('c230_t32_f0001.txt')
t34 = np.loadtxt('c230_t34_f0001.txt')
t36 = np.loadtxt('c230_t36_f0001.txt')
t38 = np.loadtxt('c230_t38_f0001.txt')
t40 = np.loadtxt('c230_t40_f0001.txt')
t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]


out = grab(t)

dcalcite_f001 = out[2]# + out[5]*out[4]
dcalcite_f001[-2] = dcalcite_f001[-3]
dcalcite_f001[-1] = dcalcite_f001[-3]
dglass_f001 = out[6]



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

