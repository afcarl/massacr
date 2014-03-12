# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,42,2)
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
    for i in range(len(temps)):
        bit = np.asarray(t[i])
        t_alk[i] = np.max(bit[0,:,3])
        gx, gy = np.gradient(bit[0,:,:])
        t_alkflux[i] = gx[m,3]*1000.0
        t_alkflux[i] = bit[0,m,3]

        
        t_alkflux0[i] = bit[0,m,3]*1000.0 
        t_glass[i] = np.max(abs(gx[:,59])) # [mol / kyr]
        t_glass[i] = gx[n,59]

        t_water[i] = bit[0,np.argmax(abs(gx[:,59])),63]
        t_water[i] = bit[0,n,63]
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
        t_na[i] = bit[0,m,6]
        t_si[i] = bit[0,m,10]

        t_kaolinite0[i] = bit[0,n,19]
        t_albite0[i] = bit[0,n,21]

        t_mg[i] = bit[0,m,5]
        t_montna[i] = bit[0,m,31]

    ##    t_kaolinite[i] = np.max(abs(gx[:,19]))
    ##    t_stilbite[i] = np.max(abs(gx[:,15]))
    ##    t90c475_saponite[i] = np.max(abs(gx[:,23]))
    ##    t_albite[i] = np.max(abs(gx[:,21]))
    ##    t_celadonite[i] = np.max(abs(gx[:,25]))
    ##    t_quartz[i] = np.max(abs(gx[:,47]))


    # H+ concentration change from minerals
    t_dH_clay = ((t_kaolinite*6.0 + t_stilbite*8.72 + t_saponite*7.32 + \
              t_albite*4.0 + t_celadonite*6.0) / t_water) 
    # H+ consumption by basalt dissolution
    t_dH_diss = -t_glass * .5 / t_water

    out = [t_alk, t_alkflux, t_alkflux0, t_glass, t_water, t_HCO3, t_CO3,
           t_kaolinite, t_stilbite, t_saponite, t_albite, t_celadonite, t_quartz,
           t_dH_clay, t_dH_diss, t_ph, t_al, t_na, t_si, t_kaolinite0, t_albite0,
           t_mg, t_montna]
    
    return out




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 350 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c350pwa.txt')
t04 = np.loadtxt('r90t04c350pwa.txt')
t06 = np.loadtxt('r90t06c350pwa.txt')
t08 = np.loadtxt('r90t08c350pwa.txt')
t10 = np.loadtxt('r90t10c350pwa.txt')
t12 = np.loadtxt('r90t12c350pwa.txt')
t14 = np.loadtxt('r90t14c350pwa.txt')
t16 = np.loadtxt('r90t16c350pwa.txt')
t18 = np.loadtxt('r90t18c350pwa.txt')
t20 = np.loadtxt('r90t20c350pwa.txt')
t22 = np.loadtxt('r90t22c350pwa.txt')
t24 = np.loadtxt('r90t24c350pwa.txt')
t26 = np.loadtxt('r90t26c350pwa.txt')
t28 = np.loadtxt('r90t28c350pwa.txt')
t30 = np.loadtxt('r90t30c350pwa.txt')
t32 = np.loadtxt('r90t32c350pwa.txt')
t34 = np.loadtxt('r90t34c350pwa.txt')
t36 = np.loadtxt('r90t36c350pwa.txt')
t38 = np.loadtxt('r90t38c350pwa.txt')
t40 = np.loadtxt('r90t40c350pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90_alk = out[0]
t90_alkflux = out[1]
t90_alkflux0 = out[2]
t90_glass = out[3]
t90_water = out[4]
t90_HCO3 = out[5]
t90_CO3 = out[6]
t90_kaolinite = out[7]
t90_stilbite = out[8]
t90_saponite = out[9]
t90_albite = out[10]
t90_celadonite = out[11]
t90_quartz = out[12]
t90_dH_clay = out[13]
t90_dH_diss = out[14]
t90_ph = out[15]
t90_al = out[16]
t90_na = out[17]
t90_si = out[18]
t90_kaolinite0 = out[19]
t90_albite0 = out[20]
t90_mg = out[21]
t90_montna = out[22]


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 550 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c550pwa.txt')
t04 = np.loadtxt('r90t04c550pwa.txt')
t06 = np.loadtxt('r90t06c550pwa.txt')
t08 = np.loadtxt('r90t08c550pwa.txt')
t10 = np.loadtxt('r90t10c550pwa.txt')
t12 = np.loadtxt('r90t12c550pwa.txt')
t14 = np.loadtxt('r90t14c550pwa.txt')
t16 = np.loadtxt('r90t16c550pwa.txt')
t18 = np.loadtxt('r90t18c550pwa.txt')
t20 = np.loadtxt('r90t20c550pwa.txt')
t22 = np.loadtxt('r90t22c550pwa.txt')
t24 = np.loadtxt('r90t24c550pwa.txt')
t26 = np.loadtxt('r90t26c550pwa.txt')
t28 = np.loadtxt('r90t28c550pwa.txt')
t30 = np.loadtxt('r90t30c550pwa.txt')
t32 = np.loadtxt('r90t32c550pwa.txt')
t34 = np.loadtxt('r90t34c550pwa.txt')
t36 = np.loadtxt('r90t36c550pwa.txt')
t38 = np.loadtxt('r90t38c550pwa.txt')
t40 = np.loadtxt('r90t40c550pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c550_alk = out[0]
t90c550_alkflux = out[1]
t90c550_alkflux0 = out[2]
t90c550_glass = out[3]
t90c550_water = out[4]
t90c550_HCO3 = out[5]
t90c550_CO3 = out[6]
t90c550_kaolinite = out[7]
t90c550_stilbite = out[8]
t90c550_saponite = out[9]
t90c550_albite = out[10]
t90c550_celadonite = out[11]
t90c550_quartz = out[12]
t90c550_dH_clay = out[13]
t90c550_dH_diss = out[14]
t90c550_ph = out[15]
t90c550_al = out[16]
t90c550_na = out[17]
t90c550_si = out[18]
t90c550_kaolinite0 = out[19]
t90c550_albite0 = out[20]
t90c550_mg = out[21]
t90c550_montna = out[22]



#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 500 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c500pwa.txt')
t04 = np.loadtxt('r90t04c500pwa.txt')
t06 = np.loadtxt('r90t06c500pwa.txt')
t08 = np.loadtxt('r90t08c500pwa.txt')
t10 = np.loadtxt('r90t10c500pwa.txt')
t12 = np.loadtxt('r90t12c500pwa.txt')
t14 = np.loadtxt('r90t14c500pwa.txt')
t16 = np.loadtxt('r90t16c500pwa.txt')
t18 = np.loadtxt('r90t18c500pwa.txt')
t20 = np.loadtxt('r90t20c500pwa.txt')
t22 = np.loadtxt('r90t22c500pwa.txt')
t24 = np.loadtxt('r90t24c500pwa.txt')
t26 = np.loadtxt('r90t26c500pwa.txt')
t28 = np.loadtxt('r90t28c500pwa.txt')
t30 = np.loadtxt('r90t30c500pwa.txt')
t32 = np.loadtxt('r90t32c500pwa.txt')
t34 = np.loadtxt('r90t34c500pwa.txt')
t36 = np.loadtxt('r90t36c500pwa.txt')
t38 = np.loadtxt('r90t38c500pwa.txt')
t40 = np.loadtxt('r90t40c500pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c500_alk = out[0]
t90c500_alkflux = out[1]
t90c500_alkflux0 = out[2]
t90c500_glass = out[3]
t90c500_water = out[4]
t90c500_HCO3 = out[5]
t90c500_CO3 = out[6]
t90c500_kaolinite = out[7]
t90c500_stilbite = out[8]
t90c500_saponite = out[9]
t90c500_albite = out[10]
t90c500_celadonite = out[11]
t90c500_quartz = out[12]
t90c500_dH_clay = out[13]
t90c500_dH_diss = out[14]
t90c500_ph = out[15]
t90c500_al = out[16]
t90c500_na = out[17]
t90c500_si = out[18]
t90c500_kaolinite0 = out[19]
t90c500_albite0 = out[20]
t90c500_mg = out[21]
t90c500_montna = out[22]



#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 450 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c450pwa.txt')
t04 = np.loadtxt('r90t04c450pwa.txt')
t06 = np.loadtxt('r90t06c450pwa.txt')
t08 = np.loadtxt('r90t08c450pwa.txt')
t10 = np.loadtxt('r90t10c450pwa.txt')
t12 = np.loadtxt('r90t12c450pwa.txt')
t14 = np.loadtxt('r90t14c450pwa.txt')
t16 = np.loadtxt('r90t16c450pwa.txt')
t18 = np.loadtxt('r90t18c450pwa.txt')
t20 = np.loadtxt('r90t20c450pwa.txt')
t22 = np.loadtxt('r90t22c450pwa.txt')
t24 = np.loadtxt('r90t24c450pwa.txt')
t26 = np.loadtxt('r90t26c450pwa.txt')
t28 = np.loadtxt('r90t28c450pwa.txt')
t30 = np.loadtxt('r90t30c450pwa.txt')
t32 = np.loadtxt('r90t32c450pwa.txt')
t34 = np.loadtxt('r90t34c450pwa.txt')
t36 = np.loadtxt('r90t36c450pwa.txt')
t38 = np.loadtxt('r90t38c450pwa.txt')
t40 = np.loadtxt('r90t40c450pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c450_alk = out[0]
t90c450_alkflux = out[1]
t90c450_alkflux0 = out[2]
t90c450_glass = out[3]
t90c450_water = out[4]
t90c450_HCO3 = out[5]
t90c450_CO3 = out[6]
t90c450_kaolinite = out[7]
t90c450_stilbite = out[8]
t90c450_saponite = out[9]
t90c450_albite = out[10]
t90c450_celadonite = out[11]
t90c450_quartz = out[12]
t90c450_dH_clay = out[13]
t90c450_dH_diss = out[14]
t90c450_ph = out[15]
t90c450_al = out[16]
t90c450_na = out[17]
t90c450_si = out[18]
t90c450_kaolinite0 = out[19]
t90c450_albite0 = out[20]
t90c450_mg = out[21]
t90c450_montna = out[22]




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 400 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c400pwa.txt')
t04 = np.loadtxt('r90t04c400pwa.txt')
t06 = np.loadtxt('r90t06c400pwa.txt')
t08 = np.loadtxt('r90t08c400pwa.txt')
t10 = np.loadtxt('r90t10c400pwa.txt')
t12 = np.loadtxt('r90t12c400pwa.txt')
t14 = np.loadtxt('r90t14c400pwa.txt')
t16 = np.loadtxt('r90t16c400pwa.txt')
t18 = np.loadtxt('r90t18c400pwa.txt')
t20 = np.loadtxt('r90t20c400pwa.txt')
t22 = np.loadtxt('r90t22c400pwa.txt')
t24 = np.loadtxt('r90t24c400pwa.txt')
t26 = np.loadtxt('r90t26c400pwa.txt')
t28 = np.loadtxt('r90t28c400pwa.txt')
t30 = np.loadtxt('r90t30c400pwa.txt')
t32 = np.loadtxt('r90t32c400pwa.txt')
t34 = np.loadtxt('r90t34c400pwa.txt')
t36 = np.loadtxt('r90t36c400pwa.txt')
t38 = np.loadtxt('r90t38c400pwa.txt')
t40 = np.loadtxt('r90t40c400pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c400_alk = out[0]
t90c400_alkflux = out[1]
t90c400_alkflux0 = out[2]
t90c400_glass = out[3]
t90c400_water = out[4]
t90c400_HCO3 = out[5]
t90c400_CO3 = out[6]
t90c400_kaolinite = out[7]
t90c400_stilbite = out[8]
t90c400_saponite = out[9]
t90c400_albite = out[10]
t90c400_celadonite = out[11]
t90c400_quartz = out[12]
t90c400_dH_clay = out[13]
t90c400_dH_diss = out[14]
t90c400_ph = out[15]
t90c400_al = out[16]
t90c400_na = out[17]
t90c400_si = out[18]
t90c400_kaolinite0 = out[19]
t90c400_albite0 = out[20]
t90c400_mg = out[21]
t90c400_montna = out[22]



#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 300 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c300pwa.txt')
t04 = np.loadtxt('r90t04c300pwa.txt')
t06 = np.loadtxt('r90t06c300pwa.txt')
t08 = np.loadtxt('r90t08c300pwa.txt')
t10 = np.loadtxt('r90t10c300pwa.txt')
t12 = np.loadtxt('r90t12c300pwa.txt')
t14 = np.loadtxt('r90t14c300pwa.txt')
t16 = np.loadtxt('r90t16c300pwa.txt')
t18 = np.loadtxt('r90t18c300pwa.txt')
t20 = np.loadtxt('r90t20c300pwa.txt')
t22 = np.loadtxt('r90t22c300pwa.txt')
t24 = np.loadtxt('r90t24c300pwa.txt')
t26 = np.loadtxt('r90t26c300pwa.txt')
t28 = np.loadtxt('r90t28c300pwa.txt')
t30 = np.loadtxt('r90t30c300pwa.txt')
t32 = np.loadtxt('r90t32c300pwa.txt')
t34 = np.loadtxt('r90t34c300pwa.txt')
t36 = np.loadtxt('r90t36c300pwa.txt')
t38 = np.loadtxt('r90t38c300pwa.txt')
t40 = np.loadtxt('r90t40c300pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c300_alk = out[0]
t90c300_alkflux = out[1]
t90c300_alkflux0 = out[2]
t90c300_glass = out[3]
t90c300_water = out[4]
t90c300_HCO3 = out[5]
t90c300_CO3 = out[6]
t90c300_kaolinite = out[7]
t90c300_stilbite = out[8]
t90c300_saponite = out[9]
t90c300_albite = out[10]
t90c300_celadonite = out[11]
t90c300_quartz = out[12]
t90c300_dH_clay = out[13]
t90c300_dH_diss = out[14]
t90c300_ph = out[15]
t90c300_al = out[16]
t90c300_na = out[17]
t90c300_si = out[18]
t90c300_kaolinite0 = out[19]
t90c300_albite0 = out[20]
t90c300_mg = out[21]
t90c300_montna = out[22]




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 250 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c250pwa.txt')
t04 = np.loadtxt('r90t04c250pwa.txt')
t06 = np.loadtxt('r90t06c250pwa.txt')
t08 = np.loadtxt('r90t08c250pwa.txt')
t10 = np.loadtxt('r90t10c250pwa.txt')
t12 = np.loadtxt('r90t12c250pwa.txt')
t14 = np.loadtxt('r90t14c250pwa.txt')
t16 = np.loadtxt('r90t16c250pwa.txt')
t18 = np.loadtxt('r90t18c250pwa.txt')
t20 = np.loadtxt('r90t20c250pwa.txt')
t22 = np.loadtxt('r90t22c250pwa.txt')
t24 = np.loadtxt('r90t24c250pwa.txt')
t26 = np.loadtxt('r90t26c250pwa.txt')
t28 = np.loadtxt('r90t28c250pwa.txt')
t30 = np.loadtxt('r90t30c250pwa.txt')
t32 = np.loadtxt('r90t32c250pwa.txt')
t34 = np.loadtxt('r90t34c250pwa.txt')
t36 = np.loadtxt('r90t36c250pwa.txt')
t38 = np.loadtxt('r90t36c250pwa.txt') #####
t40 = np.loadtxt('r90t40c250pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c250_alk = out[0]
t90c250_alkflux = out[1]
t90c250_alkflux0 = out[2]
t90c250_glass = out[3]
t90c250_water = out[4]
t90c250_HCO3 = out[5]
t90c250_CO3 = out[6]
t90c250_kaolinite = out[7]
t90c250_stilbite = out[8]
t90c250_saponite = out[9]
t90c250_albite = out[10]
t90c250_celadonite = out[11]
t90c250_quartz = out[12]
t90c250_dH_clay = out[13]
t90c250_dH_diss = out[14]
t90c250_ph = out[15]
t90c250_al = out[16]
t90c250_na = out[17]
t90c250_si = out[18]
t90c250_kaolinite0 = out[19]
t90c250_albite0 = out[20]
t90c250_mg = out[21]
t90c250_montna = out[22]





#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 200 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c200pwa.txt')
t04 = np.loadtxt('r90t04c200pwa.txt')
t06 = np.loadtxt('r90t06c200pwa.txt')
t08 = np.loadtxt('r90t08c200pwa.txt')
t10 = np.loadtxt('r90t10c200pwa.txt')
t12 = np.loadtxt('r90t12c200pwa.txt')
t14 = np.loadtxt('r90t14c200pwa.txt')
t16 = np.loadtxt('r90t16c200pwa.txt')
t18 = np.loadtxt('r90t18c200pwa.txt')
t20 = np.loadtxt('r90t20c200pwa.txt')
t22 = np.loadtxt('r90t22c200pwa.txt')
t24 = np.loadtxt('r90t24c200pwa.txt')
t26 = np.loadtxt('r90t26c200pwa.txt')
t28 = np.loadtxt('r90t28c200pwa.txt')
t30 = np.loadtxt('r90t30c200pwa.txt')
t32 = np.loadtxt('r90t32c200pwa.txt')
t34 = np.loadtxt('r90t34c200pwa.txt')
t36 = np.loadtxt('r90t36c200pwa.txt')
t38 = np.loadtxt('r90t38c200pwa.txt')
t40 = np.loadtxt('r90t40c200pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c200_alk = out[0]
t90c200_alkflux = out[1]
t90c200_alkflux0 = out[2]
t90c200_glass = out[3]
t90c200_water = out[4]
t90c200_HCO3 = out[5]
t90c200_CO3 = out[6]
t90c200_kaolinite = out[7]
t90c200_stilbite = out[8]
t90c200_saponite = out[9]
t90c200_albite = out[10]
t90c200_celadonite = out[11]
t90c200_quartz = out[12]
t90c200_dH_clay = out[13]
t90c200_dH_diss = out[14]
t90c200_ph = out[15]
t90c200_al = out[16]
t90c200_na = out[17]
t90c200_si = out[18]
t90c200_kaolinite0 = out[19]
t90c200_albite0 = out[20]
t90c200_mg = out[21]
t90c200_montna = out[22]




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 150 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c150pwa.txt')
t04 = np.loadtxt('r90t04c150pwa.txt')
t06 = np.loadtxt('r90t06c150pwa.txt')
t08 = np.loadtxt('r90t08c150pwa.txt')
t10 = np.loadtxt('r90t10c150pwa.txt')
t12 = np.loadtxt('r90t12c150pwa.txt')
t14 = np.loadtxt('r90t14c150pwa.txt')
t16 = np.loadtxt('r90t16c150pwa.txt')
t18 = np.loadtxt('r90t18c150pwa.txt')
t20 = np.loadtxt('r90t20c150pwa.txt')
t22 = np.loadtxt('r90t22c150pwa.txt')
t24 = np.loadtxt('r90t24c150pwa.txt')
t26 = np.loadtxt('r90t26c150pwa.txt')
t28 = np.loadtxt('r90t28c150pwa.txt')
t30 = np.loadtxt('r90t30c150pwa.txt')
t32 = np.loadtxt('r90t32c150pwa.txt')
t34 = np.loadtxt('r90t34c150pwa.txt')
t36 = np.loadtxt('r90t36c150pwa.txt')
t38 = np.loadtxt('r90t38c150pwa.txt')
t40 = np.loadtxt('r90t40c150pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c150_alk = out[0]
t90c150_alkflux = out[1]
t90c150_alkflux0 = out[2]
t90c150_glass = out[3]
t90c150_water = out[4]
t90c150_HCO3 = out[5]
t90c150_CO3 = out[6]
t90c150_kaolinite = out[7]
t90c150_stilbite = out[8]
t90c150_saponite = out[9]
t90c150_albite = out[10]
t90c150_celadonite = out[11]
t90c150_quartz = out[12]
t90c150_dH_clay = out[13]
t90c150_dH_diss = out[14]
t90c150_ph = out[15]
t90c150_al = out[16]
t90c150_na = out[17]
t90c150_si = out[18]
t90c150_kaolinite0 = out[19]
t90c150_albite0 = out[20]
t90c150_mg = out[21]
t90c150_montna = out[22]




#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 100 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c100pwa.txt')
t04 = np.loadtxt('r90t04c100pwa.txt')
t06 = np.loadtxt('r90t06c100pwa.txt')
t08 = np.loadtxt('r90t08c100pwa.txt')
t10 = np.loadtxt('r90t10c100pwa.txt')
t12 = np.loadtxt('r90t12c100pwa.txt')
t14 = np.loadtxt('r90t14c100pwa.txt')
t16 = np.loadtxt('r90t16c100pwa.txt')
t18 = np.loadtxt('r90t18c100pwa.txt')
t20 = np.loadtxt('r90t20c100pwa.txt')
t22 = np.loadtxt('r90t22c100pwa.txt')
t24 = np.loadtxt('r90t24c100pwa.txt')
t26 = np.loadtxt('r90t26c100pwa.txt')
t28 = np.loadtxt('r90t28c100pwa.txt')
t30 = np.loadtxt('r90t30c100pwa.txt')
t32 = np.loadtxt('r90t32c100pwa.txt')
t34 = np.loadtxt('r90t34c100pwa.txt')
t36 = np.loadtxt('r90t36c100pwa.txt')
t38 = np.loadtxt('r90t38c100pwa.txt')
t40 = np.loadtxt('r90t40c100pwa.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c100_alk = out[0]
t90c100_alkflux = out[1]
t90c100_alkflux0 = out[2]
t90c100_glass = out[3]
t90c100_water = out[4]
t90c100_HCO3 = out[5]
t90c100_CO3 = out[6]
t90c100_kaolinite = out[7]
t90c100_stilbite = out[8]
t90c100_saponite = out[9]
t90c100_albite = out[10]
t90c100_celadonite = out[11]
t90c100_quartz = out[12]
t90c100_dH_clay = out[13]
t90c100_dH_diss = out[14]
t90c100_ph = out[15]
t90c100_al = out[16]
t90c100_na = out[17]
t90c100_si = out[18]
t90c100_kaolinite0 = out[19]
t90c100_albite0 = out[20]
t90c100_mg = out[21]
t90c100_montna = out[22]




##############################
# PLOT PCOLOR 2D PARAM SPACE #
##############################

fig=plt.figure()

ax = plt.subplot(1,1,1)

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

pco2 = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0])

grid = np.array([t90c100_alkflux, t90c150_alkflux, t90c200_alkflux,
                 t90c250_alkflux, t90c300_alkflux, t90_alkflux,
                 t90c400_alkflux, t90c450_alkflux, t90c500_alkflux,
                 t90c550_alkflux])

p = plt.pcolor(pco2, temps, np.transpose(grid),
               cmap=cm.Spectral_r, edgecolors='#444444', linewidth=2)
#p = plt.contourf(pco2, temps, np.transpose(grid),20)

cbar = fig.colorbar(p, orientation='horizontal')

plt.xlabel('pCO2 [ppm]',fontsize=10)
plt.ylabel('T [$^{\circ}$C]',fontsize=10)

plt.xticks(pco2+25.0,pco2)
plt.yticks(temps[::-1]+1.0,temps[::-1])

##plt.xticks(pco2,pco2)
##plt.yticks(temps[::-1],temps[::-1])

plt.xlim([100.0,600.0])
plt.ylim([40.0,2.0])

plt.colorbar()

plt.savefig('pcolor.png')




###############################
# PLOT CONTOUR 2D PARAM SPACE #
###############################

fig=plt.figure()

ax = plt.subplot(1,1,1)

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

pco2 = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0])

grid = np.array([t90c100_alkflux, t90c150_alkflux, t90c200_alkflux,
                 t90c250_alkflux, t90c300_alkflux, t90_alkflux,
                 t90c400_alkflux, t90c450_alkflux, t90c500_alkflux,
                 t90c550_alkflux])

#p = plt.pcolor(pco2, temps, np.transpose(grid),
#               cmap=cm.Spectral_r, edgecolors='#444444', linewidth=2)
p = plt.contourf(pco2, temps, np.transpose(grid),20, cm=cm.rainbow)

plt.xlabel('pCO2 [ppm]',fontsize=10)
plt.ylabel('T [$^{\circ}$C]',fontsize=10)
plt.title('ALKALINITY IN SYSTEM [eq kgw$^{-1}$]')

##plt.xticks(pco2+25.0,pco2)
##plt.yticks(temps[::-1]+1.0,temps[::-1])

plt.xticks(pco2,pco2)
plt.yticks(temps[::-1],temps[::-1])

plt.xlim([100.0,550.0])
plt.ylim([40.0,2.0])

plt.colorbar()

plt.savefig('contour.png')






########################
# PLOT ALK FLUX (DIFF) #
########################


fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)


ax = plt.subplot(1,1,1)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,t90c100_alkflux,'r',linewidth=2,label='100ppm CO$_2$')
p = plt.plot(temps,t90c150_alkflux,'orange',linewidth=2,label='150ppm CO$_2$')
p = plt.plot(temps,t90c200_alkflux,'gold',linewidth=2,label='200ppm CO$_2$')
p = plt.plot(temps,t90c250_alkflux,'yellow',linewidth=2,label='250ppm CO$_2$')
p = plt.plot(temps,t90c300_alkflux,'g',linewidth=2,label='300ppm CO$_2$')
p = plt.plot(temps,t90_alkflux,'c',linewidth=2,label='350ppm CO$_2$')
p = plt.plot(temps,t90c400_alkflux,'b',linewidth=2,label='400ppm CO$_2$')
p = plt.plot(temps,t90c450_alkflux,'purple',linewidth=2,label='450ppm CO$_2$')
p = plt.plot(temps,t90c500_alkflux,'m',linewidth=2,label='500ppm CO$_2$')
p = plt.plot(temps,t90c550_alkflux,'k',linewidth=2,label='550ppm CO$_2$')




#plt.text(6, .023, r"+ ALK TO OCEAN", horizontalalignment='center', fontsize=12)
#plt.text(23, -.02, r"- ALK TO OCEAN", horizontalalignment='center', fontsize=12)


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

ax = plt.subplot(1,1,1)


p = plt.plot(np.log10(t90c100_si),np.log10(t90c100_na/(10.0**-t90c100_ph)),
             'red',linewidth=2,label='100 ppm')

p = plt.plot(np.log10(t90c150_si),np.log10(t90c150_na/(10.0**-t90c150_ph)),
             'orange',linewidth=2,label='150 ppm')

p = plt.plot(np.log10(t90c200_si),np.log10(t90c200_na/(10.0**-t90c200_ph)),
             'gold',linewidth=2,label='200 ppm')

p = plt.plot(np.log10(t90c250_si),np.log10(t90c250_na/(10.0**-t90c250_ph)),
             'yellow',linewidth=2,label='250 ppm')

p = plt.plot(np.log10(t90c300_si),np.log10(t90c300_na/(10.0**-t90c300_ph)),
             'g',linewidth=2,label='300 ppm')

p = plt.plot(np.log10(t90_si),np.log10(t90_na/(10.0**-t90_ph)),
             'c',linewidth=2,label='350 ppm')

p = plt.plot(np.log10(t90c400_si),np.log10(t90c400_na/(10.0**-t90c400_ph)),
             'b',linewidth=2,label='400 ppm')

p = plt.plot(np.log10(t90c450_si),np.log10(t90c450_na/(10.0**-t90c450_ph)),
             'purple',linewidth=2,label='450 ppm')

p = plt.plot(np.log10(t90c500_si),np.log10(t90c500_na/(10.0**-t90c500_ph)),
             'm',linewidth=2,label='500 ppm')

p = plt.plot(np.log10(t90c550_si),np.log10(t90c550_na/(10.0**-t90c550_ph)),
             'k',linewidth=2,label='550 ppm')


plt.title('ALUMINUM MINERAL PHASE DIAGRAM',fontsize=16)
plt.ylabel('log$_{10}$([Na$^+$] / [H$^+$])',fontsize=16)
plt.xlabel('log$_{10}$([H$_2$SiO$_4$]) (TEMPERATURE)',fontsize=16)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':10}, ncol=1)



plt.savefig('flushPhase.png')






fig=plt.figure()

plt.rc('xtick', labelsize=7) 
plt.rc('ytick', labelsize=7)



##################
# PLOT dH COMBOS #
##################


ax = plt.subplot(2,2,2)

p = plt.plot(temps,((10.0**-t90c400_ph)**3.0)/t90c400_al,'c--',linewidth=1,label='p400 clay')
p = plt.plot(temps,((10.0**-t90c300_ph)**3.0)/t90c300_al,'k--',linewidth=1,label='p300 clay')
p = plt.plot(temps,((10.0**-t90c200_ph)**3.0)/t90c200_al,'b--',linewidth=1,label='p200')
p = plt.plot(temps,((10.0**-t90c100_ph)**3.0)/t90c100_al,'r--',linewidth=1,label='p100')



plt.title('pH',fontsize=8)
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

# [5:11]

p = plt.plot(temps[0:20],t90c400_dH_clay[0:20]-t90c400_dH_diss[0:20],
             'c-',linewidth=1,label='p475 both')
p = plt.plot(temps[0:20],t90c300_dH_clay[0:20]-t90c300_dH_diss[0:20],
             'k-',linewidth=1,label='p350')
p = plt.plot(temps[0:20],t90c200_dH_clay[0:20]-t90c200_dH_diss[0:20],
             'b-',linewidth=1,label='p225')
p = plt.plot(temps[0:20],t90c100_dH_clay[0:20]-t90c100_dH_diss[0:20],
             'r-',linewidth=1,label='p100')


plt.title('PRODUCTION - H+ CONSUMPTION',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

#(10.0**-t90c475_ph)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





##############
# PLOT CLAYS #
##############

ax = plt.subplot(2,2,3)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,t90c300_montna,'-',linewidth=1,label='400ppm mont-na')
p = plt.plot(temps,t90c300_albite,'-',linewidth=1,label='400ppm albite')



plt.title('mineral production rate', fontsize=8)
plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
