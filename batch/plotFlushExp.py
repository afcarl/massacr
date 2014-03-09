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
           t_dH_clay, t_dH_diss, t_ph, t_al, t_na, t_si]
    
    return out

    


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 350 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c475pw.txt')
t04 = np.loadtxt('r90t04c475pw.txt')
t06 = np.loadtxt('r90t06c475pw.txt')
t08 = np.loadtxt('r90t08c475pw.txt')
t10 = np.loadtxt('r90t10c475pw.txt')
t12 = np.loadtxt('r90t12c475pw.txt')
t14 = np.loadtxt('r90t14c475pw.txt')
t16 = np.loadtxt('r90t16c475pw.txt')
t18 = np.loadtxt('r90t18c475pw.txt')
t20 = np.loadtxt('r90t20c475pw.txt')
t22 = np.loadtxt('r90t22c475pw.txt')
t24 = np.loadtxt('r90t24c475pw.txt')
t26 = np.loadtxt('r90t26c475pw.txt')
t28 = np.loadtxt('r90t28c475pw.txt')
t30 = np.loadtxt('r90t30c475pw.txt')
t32 = np.loadtxt('r90t32c475pw.txt')
t34 = np.loadtxt('r90t34c475pw.txt')
t36 = np.loadtxt('r90t36c475pw.txt')
t38 = np.loadtxt('r90t38c475pw.txt')
t40 = np.loadtxt('r90t40c475pw.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c475_alk = out[0]
t90c475_alkflux = out[1]
t90c475_alkflux0 = out[2]
t90c475_glass = out[3]
t90c475_water = out[4]
t90c475_HCO3 = out[5]
t90c475_CO3 = out[6]
t90c475_kaolinite = out[7]
t90c475_stilbite = out[8]
t90c475_saponite = out[9]
t90c475_albite = out[10]
t90c475_celadonite = out[11]
t90c475_quartz = out[12]
t90c475_dH_clay = out[13]
t90c475_dH_diss = out[14]
t90c475_ph = out[15]
t90c475_al = out[16]
t90c475_na = out[17]
t90c475_si = out[18]





#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 350 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c350pw.txt')
t04 = np.loadtxt('r90t04c350pw.txt')
t06 = np.loadtxt('r90t06c350pw.txt')
t08 = np.loadtxt('r90t08c350pw.txt')
t10 = np.loadtxt('r90t10c350pw.txt')
t12 = np.loadtxt('r90t12c350pw.txt')
t14 = np.loadtxt('r90t14c350pw.txt')
t16 = np.loadtxt('r90t16c350pw.txt')
t18 = np.loadtxt('r90t18c350pw.txt')
t20 = np.loadtxt('r90t20c350pw.txt')
t22 = np.loadtxt('r90t22c350pw.txt')
t24 = np.loadtxt('r90t24c350pw.txt')
t26 = np.loadtxt('r90t26c350pw.txt')
t28 = np.loadtxt('r90t28c350pw.txt')
t30 = np.loadtxt('r90t30c350pw.txt')
t32 = np.loadtxt('r90t32c350pw.txt')
t34 = np.loadtxt('r90t34c350pw.txt')
t36 = np.loadtxt('r90t36c350pw.txt')
t38 = np.loadtxt('r90t38c350pw.txt')
t40 = np.loadtxt('r90t40c350pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 225 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c225pw.txt')
t04 = np.loadtxt('r90t04c225pw.txt')
t06 = np.loadtxt('r90t06c225pw.txt')
t08 = np.loadtxt('r90t08c225pw.txt')
t10 = np.loadtxt('r90t10c225pw.txt')
t12 = np.loadtxt('r90t12c225pw.txt')
t14 = np.loadtxt('r90t14c225pw.txt')
t16 = np.loadtxt('r90t16c225pw.txt')
t18 = np.loadtxt('r90t18c225pw.txt')
t20 = np.loadtxt('r90t20c225pw.txt')
t22 = np.loadtxt('r90t22c225pw.txt')
t24 = np.loadtxt('r90t24c225pw.txt')
t26 = np.loadtxt('r90t26c225pw.txt')
t28 = np.loadtxt('r90t28c225pw.txt')
t30 = np.loadtxt('r90t30c225pw.txt')
t32 = np.loadtxt('r90t32c225pw.txt')
t34 = np.loadtxt('r90t34c225pw.txt')
t36 = np.loadtxt('r90t36c225pw.txt')
t38 = np.loadtxt('r90t38c225pw.txt')
t40 = np.loadtxt('r90t40c225pw.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

t90c225_alk = out[0]
t90c225_alkflux = out[1]
t90c225_alkflux0 = out[2]
t90c225_glass = out[3]
t90c225_water = out[4]
t90c225_HCO3 = out[5]
t90c475_CO3 = out[6]
t90c225_kaolinite = out[7]
t90c225_stilbite = out[8]
t90c225_saponite = out[9]
t90c225_albite = out[10]
t90c225_celadonite = out[11]
t90c225_quartz = out[12]
t90c225_dH_clay = out[13]
t90c225_dH_diss = out[14]
t90c225_ph = out[15]
t90c225_al = out[16]
t90c225_na = out[17]
t90c225_si = out[18]


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 550 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c550pw.txt')
t04 = np.loadtxt('r90t04c550pw.txt')
t06 = np.loadtxt('r90t06c550pw.txt')
t08 = np.loadtxt('r90t08c550pw.txt')
t10 = np.loadtxt('r90t10c550pw.txt')
t12 = np.loadtxt('r90t12c550pw.txt')
t14 = np.loadtxt('r90t14c550pw.txt')
t16 = np.loadtxt('r90t16c550pw.txt')
t18 = np.loadtxt('r90t18c550pw.txt')
t20 = np.loadtxt('r90t20c550pw.txt')
t22 = np.loadtxt('r90t22c550pw.txt')
t24 = np.loadtxt('r90t24c550pw.txt')
t26 = np.loadtxt('r90t26c550pw.txt')
t28 = np.loadtxt('r90t28c550pw.txt')
t30 = np.loadtxt('r90t30c550pw.txt')
t32 = np.loadtxt('r90t32c550pw.txt')
t34 = np.loadtxt('r90t34c550pw.txt')
t36 = np.loadtxt('r90t36c550pw.txt')
t38 = np.loadtxt('r90t38c550pw.txt')
t40 = np.loadtxt('r90t40c550pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 500 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c500pw.txt')
t04 = np.loadtxt('r90t04c500pw.txt')
t06 = np.loadtxt('r90t06c500pw.txt')
t08 = np.loadtxt('r90t08c500pw.txt')
t10 = np.loadtxt('r90t10c500pw.txt')
t12 = np.loadtxt('r90t12c500pw.txt')
t14 = np.loadtxt('r90t14c500pw.txt')
t16 = np.loadtxt('r90t16c500pw.txt')
t18 = np.loadtxt('r90t18c500pw.txt')
t20 = np.loadtxt('r90t20c500pw.txt')
t22 = np.loadtxt('r90t22c500pw.txt')
t24 = np.loadtxt('r90t24c500pw.txt')
t26 = np.loadtxt('r90t26c500pw.txt')
t28 = np.loadtxt('r90t28c500pw.txt')
t30 = np.loadtxt('r90t30c500pw.txt')
t32 = np.loadtxt('r90t32c500pw.txt')
t34 = np.loadtxt('r90t34c500pw.txt')
t36 = np.loadtxt('r90t36c500pw.txt')
t38 = np.loadtxt('r90t38c500pw.txt')
t40 = np.loadtxt('r90t40c500pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 450 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c450pw.txt')
t04 = np.loadtxt('r90t04c450pw.txt')
t06 = np.loadtxt('r90t06c450pw.txt')
t08 = np.loadtxt('r90t08c450pw.txt')
t10 = np.loadtxt('r90t10c450pw.txt')
t12 = np.loadtxt('r90t12c450pw.txt')
t14 = np.loadtxt('r90t14c450pw.txt')
t16 = np.loadtxt('r90t16c450pw.txt')
t18 = np.loadtxt('r90t18c450pw.txt')
t20 = np.loadtxt('r90t20c450pw.txt')
t22 = np.loadtxt('r90t22c450pw.txt')
t24 = np.loadtxt('r90t24c450pw.txt')
t26 = np.loadtxt('r90t26c450pw.txt')
t28 = np.loadtxt('r90t28c450pw.txt')
t30 = np.loadtxt('r90t30c450pw.txt')
t32 = np.loadtxt('r90t32c450pw.txt')
t34 = np.loadtxt('r90t34c450pw.txt')
t36 = np.loadtxt('r90t36c450pw.txt')
t38 = np.loadtxt('r90t38c450pw.txt')
t40 = np.loadtxt('r90t40c450pw.txt')


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



#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 400 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c400pw.txt')
t04 = np.loadtxt('r90t04c400pw.txt')
t06 = np.loadtxt('r90t06c400pw.txt')
t08 = np.loadtxt('r90t08c400pw.txt')
t10 = np.loadtxt('r90t10c400pw.txt')
t12 = np.loadtxt('r90t12c400pw.txt')
t14 = np.loadtxt('r90t14c400pw.txt')
t16 = np.loadtxt('r90t16c400pw.txt')
t18 = np.loadtxt('r90t18c400pw.txt')
t20 = np.loadtxt('r90t20c400pw.txt')
t22 = np.loadtxt('r90t22c400pw.txt')
t24 = np.loadtxt('r90t24c400pw.txt')
t26 = np.loadtxt('r90t26c400pw.txt')
t28 = np.loadtxt('r90t28c400pw.txt')
t30 = np.loadtxt('r90t30c400pw.txt')
t32 = np.loadtxt('r90t32c400pw.txt')
t34 = np.loadtxt('r90t34c400pw.txt')
t36 = np.loadtxt('r90t36c400pw.txt')
t38 = np.loadtxt('r90t38c400pw.txt')
t40 = np.loadtxt('r90t40c400pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 300 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c300pw.txt')
t04 = np.loadtxt('r90t04c300pw.txt')
t06 = np.loadtxt('r90t06c300pw.txt')
t08 = np.loadtxt('r90t08c300pw.txt')
t10 = np.loadtxt('r90t10c300pw.txt')
t12 = np.loadtxt('r90t12c300pw.txt')
t14 = np.loadtxt('r90t14c300pw.txt')
t16 = np.loadtxt('r90t16c300pw.txt')
t18 = np.loadtxt('r90t18c300pw.txt')
t20 = np.loadtxt('r90t20c300pw.txt')
t22 = np.loadtxt('r90t22c300pw.txt')
t24 = np.loadtxt('r90t24c300pw.txt')
t26 = np.loadtxt('r90t26c300pw.txt')
t28 = np.loadtxt('r90t28c300pw.txt')
t30 = np.loadtxt('r90t30c300pw.txt')
t32 = np.loadtxt('r90t32c300pw.txt')
t34 = np.loadtxt('r90t34c300pw.txt')
t36 = np.loadtxt('r90t36c300pw.txt')
t38 = np.loadtxt('r90t38c300pw.txt')
t40 = np.loadtxt('r90t40c300pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 250 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c250pw.txt')
t04 = np.loadtxt('r90t04c250pw.txt')
t06 = np.loadtxt('r90t06c250pw.txt')
t08 = np.loadtxt('r90t08c250pw.txt')
t10 = np.loadtxt('r90t10c250pw.txt')
t12 = np.loadtxt('r90t12c250pw.txt')
t14 = np.loadtxt('r90t14c250pw.txt')
t16 = np.loadtxt('r90t16c250pw.txt')
t18 = np.loadtxt('r90t18c250pw.txt')
t20 = np.loadtxt('r90t20c250pw.txt')
t22 = np.loadtxt('r90t22c250pw.txt')
t24 = np.loadtxt('r90t24c250pw.txt')
t26 = np.loadtxt('r90t26c250pw.txt')
t28 = np.loadtxt('r90t28c250pw.txt')
t30 = np.loadtxt('r90t30c250pw.txt')
t32 = np.loadtxt('r90t32c250pw.txt')
t34 = np.loadtxt('r90t34c250pw.txt')
t36 = np.loadtxt('r90t36c250pw.txt')
t38 = np.loadtxt('r90t38c250pw.txt')
t40 = np.loadtxt('r90t40c250pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 200 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c200pw.txt')
t04 = np.loadtxt('r90t04c200pw.txt')
t06 = np.loadtxt('r90t06c200pw.txt')
t08 = np.loadtxt('r90t08c200pw.txt')
t10 = np.loadtxt('r90t10c200pw.txt')
t12 = np.loadtxt('r90t12c200pw.txt')
t14 = np.loadtxt('r90t14c200pw.txt')
t16 = np.loadtxt('r90t16c200pw.txt')
t18 = np.loadtxt('r90t18c200pw.txt')
t20 = np.loadtxt('r90t20c200pw.txt')
t22 = np.loadtxt('r90t22c200pw.txt')
t24 = np.loadtxt('r90t24c200pw.txt')
t26 = np.loadtxt('r90t26c200pw.txt')
t28 = np.loadtxt('r90t28c200pw.txt')
t30 = np.loadtxt('r90t30c200pw.txt')
t32 = np.loadtxt('r90t32c200pw.txt')
t34 = np.loadtxt('r90t34c200pw.txt')
t36 = np.loadtxt('r90t36c200pw.txt')
t38 = np.loadtxt('r90t38c200pw.txt')
t40 = np.loadtxt('r90t40c200pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 150 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c150pw.txt')
t04 = np.loadtxt('r90t04c150pw.txt')
t06 = np.loadtxt('r90t06c150pw.txt')
t08 = np.loadtxt('r90t08c150pw.txt')
t10 = np.loadtxt('r90t10c150pw.txt')
t12 = np.loadtxt('r90t12c150pw.txt')
t14 = np.loadtxt('r90t14c150pw.txt')
t16 = np.loadtxt('r90t16c150pw.txt')
t18 = np.loadtxt('r90t18c150pw.txt')
t20 = np.loadtxt('r90t20c150pw.txt')
t22 = np.loadtxt('r90t22c150pw.txt')
t24 = np.loadtxt('r90t24c150pw.txt')
t26 = np.loadtxt('r90t26c150pw.txt')
t28 = np.loadtxt('r90t28c150pw.txt')
t30 = np.loadtxt('r90t30c150pw.txt')
t32 = np.loadtxt('r90t32c150pw.txt')
t34 = np.loadtxt('r90t34c150pw.txt')
t36 = np.loadtxt('r90t36c150pw.txt')
t38 = np.loadtxt('r90t38c150pw.txt')
t40 = np.loadtxt('r90t40c150pw.txt')


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


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
# A.K.A.                    #
# t_RES = 3.14e11 (10kyr)   #
# pCO2 = 100 ppm            #
#############################

# load temps
t02 = np.loadtxt('r90t02c100pw.txt')
t04 = np.loadtxt('r90t04c100pw.txt')
t06 = np.loadtxt('r90t06c100pw.txt')
t08 = np.loadtxt('r90t08c100pw.txt')
t10 = np.loadtxt('r90t10c100pw.txt')
t12 = np.loadtxt('r90t12c100pw.txt')
t14 = np.loadtxt('r90t14c100pw.txt')
t16 = np.loadtxt('r90t16c100pw.txt')
t18 = np.loadtxt('r90t18c100pw.txt')
t20 = np.loadtxt('r90t20c100pw.txt')
t22 = np.loadtxt('r90t22c100pw.txt')
t24 = np.loadtxt('r90t24c100pw.txt')
t26 = np.loadtxt('r90t26c100pw.txt')
t28 = np.loadtxt('r90t28c100pw.txt')
t30 = np.loadtxt('r90t30c100pw.txt')
t32 = np.loadtxt('r90t32c100pw.txt')
t34 = np.loadtxt('r90t34c100pw.txt')
t36 = np.loadtxt('r90t36c100pw.txt')
t38 = np.loadtxt('r90t38c100pw.txt')
t40 = np.loadtxt('r90t40c100pw.txt')


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



##############################
# PLOT PCOLOR 2D PARAM SPACE #
##############################

fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

pco2 = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 650.0])

grid = np.array([t90c100_alkflux, t90c150_alkflux, t90c200_alkflux,
                 t90c250_alkflux, t90c300_alkflux, t90_alkflux,
                 t90c400_alkflux, t90c450_alkflux, t90c500_alkflux,
                 t90c550_alkflux, t90c475_alkflux])

p = plt.pcolor(pco2, temps, np.transpose(grid),
               cmap=cm.Spectral_r, edgecolors='#444444', linewidth=2)
#p = plt.contourf(pco2, temps, np.transpose(grid),20)

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









########################
# PLOT ALK FLUX (DIFF) #
########################


fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)


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

p = plt.plot(temps,((10.0**-t90c475_ph)**3.0)/t90c475_al,'c--',linewidth=1,label='p475 clay')
p = plt.plot(temps,((10.0**-t90_ph)**3.0)/t90_al,'k--',linewidth=1,label='p350 clay')
p = plt.plot(temps,((10.0**-t90c225_ph)**3.0)/t90c225_al,'b--',linewidth=1,label='p225')
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

#(10.0**-t90c475_ph)

p = plt.plot(np.log10(t90c475_si),np.log10(t90c475_na/(10.0**-t90c475_ph)),'c-',linewidth=1,label='p475 clay')
p = plt.plot(np.log10(t90_si),np.log10(t90_na/(10.0**-t90_ph)),'k-',linewidth=1,label='p350 clay')
p = plt.plot(np.log10(t90c225_si),np.log10(t90c225_na/(10.0**-t90c225_ph)),'b-',linewidth=1,label='p225')
p = plt.plot(np.log10(t90c100_si),np.log10(t90c100_na/(10.0**-t90c100_ph)),'r-',linewidth=1,label='p100')


plt.title('AL3+',fontsize=8)
plt.ylabel('log$_{10}$(a$_{Na+}$ / a$_{H+}$)',fontsize=10)
plt.xlabel('log$_{10}$(H$_4$SiO$_4$) (TEMPERATURE)',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





##############
# PLOT CLAYS #
##############

ax = plt.subplot(2,2,3)


p = plt.plot(temps,t90c475_kaolinite-t90c475_albite,'c*-',linewidth=1,label='albite')
p = plt.plot(temps,t90_kaolinite-t90_albite,'k*-',linewidth=1,label='albite')
p = plt.plot(temps,t90c225_kaolinite-t90c225_albite,'b*-',linewidth=1,label='albite')
p = plt.plot(temps,t90c100_kaolinite-t90c100_albite,'r*-',linewidth=1,label='albite')


plt.title('mineral production rate', fontsize=8)
plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)





plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
