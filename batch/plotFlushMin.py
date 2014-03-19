# plotFlushMin.py
#
# same as plotFlushExp.py for batchControl experiments,
# except with no positive saturations (mo minerals mo problems)

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
# DIC = 1.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c01.txt')
t04 = np.loadtxt('r10kt04c01.txt')
t06 = np.loadtxt('r10kt06c01.txt')
t08 = np.loadtxt('r10kt08c01.txt')
t10 = np.loadtxt('r10kt10c01.txt')
t12 = np.loadtxt('r10kt12c01.txt')
t14 = np.loadtxt('r10kt14c01.txt')
t16 = np.loadtxt('r10kt16c01.txt')
t18 = np.loadtxt('r10kt18c01.txt')
t20 = np.loadtxt('r10kt20c01.txt')
t22 = np.loadtxt('r10kt22c01.txt')
t24 = np.loadtxt('r10kt24c01.txt')
t26 = np.loadtxt('r10kt26c01.txt')
t28 = np.loadtxt('r10kt28c01.txt')
t30 = np.loadtxt('r10kt30c01.txt')
t32 = np.loadtxt('r10kt32c01.txt')
t34 = np.loadtxt('r10kt34c01.txt')
t36 = np.loadtxt('r10kt36c01.txt')
t38 = np.loadtxt('r10kt38c01.txt')
t40 = np.loadtxt('r10kt40c01.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c01_alk = out[0]
c01_alkflux = out[1]
c01_alkflux0 = out[2]
c01_glass = out[3]
c01_water = out[4]
c01_HCO3 = out[5]
c01_CO3 = out[6]
c01_kaolinite = out[7]
c01_stilbite = out[8]
c01_saponite = out[9]
c01_albite = out[10]
c01_celadonite = out[11]
c01_quartz = out[12]
c01_dH_clay = out[13]
c01_dH_diss = out[14]
c01_ph = out[15]
c01_al = out[16]
c01_na = out[17]
c01_si = out[18]
c01_kaolinite0 = out[19]
c01_albite0 = out[20]
c01_mg = out[21]
c01_montna = out[22]
c01_dolomite = out[23]
c01_calcite = out[24]
c01_ca = out[25]



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
t32 = np.loadtxt('r10kt32c02.txt')
t34 = np.loadtxt('r10kt34c02.txt')
t36 = np.loadtxt('r10kt36c02.txt')
t38 = np.loadtxt('r10kt38c02.txt')
t40 = np.loadtxt('r10kt40c02.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

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
t32 = np.loadtxt('r10kt32c03.txt')
t34 = np.loadtxt('r10kt34c03.txt')
t36 = np.loadtxt('r10kt36c03.txt')
t38 = np.loadtxt('r10kt38c03.txt')
t40 = np.loadtxt('r10kt40c03.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

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
t32 = np.loadtxt('r10kt32c04.txt')
t34 = np.loadtxt('r10kt34c04.txt')
t36 = np.loadtxt('r10kt36c04.txt')
t38 = np.loadtxt('r10kt38c04.txt')
t40 = np.loadtxt('r10kt40c04.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

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


#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 5.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c05.txt')
t04 = np.loadtxt('r10kt04c05.txt')
t06 = np.loadtxt('r10kt06c05.txt')
t08 = np.loadtxt('r10kt08c05.txt')
t10 = np.loadtxt('r10kt10c05.txt')
t12 = np.loadtxt('r10kt12c05.txt')
t14 = np.loadtxt('r10kt14c05.txt')
t16 = np.loadtxt('r10kt16c05.txt')
t18 = np.loadtxt('r10kt18c05.txt')
t20 = np.loadtxt('r10kt20c05.txt')
t22 = np.loadtxt('r10kt22c05.txt')
t24 = np.loadtxt('r10kt24c05.txt')
t26 = np.loadtxt('r10kt26c05.txt')
t28 = np.loadtxt('r10kt28c05.txt')
t30 = np.loadtxt('r10kt30c05.txt')
t32 = np.loadtxt('r10kt32c05.txt')
t34 = np.loadtxt('r10kt34c05.txt')
t36 = np.loadtxt('r10kt36c05.txt')
t38 = np.loadtxt('r10kt38c05.txt')
t40 = np.loadtxt('r10kt40c05.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c05_alk = out[0]
c05_alkflux = out[1]
c05_alkflux0 = out[2]
c05_glass = out[3]
c05_water = out[4]
c05_HCO3 = out[5]
c05_CO3 = out[6]
c05_kaolinite = out[7]
c05_stilbite = out[8]
c05_saponite = out[9]
c05_albite = out[10]
c05_celadonite = out[11]
c05_quartz = out[12]
c05_dH_clay = out[13]
c05_dH_diss = out[14]
c05_ph = out[15]
c05_al = out[16]
c05_na = out[17]
c05_si = out[18]
c05_kaolinite0 = out[19]
c05_albite0 = out[20]
c05_mg = out[21]
c05_montna = out[22]
c05_dolomite = out[23]
c05_calcite = out[24]
c05_ca = out[25]



#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 6.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c06.txt')
t04 = np.loadtxt('r10kt04c06.txt')
t06 = np.loadtxt('r10kt06c06.txt')
t08 = np.loadtxt('r10kt08c06.txt')
t10 = np.loadtxt('r10kt10c06.txt')
t12 = np.loadtxt('r10kt12c06.txt')
t14 = np.loadtxt('r10kt14c06.txt')
t16 = np.loadtxt('r10kt16c06.txt')
t18 = np.loadtxt('r10kt18c06.txt')
t20 = np.loadtxt('r10kt20c06.txt')
t22 = np.loadtxt('r10kt22c06.txt')
t24 = np.loadtxt('r10kt24c06.txt')
t26 = np.loadtxt('r10kt26c06.txt')
t28 = np.loadtxt('r10kt28c06.txt')
t30 = np.loadtxt('r10kt30c06.txt')
t32 = np.loadtxt('r10kt32c06.txt')
t34 = np.loadtxt('r10kt34c06.txt')
t36 = np.loadtxt('r10kt36c06.txt')
t38 = np.loadtxt('r10kt38c06.txt')
t40 = np.loadtxt('r10kt40c06.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c06_alk = out[0]
c06_alkflux = out[1]
c06_alkflux0 = out[2]
c06_glass = out[3]
c06_water = out[4]
c06_HCO3 = out[5]
c06_CO3 = out[6]
c06_kaolinite = out[7]
c06_stilbite = out[8]
c06_saponite = out[9]
c06_albite = out[10]
c06_celadonite = out[11]
c06_quartz = out[12]
c06_dH_clay = out[13]
c06_dH_diss = out[14]
c06_ph = out[15]
c06_al = out[16]
c06_na = out[17]
c06_si = out[18]
c06_kaolinite0 = out[19]
c06_albite0 = out[20]
c06_mg = out[21]
c06_montna = out[22]
c06_dolomite = out[23]
c06_calcite = out[24]
c06_ca = out[25]


#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 7.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c07.txt')
t04 = np.loadtxt('r10kt04c07.txt')
t06 = np.loadtxt('r10kt06c07.txt')
t08 = np.loadtxt('r10kt08c07.txt')
t10 = np.loadtxt('r10kt10c07.txt')
t12 = np.loadtxt('r10kt12c07.txt')
t14 = np.loadtxt('r10kt14c07.txt')
t16 = np.loadtxt('r10kt16c07.txt')
t18 = np.loadtxt('r10kt18c07.txt')
t20 = np.loadtxt('r10kt20c07.txt')
t22 = np.loadtxt('r10kt22c07.txt')
t24 = np.loadtxt('r10kt24c07.txt')
t26 = np.loadtxt('r10kt26c07.txt')
t28 = np.loadtxt('r10kt28c07.txt')
t30 = np.loadtxt('r10kt30c07.txt')
t32 = np.loadtxt('r10kt32c07.txt')
t34 = np.loadtxt('r10kt34c07.txt')
t36 = np.loadtxt('r10kt36c07.txt')
t38 = np.loadtxt('r10kt38c07.txt')
t40 = np.loadtxt('r10kt40c07.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c07_alk = out[0]
c07_alkflux = out[1]
c07_alkflux0 = out[2]
c07_glass = out[3]
c07_water = out[4]
c07_HCO3 = out[5]
c07_CO3 = out[6]
c07_kaolinite = out[7]
c07_stilbite = out[8]
c07_saponite = out[9]
c07_albite = out[10]
c07_celadonite = out[11]
c07_quartz = out[12]
c07_dH_clay = out[13]
c07_dH_diss = out[14]
c07_ph = out[15]
c07_al = out[16]
c07_na = out[17]
c07_si = out[18]
c07_kaolinite0 = out[19]
c07_albite0 = out[20]
c07_mg = out[21]
c07_montna = out[22]
c07_dolomite = out[23]
c07_calcite = out[24]
c07_ca = out[25]




#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 8.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c08.txt')
t04 = np.loadtxt('r10kt04c08.txt')
t06 = np.loadtxt('r10kt06c08.txt')
t08 = np.loadtxt('r10kt08c08.txt')
t10 = np.loadtxt('r10kt10c08.txt')
t12 = np.loadtxt('r10kt12c08.txt')
t14 = np.loadtxt('r10kt14c08.txt')
t16 = np.loadtxt('r10kt16c08.txt')
t18 = np.loadtxt('r10kt18c08.txt')
t20 = np.loadtxt('r10kt20c08.txt')
t22 = np.loadtxt('r10kt22c08.txt')
t24 = np.loadtxt('r10kt24c08.txt')
t26 = np.loadtxt('r10kt26c08.txt')
t28 = np.loadtxt('r10kt28c08.txt')
t30 = np.loadtxt('r10kt30c08.txt')
t32 = np.loadtxt('r10kt32c08.txt')
t34 = np.loadtxt('r10kt34c08.txt')
t36 = np.loadtxt('r10kt36c08.txt')
t38 = np.loadtxt('r10kt38c08.txt')
t40 = np.loadtxt('r10kt40c08.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c08_alk = out[0]
c08_alkflux = out[1]
c08_alkflux0 = out[2]
c08_glass = out[3]
c08_water = out[4]
c08_HCO3 = out[5]
c08_CO3 = out[6]
c08_kaolinite = out[7]
c08_stilbite = out[8]
c08_saponite = out[9]
c08_albite = out[10]
c08_celadonite = out[11]
c08_quartz = out[12]
c08_dH_clay = out[13]
c08_dH_diss = out[14]
c08_ph = out[15]
c08_al = out[16]
c08_na = out[17]
c08_si = out[18]
c08_kaolinite0 = out[19]
c08_albite0 = out[20]
c08_mg = out[21]
c08_montna = out[22]
c08_dolomite = out[23]
c08_calcite = out[24]
c08_ca = out[25]




#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 9.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r10kt02c09.txt')
t04 = np.loadtxt('r10kt04c09.txt')
t06 = np.loadtxt('r10kt06c09.txt')
t08 = np.loadtxt('r10kt08c09.txt')
t10 = np.loadtxt('r10kt10c09.txt')
t12 = np.loadtxt('r10kt12c09.txt')
t14 = np.loadtxt('r10kt14c09.txt')
t16 = np.loadtxt('r10kt16c09.txt')
t18 = np.loadtxt('r10kt18c09.txt')
t20 = np.loadtxt('r10kt20c09.txt')
t22 = np.loadtxt('r10kt22c09.txt')
t24 = np.loadtxt('r10kt24c09.txt')
t26 = np.loadtxt('r10kt26c09.txt')
t28 = np.loadtxt('r10kt28c09.txt')
t30 = np.loadtxt('r10kt30c09.txt')
t32 = np.loadtxt('r10kt32c09.txt')
t34 = np.loadtxt('r10kt34c09.txt')
t36 = np.loadtxt('r10kt36c09.txt')
t38 = np.loadtxt('r10kt38c09.txt')
t40 = np.loadtxt('r10kt40c09.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c09_alk = out[0]
c09_alkflux = out[1]
c09_alkflux0 = out[2]
c09_glass = out[3]
c09_water = out[4]
c09_HCO3 = out[5]
c09_CO3 = out[6]
c09_kaolinite = out[7]
c09_stilbite = out[8]
c09_saponite = out[9]
c09_albite = out[10]
c09_celadonite = out[11]
c09_quartz = out[12]
c09_dH_clay = out[13]
c09_dH_diss = out[14]
c09_ph = out[15]
c09_al = out[16]
c09_na = out[17]
c09_si = out[18]
c09_kaolinite0 = out[19]
c09_albite0 = out[20]
c09_mg = out[21]
c09_montna = out[22]
c09_dolomite = out[23]
c09_calcite = out[24]
c09_ca = out[25]




#############################
# t_RES = 3.14e10 (10kyr)   #
# DIC = 10.0 mmol/kgw       #
#############################

# load temps
t02 = np.loadtxt('r10kt02c10.txt')
t04 = np.loadtxt('r10kt04c10.txt')
t06 = np.loadtxt('r10kt06c10.txt')
t08 = np.loadtxt('r10kt08c10.txt')
t10 = np.loadtxt('r10kt10c10.txt')
t12 = np.loadtxt('r10kt12c10.txt')
t14 = np.loadtxt('r10kt14c10.txt')
t16 = np.loadtxt('r10kt16c10.txt')
t18 = np.loadtxt('r10kt18c10.txt')
t20 = np.loadtxt('r10kt20c10.txt')
t22 = np.loadtxt('r10kt22c10.txt')
t24 = np.loadtxt('r10kt24c10.txt')
t26 = np.loadtxt('r10kt26c10.txt')
t28 = np.loadtxt('r10kt28c10.txt')
t30 = np.loadtxt('r10kt30c10.txt')
t32 = np.loadtxt('r10kt32c10.txt')
t34 = np.loadtxt('r10kt34c10.txt')
t36 = np.loadtxt('r10kt36c10.txt')
t38 = np.loadtxt('r10kt38c10.txt')
t40 = np.loadtxt('r10kt40c10.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c10_alk = out[0]
c10_alkflux = out[1]
c10_alkflux0 = out[2]
c10_glass = out[3]
c10_water = out[4]
c10_HCO3 = out[5]
c10_CO3 = out[6]
c10_kaolinite = out[7]
c10_stilbite = out[8]
c10_saponite = out[9]
c10_albite = out[10]
c10_celadonite = out[11]
c10_quartz = out[12]
c10_dH_clay = out[13]
c10_dH_diss = out[14]
c10_ph = out[15]
c10_al = out[16]
c10_na = out[17]
c10_si = out[18]
c10_kaolinite0 = out[19]
c10_albite0 = out[20]
c10_mg = out[21]
c10_montna = out[22]
c10_dolomite = out[23]
c10_calcite = out[24]
c10_ca = out[25]







##############################
# PLOT PCOLOR 2D PARAM SPACE #
##############################

fig=plt.figure()

ax = plt.subplot(1,1,1)

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

dic = np.array([.001, .002, .003, .004, .005, .006, .007, .008, .009, .010])

grid = np.array([c01_calcite, c02_calcite, c03_calcite, c04_calcite,
                 c05_calcite, c06_calcite, c07_calcite, c08_calcite,
                 c09_calcite, c10_calcite])

p = plt.contourf(dic, temps, np.transpose(grid),
               cmap=cm.Spectral_r, edgecolors='#444444', linewidth=2)


plt.xlabel('SEAWATER DIC [mol kgw$^{-1}$]',fontsize=10)
plt.ylabel('T [$^{\circ}$C]',fontsize=10)

#plt.xticks(dic+.0005,dic)
#plt.yticks(temps[::-1]+1.0,temps[::-1])

plt.xticks(dic,dic)
plt.yticks(temps[::-1],temps[::-1])

#plt.xlim([.002,.005])
plt.ylim([40.0,2.0])

cbar = plt.colorbar(orientation='horizontal')
cbar.ax.set_xlabel('CALCITE PRECIPITATION RATE [mol yr$^{-1}$]')

plt.savefig('pcolor0.png')









############################
# PLOT CALCITE GROWTH RATE #
############################


fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)


ax = plt.subplot(2,2,1)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,c01_calcite*100.0,linewidth=2,label='DIC = 1.0 mmol/kgw')
p = plt.plot(temps,c02_calcite*100.0,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03_calcite*100.0,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04_calcite*100.0,linewidth=2,label='DIC = 4.0 mmol/kgw')
p = plt.plot(temps,c05_calcite*100.0,linewidth=2,label='DIC = 5.0 mmol/kgw')
p = plt.plot(temps,c06_calcite*100.0,linewidth=2,label='DIC = 6.0 mmol/kgw')


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

plt.title('GLASS',fontsize=10)
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
