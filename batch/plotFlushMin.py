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
        t_na[i] = bit[0,m,6]
        t_si[i] = bit[0,m,10]

        t_kaolinite0[i] = bit[0,n,19]
        t_albite0[i] = bit[0,n,21]

        t_mg[i] = bit[0,m,5]
        t_montna[i] = bit[0,m,31]
        
        t_dolomite[i] = bit[0,m,35]
        t_calcite[i] = bit[0,m,45]

        t_ca[i] = bit[0,m,4]

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
# t_RES = 3.14e9 (109yr)    #
# DIC = 1.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c01mm.txt')
t04 = np.loadtxt('r100t04c01mm.txt')
t06 = np.loadtxt('r100t06c01mm.txt')
t08 = np.loadtxt('r100t08c01mm.txt')
t10 = np.loadtxt('r100t10c01mm.txt')
t12 = np.loadtxt('r100t12c01mm.txt')
t14 = np.loadtxt('r100t14c01mm.txt')
t16 = np.loadtxt('r100t16c01mm.txt')
t18 = np.loadtxt('r100t18c01mm.txt')
t20 = np.loadtxt('r100t20c01mm.txt')
t22 = np.loadtxt('r100t22c01mm.txt')
t24 = np.loadtxt('r100t24c01mm.txt')
t26 = np.loadtxt('r100t26c01mm.txt')
t28 = np.loadtxt('r100t28c01mm.txt')
t30 = np.loadtxt('r100t30c01mm.txt')
t32 = np.loadtxt('r100t32c01mm.txt')
t34 = np.loadtxt('r100t34c01mm.txt')
t36 = np.loadtxt('r100t36c01mm.txt')
t38 = np.loadtxt('r100t38c01mm.txt')
t40 = np.loadtxt('r100t40c01mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c01mm_alk = out[0]
c01mm_alkflux = out[1]
c01mm_alkflux0 = out[2]
c01mm_glass = out[3]
c01mm_water = out[4]
c01mm_HCO3 = out[5]
c01mm_CO3 = out[6]
c01mm_kaolinite = out[7]
c01mm_stilbite = out[8]
c01mm_saponite = out[9]
c01mm_albite = out[10]
c01mm_celadonite = out[11]
c01mm_quartz = out[12]
c01mm_dH_clay = out[13]
c01mm_dH_diss = out[14]
c01mm_ph = out[15]
c01mm_al = out[16]
c01mm_na = out[17]
c01mm_si = out[18]
c01mm_kaolinite0 = out[19]
c01mm_albite0 = out[20]
c01mm_mg = out[21]
c01mm_montna = out[22]
c01mm_dolomite = out[23]
c01mm_calcite = out[24]
c01mm_ca = out[25]




#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 2.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c02mm.txt')
t04 = np.loadtxt('r100t04c02mm.txt')
t06 = np.loadtxt('r100t06c02mm.txt')
t08 = np.loadtxt('r100t08c02mm.txt')
t10 = np.loadtxt('r100t10c02mm.txt')
t12 = np.loadtxt('r100t12c02mm.txt')
t14 = np.loadtxt('r100t14c02mm.txt')
t16 = np.loadtxt('r100t16c02mm.txt')
t18 = np.loadtxt('r100t18c02mm.txt')
t20 = np.loadtxt('r100t20c02mm.txt')
t22 = np.loadtxt('r100t22c02mm.txt')
t24 = np.loadtxt('r100t24c02mm.txt')
t26 = np.loadtxt('r100t26c02mm.txt')
t28 = np.loadtxt('r100t28c02mm.txt')
t30 = np.loadtxt('r100t30c02mm.txt')
t32 = np.loadtxt('r100t32c02mm.txt')
t34 = np.loadtxt('r100t34c02mm.txt')
t36 = np.loadtxt('r100t36c02mm.txt')
t38 = np.loadtxt('r100t38c02mm.txt')
t40 = np.loadtxt('r100t40c02mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c02mm_alk = out[0]
c02mm_alkflux = out[1]
c02mm_alkflux0 = out[2]
c02mm_glass = out[3]
c02mm_water = out[4]
c02mm_HCO3 = out[5]
c02mm_CO3 = out[6]
c02mm_kaolinite = out[7]
c02mm_stilbite = out[8]
c02mm_saponite = out[9]
c02mm_albite = out[10]
c02mm_celadonite = out[11]
c02mm_quartz = out[12]
c02mm_dH_clay = out[13]
c02mm_dH_diss = out[14]
c02mm_ph = out[15]
c02mm_al = out[16]
c02mm_na = out[17]
c02mm_si = out[18]
c02mm_kaolinite0 = out[19]
c02mm_albite0 = out[20]
c02mm_mg = out[21]
c02mm_montna = out[22]
c02mm_dolomite = out[23]
c02mm_calcite = out[24]
c02mm_ca = out[25]



#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 3.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c03mm.txt')
t04 = np.loadtxt('r100t04c03mm.txt')
t06 = np.loadtxt('r100t06c03mm.txt')
t08 = np.loadtxt('r100t08c03mm.txt')
t10 = np.loadtxt('r100t10c03mm.txt')
t12 = np.loadtxt('r100t12c03mm.txt')
t14 = np.loadtxt('r100t14c03mm.txt')
t16 = np.loadtxt('r100t16c03mm.txt')
t18 = np.loadtxt('r100t18c03mm.txt')
t20 = np.loadtxt('r100t20c03mm.txt')
t22 = np.loadtxt('r100t22c03mm.txt')
t24 = np.loadtxt('r100t24c03mm.txt')
t26 = np.loadtxt('r100t26c03mm.txt')
t28 = np.loadtxt('r100t28c03mm.txt')
t30 = np.loadtxt('r100t30c03mm.txt')
t32 = np.loadtxt('r100t32c03mm.txt')
t34 = np.loadtxt('r100t32c03mm.txt') ####### 32 repeat
t36 = np.loadtxt('r100t36c03mm.txt')
t38 = np.loadtxt('r100t38c03mm.txt')
t40 = np.loadtxt('r100t40c03mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c03mm_alk = out[0]
c03mm_alkflux = out[1]
c03mm_alkflux0 = out[2]
c03mm_glass = out[3]
c03mm_water = out[4]
c03mm_HCO3 = out[5]
c03mm_CO3 = out[6]
c03mm_kaolinite = out[7]
c03mm_stilbite = out[8]
c03mm_saponite = out[9]
c03mm_albite = out[10]
c03mm_celadonite = out[11]
c03mm_quartz = out[12]
c03mm_dH_clay = out[13]
c03mm_dH_diss = out[14]
c03mm_ph = out[15]
c03mm_al = out[16]
c03mm_na = out[17]
c03mm_si = out[18]
c03mm_kaolinite0 = out[19]
c03mm_albite0 = out[20]
c03mm_mg = out[21]
c03mm_montna = out[22]
c03mm_dolomite = out[23]
c03mm_calcite = out[24]
c03mm_ca = out[25]



#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 4.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c04mm.txt')
t04 = np.loadtxt('r100t04c04mm.txt')
t06 = np.loadtxt('r100t06c04mm.txt')
t08 = np.loadtxt('r100t08c04mm.txt')
t10 = np.loadtxt('r100t10c04mm.txt')
t12 = np.loadtxt('r100t12c04mm.txt')
t14 = np.loadtxt('r100t14c04mm.txt')
t16 = np.loadtxt('r100t16c04mm.txt')
t18 = np.loadtxt('r100t18c04mm.txt')
t20 = np.loadtxt('r100t20c04mm.txt')
t22 = np.loadtxt('r100t22c04mm.txt')
t24 = np.loadtxt('r100t24c04mm.txt')
t26 = np.loadtxt('r100t26c04mm.txt')
t28 = np.loadtxt('r100t28c04mm.txt')
t30 = np.loadtxt('r100t30c04mm.txt')
t32 = np.loadtxt('r100t32c04mm.txt')
t34 = np.loadtxt('r100t34c04mm.txt')
t36 = np.loadtxt('r100t36c04mm.txt')
t38 = np.loadtxt('r100t38c04mm.txt')
t40 = np.loadtxt('r100t40c04mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c04mm_alk = out[0]
c04mm_alkflux = out[1]
c04mm_alkflux0 = out[2]
c04mm_glass = out[3]
c04mm_water = out[4]
c04mm_HCO3 = out[5]
c04mm_CO3 = out[6]
c04mm_kaolinite = out[7]
c04mm_stilbite = out[8]
c04mm_saponite = out[9]
c04mm_albite = out[10]
c04mm_celadonite = out[11]
c04mm_quartz = out[12]
c04mm_dH_clay = out[13]
c04mm_dH_diss = out[14]
c04mm_ph = out[15]
c04mm_al = out[16]
c04mm_na = out[17]
c04mm_si = out[18]
c04mm_kaolinite0 = out[19]
c04mm_albite0 = out[20]
c04mm_mg = out[21]
c04mm_montna = out[22]
c04mm_dolomite = out[23]
c04mm_calcite = out[24]
c04mm_ca = out[25]


#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 5.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c05mm.txt')
t04 = np.loadtxt('r100t04c05mm.txt')
t06 = np.loadtxt('r100t06c05mm.txt')
t08 = np.loadtxt('r100t08c05mm.txt')
t10 = np.loadtxt('r100t10c05mm.txt')
t12 = np.loadtxt('r100t12c05mm.txt')
t14 = np.loadtxt('r100t14c05mm.txt')
t16 = np.loadtxt('r100t16c05mm.txt')
t18 = np.loadtxt('r100t18c05mm.txt')
t20 = np.loadtxt('r100t20c05mm.txt')
t22 = np.loadtxt('r100t22c05mm.txt')
t24 = np.loadtxt('r100t24c05mm.txt')
t26 = np.loadtxt('r100t26c05mm.txt')
t28 = np.loadtxt('r100t28c05mm.txt')
t30 = np.loadtxt('r100t30c05mm.txt')
t32 = np.loadtxt('r100t32c05mm.txt')
t34 = np.loadtxt('r100t34c05mm.txt')
t36 = np.loadtxt('r100t36c05mm.txt')
t38 = np.loadtxt('r100t36c05mm.txt') ### REPEAT?
t40 = np.loadtxt('r100t40c05mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c05mm_alk = out[0]
c05mm_alkflux = out[1]
c05mm_alkflux0 = out[2]
c05mm_glass = out[3]
c05mm_water = out[4]
c05mm_HCO3 = out[5]
c05mm_CO3 = out[6]
c05mm_kaolinite = out[7]
c05mm_stilbite = out[8]
c05mm_saponite = out[9]
c05mm_albite = out[10]
c05mm_celadonite = out[11]
c05mm_quartz = out[12]
c05mm_dH_clay = out[13]
c05mm_dH_diss = out[14]
c05mm_ph = out[15]
c05mm_al = out[16]
c05mm_na = out[17]
c05mm_si = out[18]
c05mm_kaolinite0 = out[19]
c05mm_albite0 = out[20]
c05mm_mg = out[21]
c05mm_montna = out[22]
c05mm_dolomite = out[23]
c05mm_calcite = out[24]
c05mm_ca = out[25]




#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 6.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c06mm.txt')
t04 = np.loadtxt('r100t04c06mm.txt')
t06 = np.loadtxt('r100t06c06mm.txt')
t08 = np.loadtxt('r100t08c06mm.txt')
t10 = np.loadtxt('r100t10c06mm.txt')
t12 = np.loadtxt('r100t12c06mm.txt')
t14 = np.loadtxt('r100t14c06mm.txt')
t16 = np.loadtxt('r100t16c06mm.txt')
t18 = np.loadtxt('r100t18c06mm.txt')
t20 = np.loadtxt('r100t20c06mm.txt')
t22 = np.loadtxt('r100t22c06mm.txt')
t24 = np.loadtxt('r100t24c06mm.txt')
t26 = np.loadtxt('r100t26c06mm.txt')
t28 = np.loadtxt('r100t28c06mm.txt')
t30 = np.loadtxt('r100t30c06mm.txt')
t32 = np.loadtxt('r100t32c06mm.txt')
t34 = np.loadtxt('r100t34c06mm.txt')
t36 = np.loadtxt('r100t36c06mm.txt')
t38 = np.loadtxt('r100t38c06mm.txt') 
t40 = np.loadtxt('r100t38c06mm.txt') ## REPEAT ??


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c06mm_alk = out[0]
c06mm_alkflux = out[1]
c06mm_alkflux0 = out[2]
c06mm_glass = out[3]
c06mm_water = out[4]
c06mm_HCO3 = out[5]
c06mm_CO3 = out[6]
c06mm_kaolinite = out[7]
c06mm_stilbite = out[8]
c06mm_saponite = out[9]
c06mm_albite = out[10]
c06mm_celadonite = out[11]
c06mm_quartz = out[12]
c06mm_dH_clay = out[13]
c06mm_dH_diss = out[14]
c06mm_ph = out[15]
c06mm_al = out[16]
c06mm_na = out[17]
c06mm_si = out[18]
c06mm_kaolinite0 = out[19]
c06mm_albite0 = out[20]
c06mm_mg = out[21]
c06mm_montna = out[22]
c06mm_dolomite = out[23]
c06mm_calcite = out[24]
c06mm_ca = out[25]




#############################
# t_RES = 3.14e9 (100yr)    #
# DIC = 7.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c07mm.txt')
t04 = np.loadtxt('r100t04c07mm.txt')
t06 = np.loadtxt('r100t06c07mm.txt')
t08 = np.loadtxt('r100t08c07mm.txt')
t10 = np.loadtxt('r100t10c07mm.txt')
t12 = np.loadtxt('r100t12c07mm.txt')
t14 = np.loadtxt('r100t14c07mm.txt')
t16 = np.loadtxt('r100t16c07mm.txt')
t18 = np.loadtxt('r100t18c07mm.txt')
t20 = np.loadtxt('r100t20c07mm.txt')
t22 = np.loadtxt('r100t22c07mm.txt')
t24 = np.loadtxt('r100t24c07mm.txt')
t26 = np.loadtxt('r100t26c07mm.txt')
t28 = np.loadtxt('r100t28c07mm.txt')
t30 = np.loadtxt('r100t30c07mm.txt')
t32 = np.loadtxt('r100t32c07mm.txt')
t34 = np.loadtxt('r100t34c07mm.txt')
t36 = np.loadtxt('r100t36c07mm.txt')
t38 = np.loadtxt('r100t36c07mm.txt') ### REPEAT?
t40 = np.loadtxt('r100t40c07mm.txt')


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c07mm_alk = out[0]
c07mm_alkflux = out[1]
c07mm_alkflux0 = out[2]
c07mm_glass = out[3]
c07mm_water = out[4]
c07mm_HCO3 = out[5]
c07mm_CO3 = out[6]
c07mm_kaolinite = out[7]
c07mm_stilbite = out[8]
c07mm_saponite = out[9]
c07mm_albite = out[10]
c07mm_celadonite = out[11]
c07mm_quartz = out[12]
c07mm_dH_clay = out[13]
c07mm_dH_diss = out[14]
c07mm_ph = out[15]
c07mm_al = out[16]
c07mm_na = out[17]
c07mm_si = out[18]
c07mm_kaolinite0 = out[19]
c07mm_albite0 = out[20]
c07mm_mg = out[21]
c07mm_montna = out[22]
c07mm_dolomite = out[23]
c07mm_calcite = out[24]
c07mm_ca = out[25]






#############################
# t_RES = 3.14e9 (109yr)    #
# DIC = 8.0 mmol/kgw        #
#############################

# load temps
t02 = np.loadtxt('r100t02c08mm.txt')
t04 = np.loadtxt('r100t04c08mm.txt')
t06 = np.loadtxt('r100t06c08mm.txt')
t08 = np.loadtxt('r100t08c08mm.txt')
t10 = np.loadtxt('r100t10c08mm.txt')
t12 = np.loadtxt('r100t12c08mm.txt')
t14 = np.loadtxt('r100t14c08mm.txt')
t16 = np.loadtxt('r100t16c08mm.txt')
t18 = np.loadtxt('r100t18c08mm.txt')
t20 = np.loadtxt('r100t20c08mm.txt')
t22 = np.loadtxt('r100t22c08mm.txt')
t24 = np.loadtxt('r100t24c08mm.txt')
t26 = np.loadtxt('r100t26c08mm.txt')
t28 = np.loadtxt('r100t28c08mm.txt')
t30 = np.loadtxt('r100t30c08mm.txt')
t32 = np.loadtxt('r100t32c08mm.txt')
t34 = np.loadtxt('r100t34c08mm.txt')
t36 = np.loadtxt('r100t36c08mm.txt')
t38 = np.loadtxt('r100t38c08mm.txt') 
t40 = np.loadtxt('r100t38c08mm.txt') ### REPEAT?


t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]],
     [t22[:,:]], [t24[:,:]], [t26[:,:]], [t28[:,:]], [t30[:,:]],
     [t32[:,:]], [t34[:,:]], [t36[:,:]], [t38[:,:]], [t40[:,:]]]

out = grab(t)

c08mm_alk = out[0]
c08mm_alkflux = out[1]
c08mm_alkflux0 = out[2]
c08mm_glass = out[3]
c08mm_water = out[4]
c08mm_HCO3 = out[5]
c08mm_CO3 = out[6]
c08mm_kaolinite = out[7]
c08mm_stilbite = out[8]
c08mm_saponite = out[9]
c08mm_albite = out[10]
c08mm_celadonite = out[11]
c08mm_quartz = out[12]
c08mm_dH_clay = out[13]
c08mm_dH_diss = out[14]
c08mm_ph = out[15]
c08mm_al = out[16]
c08mm_na = out[17]
c08mm_si = out[18]
c08mm_kaolinite0 = out[19]
c08mm_albite0 = out[20]
c08mm_mg = out[21]
c08mm_montna = out[22]
c08mm_dolomite = out[23]
c08mm_calcite = out[24]
c08mm_ca = out[25]





##############################
# PLOT PCOLOR 2D PARAM SPACE #
##############################

fig=plt.figure()

ax = plt.subplot(1,1,1)

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

dic = np.array([.001])

grid = np.array([c01mm_alkflux])

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









########################
# PLOT ALK FLUX (DIFF) #
########################


fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)


ax = plt.subplot(2,2,1)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,c01mm_dolomite,linewidth=2,label='DIC = 1.0 mmol/kgw')
p = plt.plot(temps,c02mm_dolomite,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03mm_dolomite,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04mm_dolomite,linewidth=2,label='DIC = 4.0 mmol/kgw')
p = plt.plot(temps,c05mm_dolomite,linewidth=2,label='DIC = 5.0 mmol/kgw')
p = plt.plot(temps,c06mm_dolomite,linewidth=2,label='DIC = 6.0 mmol/kgw')
p = plt.plot(temps,c07mm_dolomite,linewidth=2,label='DIC = 7.0 mmol/kgw')


plt.title('ALKALINITY FLUX',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



ax = plt.subplot(2,2,2)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,c01mm_alkflux,linewidth=2,label='DIC = 1.0 mmol/kgw')
p = plt.plot(temps,c02mm_alkflux,linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,c03mm_alkflux,linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,c04mm_alkflux,linewidth=2,label='DIC = 4.0 mmol/kgw')
p = plt.plot(temps,c05mm_alkflux,linewidth=2,label='DIC = 5.0 mmol/kgw')
p = plt.plot(temps,c06mm_alkflux,linewidth=2,label='DIC = 6.0 mmol/kgw')
p = plt.plot(temps,c07mm_alkflux,linewidth=2,label='DIC = 7.0 mmol/kgw')

plt.title('ALKALINITY FLUX',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
#plt.legend(handles, labels,loc=3,prop={'size':6}, ncol=1)




ax = plt.subplot(2,2,3)

p = plt.plot([0.0,40.0], [0.0,0.0], 'k:')
p = plt.plot(temps,.0077*c01mm_glass,'r',linewidth=2,label='DIC = 1.0 mmol/kgw')
p = plt.plot(temps,.0077*c02mm_glass,'g',linewidth=2,label='DIC = 2.0 mmol/kgw')
p = plt.plot(temps,.0077*c03mm_glass,'b',linewidth=2,label='DIC = 3.0 mmol/kgw')
p = plt.plot(temps,.0077*c04mm_glass,'c',linewidth=2,label='DIC = 4.0 mmol/kgw')
p = plt.plot(temps,.0077*c05mm_glass,'gold',linewidth=2,label='DIC = 5.0 mmol/kgw')

plt.title('Ca2+',fontsize=10)
plt.ylabel('[eq kgw$^{-1}$ yr$^{-1}$]',fontsize=10)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':10}, ncol=1)



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
