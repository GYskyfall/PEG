# -*- coding: utf-8 -*-
"""
This libary contains all physics parameters and intputs for atmosphere and earth

Created on Wed Feb 17 10:55:19 2016

@author: sig
"""

import numpy as np
import matplotlib.pyplot as plt

#Earth's constants
global g0; g0 = 9.80665           #std gravity asl [m/s2]
global mu; mu = 398600441800000;  #std gravity param for Earth [m3/s2]
global R; R = 6371000;            #earth radius [m]

#rocket launch trajectory calculation based on http://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html

#function to spit out air density at a certain altitude h [m] of earth
def rho( h ):  #air density [kg/m3]
    #troposphere (low atmosphere)
    if h < 11000:
        T = 15.04 - .00649 * h                    #temperature [°C]
        p = 101.29 * np.power((T + 273.1)/288.08,5.256)   #pressure [kilo-Pascals]
    #lower stratosphere
    elif ( 11000 <= h and h < 25000):
        T = -56.46
        p = 22.65 * np.exp(1.73 - .000157 * h)
    #upper atmosphere
    elif (h > 25000):
        T = -131.21 + .00299 * h 
        p = 2.488 * np.power((T + 273.1)/ 216.6,-11.388) 

     
    rho = p / (.2869 * (T + 273.1)) #air density from the equation of state [kg/m3]
   
    return rho

rho = np.vectorize(rho)

#function to spit out air density at a certain altitude h [m] of earth
def atmp( h ):  #air density [kg/m3]
    #troposphere (low atmosphere)
    if h < 11000:
        T = 15.04 - .00649 * h                    #temperature [°C]
        p = 101.29 * np.power((T + 273.1)/288.08,5.256)   #pressure [kilo-Pascals]
    #lower stratosphere
    elif ( 11000 <= h and h < 25000):
        T = -56.46
        p = 22.65 * np.exp(1.73 - .000157 * h)
    #upper atmosphere
    elif (h > 25000):
        T = -131.21 + .00299 * h 
        p = 2.488 * np.power((T + 273.1)/ 216.6,-11.388) 
       
    return p

atmp = np.vectorize(atmp)

#function to spit out air temperature in [C] at a certain altitude h [m] of earth
def atmtemp( h ):  #air temperature [C]
    #troposphere (low atmosphere)
    if h < 11000:
        T = 15.04 - .00649 * h                    #temperature [°C]
    #lower stratosphere
    elif ( 11000 <= h and h < 25000):
        T = -56.46
    #upper atmosphere
    elif (h > 25000):
        T = -131.21 + .00299 * h 
     
    atmtemp = T
   
    return atmtemp

atmtemp = np.vectorize(atmtemp)


#alt = np.linspace(0, 130000,100)
#fig = plt.figure()
#plt.plot(alt,rho(alt)*80)
#plt.plot(alt,atmp(alt))