# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:42:51 2016

@author: sig

https://smartech.gatech.edu/bitstream/handle/1853/6820/dukeman_greg_a_200505_phd.pdf

Calculate rocket launch trajectory example"""

import numpy as np
import matplotlib.pyplot as plt
import copy 
from matplotlib import _cntr as cntr

execfile('lib_physicalconstants.py')
execfile('class_vessel.py')
execfile('poweredExplicitGuidance.py')

#function definitions
def try_AscentSimulation(vessel, velocities, pitches):
  #setup meshgrid and number of simulations
  [v_, p_] = np.meshgrid(velocities, pitches)  
  a_ = np.empty_like(v_)    #apoapsis [m]
  q_ = np.empty_like(v_)    #maxq [kg/ms2]
  l_ = np.empty_like(v_)    #delta-v losses combined drag and gravity [m/s]
  m_ = np.empty_like(v_)    #vessel endmass try to maximize
  
  #iterate through all velocities and pitches and store apoapsis and maxq
  it = np.nditer(v_, flags=['multi_index'])
  while not it.finished:
    ves = copy.deepcopy(vessel)             #initiate/create new vessel NOTE this creates a copy of vessel and leaves the original vessel untouched   
    v = it[0]
    p = p_[it.multi_index]
    
    #do the simulation with initial v and p
    #to a target altitude of 200 km apoapis and circular orbit
    [trajectory, PEGparameters] = ves.AscentSimulation(v, p,  R+ 200000,0)

    #store maxq value and apoapsis from getOrbital() function    
    if ves.state == 's2':    
      q_[it.multi_index] = ves.maxq
      orbit = ves._getOrbital()
      a_[it.multi_index] = orbit[0]
      l_[it.multi_index] = np.abs(ves.vgloss) + np.abs(ves.vdloss)
      m_[it.multi_index] = ves.m
    else:
      m_[it.multi_index] = ves.s2_m0 - ves.s2_mp
    
    print it.multi_index   
    it.iternext()
  
  return v_, p_, a_, q_, l_, m_


#define ships drag curfe
Cd = np.array([[ 0, 0.122], [256.0, 0.122], [343.2, 0.583], [643.5, 0.558], [909.5, 0.554], [1673,  0.676],[9999,  0.776]])

ship = VESSEL(
   s1_m0 = 97198,             #stage 1: lauchmass [kg]
   s1_mp = 75744 ,            #stage 1: propellant mass [kg]                
   s1_thrust_asl = 1217150,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl = 230,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 250,          #stage 1: specific impulse isp in vaccum [s]
   s1_A = 7.06,               #stage 1: reference cross sectional area [m2]
   s2_m0 = 7442,              #stage 2: mass at separation [kg]
   s2_mp = 6284,              #stage 2: propellant mass [kg]
   s2_thrust = 55400,         #stage 2: thrust in vaccum [N]
   s2_isp = 340,              #stage 2: specific impulse in vaccum [s]
   s2_A = 2.0,                #stage 2: reference cross sectional area [m2]
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

ves2s = VESSEL(
   s1_m0 = 138855,            #stage 1: lauchmass [kg]
   s1_mp = 90603*0.96 ,        #stage 1: propellant mass [kg]                
   s1_thrust_asl = 1777537,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl = 252,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 290,          #stage 1: specific impulse isp in vaccum [s]
   s1_A = (3.05/2)*np.pi,      #stage 1: reference cross sectional area [m2]
   s2_m0 = 7442,              #stage 2: mass at separation [kg]
   s2_mp = 6284,              #stage 2: propellant mass [kg]
   s2_thrust = 55400,         #stage 2: thrust in vaccum [N]
   s2_isp = 340,              #stage 2: specific impulse in vaccum [s]
   s2_A = 2.0,                #stage 2: reference cross sectional area [m2]
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

ship.add_launchsite(441,45.9)    #add launchsite at altitude of 441 m and latitude of 45.9 deg

[trajectory, PEGparameters] = ship.AscentSimulation(50, 87, R+ 200000, 0,2)

#find trajectories with 150 km target apoapsis
#try_velocities = np.linspace(40,80,10)
#try_pitches = np.linspace(83,89,14)

#[V_,P_,A_,Q_, L_, M_] = try_AscentSimulation(ship, try_velocities, try_pitches)

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,1])
plt.xlabel('Time [s]')
plt.ylabel('Tangential velocity (orbit) [m/s]')

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,2])
plt.xlabel('Time [s]')
plt.ylabel('Radial (vertical) velocity [m/s]')

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,3])
plt.xlabel('Time [s]')
plt.ylabel('Altitude [m]')

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,4])
plt.xlabel('Time [s]')
plt.ylabel('Pitch [deg]')

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,6])
plt.ylabel('Dynamic pressure q [kg/ms2]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(PEGparameters[:,0],PEGparameters[:,1])
plt.ylabel('PEG parameter A [#]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(PEGparameters[:,0],PEGparameters[:,2])
plt.ylabel('PEG parameter B [#]')
plt.xlabel('Time [s]')



