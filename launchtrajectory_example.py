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
import scipy.signal as spsig

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

def writeKSbootscript(vessel, name):
  filename = 'boot_PEG_' + name + '.ks'
  f = open(filename,'w')

  #first write name in variable
  f.write('GLOBAL P_vName IS "' + name + '".\n\n') # python will convert \n to os.linesep
  
  #next write sequence
  f.write('//ignition (delay before release), booster jettison, fairings jettison, cutoff, separation, ullage, ignition, PEG activation, stage 2 maxT\n')
  f.write('GLOBAL P_seq IS LIST(' + ' -%.2f' % (vessel.s1_engT) + ', , , ' + '%.2f' % (vessel.s1_tb - 1.4) + ', , , , , ' + '%.2f' % (vessel.s2_tb) + ').//remember to fill this!\n\n')

  #next write pitchover maneuver list for stage 1
  #create a list of time vs. pitch (linear interpolation) 
  t_ = vessel.trajectory[:,0]
  t_ = np.linspace(t_[0], t_[-1], 50)
  p_ = np.copy(t_)
  
  for p in np.nditer(p_, op_flags=['readwrite']):
    p[...] = np.interp(p, vessel.trajectory[1:,0], vessel.trajectory[1:,4], left = 90.0) 

  f.write('GLOBAL P_pt IS LIST( ')
  for t in t_[:-1]:  
    f.write( ' %.1f, ' % (t) )     
  f.write(' %.1f ).\n' % t_[-1])

  f.write('GLOBAL P_pp IS LIST( ')
  for p in p_[:-1]:  
    f.write( ' %.1f, ' % (p) )     
  f.write(' %.1f ).\n\n' % p_[-1])
  
  #set gravitational constant of body earth
  f.write('GLOBAL P_umode IS .//remember to fill this!\n\n')

  f.write('TOGGLE AG10.\n\n')
  f.write('RUN pegas_loader.ks.\n')

  return 
  
#define ships drag curfe
Cd = np.array([[ 0, 0.122], [256.0, 0.122], [343.2, 0.583], [643.5, 0.558], [909.5, 0.554], [1673,  0.676],[9999,  0.776]])

ship = VESSEL(
   s1_m0 = 97198,             #stage 1: lauchmass [kg]
   s1_mp = 75744 ,            #stage 1: propellant mass [kg]                
   s1_thrust_asl = 1217150,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl = 230,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 250,          #stage 1: specific impulse isp in vaccum [s]
   s1_A = 7.06,               #stage 1: reference cross sectional area [m2]
   s1_engT = 1.5,             #stage 1: waittime until launchtower release engine preingition time [s]
   s2_m0 = 7442,              #stage 2: mass at separation [kg]
   s2_mp = 6284,              #stage 2: propellant mass [kg]
   s2_thrust = 55400,         #stage 2: thrust in vaccum [N]
   s2_isp = 340,              #stage 2: specific impulse in vaccum [s]
   s2_A = 2.0,                #stage 2: reference cross sectional area [m2]
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

ship.add_launchsite(441,45.9)    #add launchsite at altitude of 441 m and latitude of 45.9 deg

ship.AscentSimulation(50, 87, R+ 200000, 0,2)

writeKSbootscript(ship, 'test')


#find trajectories with 150 km target apoapsis
#try_velocities = np.linspace(40,80,10)
#try_pitches = np.linspace(83,89,14)

#[V_,P_,A_,Q_, L_, M_] = try_AscentSimulation(ship, try_velocities, try_pitches)

plt.figure()
plt.plot(ship.trajectory[:,0],ship.trajectory[:,1])
plt.xlabel('Time [s]')
plt.ylabel('Tangential velocity (orbit) [m/s]')

plt.figure()
plt.plot(ship.trajectory[:,0],ship.trajectory[:,2])
plt.xlabel('Time [s]')
plt.ylabel('Radial (vertical) velocity [m/s]')

plt.figure()
plt.plot(ship.trajectory[:,0],ship.trajectory[:,3])
plt.xlabel('Time [s]')
plt.ylabel('Altitude [m]')

plt.figure()
plt.plot(ship.trajectory[:,0],ship.trajectory[:,4])
plt.xlabel('Time [s]')
plt.ylabel('Pitch [deg]')

plt.figure()
plt.plot(ship.trajectory[:,0],ship.trajectory[:,6])
plt.ylabel('Dynamic pressure q [kg/ms2]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(ship.PEGparameters[:,0],ship.PEGparameters[:,1])
plt.ylabel('PEG parameter A [#]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(ship.PEGparameters[:,0],ship.PEGparameters[:,2])
plt.ylabel('PEG parameter B [#]')
plt.xlabel('Time [s]')



