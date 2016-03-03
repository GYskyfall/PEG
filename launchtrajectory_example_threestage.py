# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:42:51 2016

@author: sig

Recreation of matlab code from https://github.com/Noiredd/PEGAS
https://smartech.gatech.edu/bitstream/handle/1853/6820/dukeman_greg_a_200505_phd.pdf
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660006073.pdf
http://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance

Calculate rocket launch trajectory example
"""

import numpy as np
import matplotlib.pyplot as plt
import copy 
from matplotlib import _cntr as cntr

execfile('lib_physicalconstants.py')
execfile('class_vessel_threestage.py')


#function definitions
def try_AscentSimulation(vessel, velocities, pitches):
  #setup meshgrid and number of simulations
  [v_, p_] = np.meshgrid(velocities, pitches)  
  a_ = np.empty_like(v_)    #apoapsis [m]
  q_ = np.empty_like(v_)    #maxq [kg/ms2]
  l_ = np.empty_like(v_)    #delta-v losses combined drag and gravity [m/s]
  m_ = np.empty_like(v_)    #vessel endmass try to maximize
  b_ = np.empty_like(v_)    #total burntime of stage 3, should correlate with the vessel endmass
  
  #iterate through all velocities and pitches and store apoapsis and maxq
  it = np.nditer(v_, flags=['multi_index'])
  while not it.finished:
    ves = copy.deepcopy(vessel)             #initiate/create new vessel NOTE this creates a copy of vessel and leaves the original vessel untouched   
    v = it[0]
    p = p_[it.multi_index]
    
    #do the simulation with initial v and p
    #to a target altitude of 200 km apoapis and circular orbit
    [trajectory, PEGparameters] = ves.AscentSimulation(v, p,  R+ 200000,0,3)

    #store maxq value and apoapsis from getOrbital() function    
    if ves.state == 's3':    
      q_[it.multi_index] = ves.maxq
      orbit = ves._getOrbital()
      a_[it.multi_index] = orbit[0]
      l_[it.multi_index] = np.abs(ves.vgloss) + np.abs(ves.vdloss)
      m_[it.multi_index] = ves.m
      b_[it.multi_index] = ves.s3_totalburntime
    else:
      m_[it.multi_index] = ves.s2_m0 - ves.s2_mp
      b_[it.multi_index] = 0
    
    print it.multi_index   
    it.iternext()
  
  return v_, p_, a_, q_, l_, m_, b_

def writeKSbootscript(vessel, name):
  filename = 'boot_PEG_' + name + '.ks'
  f = open(filename,'w')
  
  #first write name in variable
  f.write('GLOBAL P_vName IS "' + name + '".\n\n') # python will convert \n to os.linesep
  
  #next write sequence
  f.write('//ignition (delay before release), booster jettison, fairings jettison, cutoff, separation, ullage, ignition, PEG activation, stage 2 maxT\n')
  f.write('GLOBAL P_seq IS LIST(' + ' -%.2f' % (1.4) + ', , , ' + '%.2f' % (vessel.s1_tb - 1.4) + ', , , , , ' + '%.2f' % (vessel.s2_tb) + ').//remember to fill this!\n')



  
  f.close()    # you can omit in most cases as the destructor will call it
  
  return


#define ships drag curfe
Cd = np.array([[ 0, 0.122], [256.0, 0.122], [343.2, 0.883], [643.5, 1.258], [909.5, 1.154], [1673,  0.676],[9999,  0.776]])

ship = VESSEL3S(
    #Saturn V, https://en.wikipedia.org/wiki/Saturn_V 
    #NOTE check masses of stages payload and subsequent stages have to be in the stage mass!
   s1_m0 = 3049200,             #stage 1: lauchmass [kg]
   s1_mp = 2160000,            #stage 1: propellant mass [kg]                
   s1_thrust_asl = 34020000,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl =  263,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 263,   #FOR NOW       #stage 1: specific impulse isp in vaccum [s]
   s1_A =  80.12,               #stage 1: reference cross sectional area [m2]
   s2_m0 = 759200,              #stage 2: mass at separation [kg]
   s2_mp =  456100,              #stage 2: propellant mass [kg]
   s2_thrust_asl = 4400000,         #stage 2: thrust in vaccum [N]
   s2_isp_asl = 421,   #FOR NOW    #stage 2: specific impulse in vaccum [s]
   s2_isp_vac = 421,
   s2_A = 80.12,                #stage 2: reference cross sectional area [m2]
   s3_m0 = 163000,
   s3_mp = 109500,
   s3_thrust = 1000000,
   s3_isp = 421,    #cannot be zero
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

ves = VESSEL3S(
    #sample 2stage rocket from PEG example code
   s1_m0 = 138855,             #stage 1: lauchmass [kg]
   s1_mp = 90603*0.96,            #stage 1: propellant mass [kg]                
   s1_thrust_asl = 1777537,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl =  252,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 290,          #stage 1: specific impulse isp in vaccum [s]
   s1_A =  (3.05/2)*np.pi,        #stage 1: reference cross sectional area [m2]
   s2_m0 = 48252,              #stage 2: mass at separation [kg]
   s2_mp =  31143*0.9,              #stage 2: propellant mass [kg]
   s2_thrust_asl =  269013,         #stage 2: thrust in vaccum [N]
   s2_isp_asl = 220,        #FOR NOW    #stage 2: specific impulse in vaccum [s]
   s2_isp_vac = 316,
   s2_A = (3.05/2)*np.pi,                #stage 2: reference cross sectional area [m2]
   s3_m0 = 17109,
   s3_mp = 14306,
   s3_thrust =  185000,
   s3_isp =  449,    #cannot be zero
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

#test a 2s rocket as 3s rocket
ves2s = VESSEL3S(
    #sample 2stage rocket from PEG example code
   s1_m0 = 97198,             #stage 1: lauchmass [kg]
   s1_mp = 75744,            #stage 1: propellant mass [kg]                
   s1_thrust_asl = 1217150,   #stage 1: thrust at sea level ASL [N]
   s1_isp_asl =  230,          #stage 1: specific impulse isp at sea level [s]
   s1_isp_vac = 250,          #stage 1: specific impulse isp in vaccum [s]
   s1_A =  7.06,              #stage 1: reference cross sectional area [m2]
   s2_m0 = 7442,              #stage 2: mass at separation [kg]
   s2_mp =  3323.0751881746983,   #stage 2: propellant mass [kg]
   s2_thrust_asl =  55400,     #stage 2: thrust in vaccum [N]
   s2_isp_asl = 340,          #FOR NOW    #stage 2: specific impulse in vaccum [s]
   s2_isp_vac = 340,
   s2_A = 2.0,                #stage 2: reference cross sectional area [m2]
   s3_m0 = 4118.924811825302,
   s3_mp =  2960.9248118253017,
   s3_thrust =  55400,
   s3_isp =  340,             #cannot be zero
   Cw = Cd,                   #drag curve [#] as function of velocity [m/s]
)

ves.add_launchsite(441,45.9)    #add launchsite at altitude of 441 m and latitude of 45.9 deg

ves.AscentSimulation(50, 86.5, R+ 200000, 0,3)

writeKSbootscript(ves, 'test')


try_velocities = np.linspace(50,50,1)
try_pitches = np.linspace(82,89,16)

#[V_,P_,A_,Q_, L_, M_, B_ ] = try_AscentSimulation(ves, try_velocities, try_pitches)


plt.figure()
plt.plot(ves.trajectory[:,0],ves.trajectory[:,1])
plt.xlabel('Time [s]')
plt.ylabel('Tangential velocity (orbit) [m/s]')

plt.figure()
plt.plot(ves.trajectory[:,0],ves.trajectory[:,2])
plt.xlabel('Time [s]')
plt.ylabel('Radial (vertical) velocity [m/s]')

plt.figure()
plt.plot(ves.trajectory[:,0],ves.trajectory[:,3])
plt.xlabel('Time [s]')
plt.ylabel('Altitude [m]')

plt.figure()
plt.plot(ves.trajectory[:,0],ves.trajectory[:,4])
plt.xlabel('Time [s]')
plt.ylabel('Pitch [deg]')

plt.figure()
plt.plot(ves.trajectory[:,0],ves.trajectory[:,6])
plt.ylabel('Dynamic pressure q [kg/ms2]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(ves.PEGparameters[:,0],ves.PEGparameters[:,1])
plt.ylabel('PEG parameter A [#]')
plt.xlabel('Time [s]')

plt.figure()
plt.plot(ves.PEGparameters[:,0],ves.PEGparameters[:,2])
plt.ylabel('PEG parameter B [#]')
plt.xlabel('Time [s]')

