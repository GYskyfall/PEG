# -*- coding: utf-8 -*-
"""
powered_explicit_guidance.py
Created on Mon Feb 22 13:47:41 2016

Copy and recreation of matlab code from https://github.com/Noiredd/PEGAS

Implementation of Powered Explicit Guidance major loop.
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660006073.pdf
http://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance

Dependencies: lib_physicalconstants.py

@author: sig
"""

import numpy as np
from decimal import *  

def poweredExplicitGuidance(cycletime, alt, vt, vr, tgt,   acc, ve,   oldA, oldB, oldT):
  #dt     - length of the major computer cycle (time between PEG calculations) [s]
  #alt    - current altitude as distance from body center [m]
  #vt     - tangential velocity (horizontal - orbital) [m/s]
  #vr     - radial velocity (veritcal - away from the Earth) [m/s]
  #tgt    - target altitude (measured as current) [m]
  #acc    - current vehicle acceleration [m/s^2]
  #ve     - effective exhaust velocity (isp*g0) [m/s]
  #basing on current state vector, vehicle params and given ABT estimates
  #new A and B, and calculates new T from them
  #also outputs a C component for the purposes of guidance
  #passing A=B=0 causes estimation of those directly from the given oldT
  
  global mu     #requires to run lib_physicalconstants.py first
    
  #if this is the first time estimate of oldA, oldB and oldT 
  #use a reasonable oldT, for instance burntime of stage 2
  #calculate A and B 
  tau = ve / acc      #A sort of normalized mass - time to burn the vehicle completely as if it were all propellant  [s]
  tgt = float(tgt)
  
  if (oldA == 0) and (oldB == 0):
    if oldT > tau:    #this is to prevent NAN from logarithm because of bad (to high) estimate of T
      oldT = 0.9*tau
    b0 = -ve*np.log(1-oldT/tau)
    b1 = b0*tau - ve*oldT
    c0 = b0*oldT - b1
    c1 = c0*tau - ve*oldT*oldT/2
    MB = np.array([-vr, tgt - alt - vr*oldT])  #NOTE for now r_dot target is zero for circular target orbit  
    MA = np.array([[b0, b1],[c0, c1]])
    MX = np.linalg.solve(MA,MB)
    A = MX[0]
    B = MX[1]    
  else: 
    A = oldA
    B = oldB
    
  #navigation vectors rotating coordinate system the z-axis is always parallel to radial direction
  #and x-axis always parallel to tangential direction
  h_vec = np.cross([0, 0, alt],[vt, 0, vr])     #current angular momentum [unit? m2/s]
  h = np.linalg.norm(h_vec)
  v_tgt = np.sqrt(mu/tgt)                       #orbital velocity at a target altitude of tgt [m/s]
  ht_vec = np.cross([0, 0, tgt],[v_tgt, 0, 0])  #zero vertical velocity target
  ht = np.linalg.norm(ht_vec)
  dh = ht - h                                   #angular momentum to gain [unit? m2/s]
  rbar = (alt + tgt)/2                          #mean radius [m]
  #omega = vt/alt
  #print omega
  
 
  #Vehicle performence
  #C =  (mu/(Decimal(alt)*Decimal(alt)) - Decimal(vt*vt)/Decimal(alt)) / Decimal(acc)     #Portion of vehicle acceleration used to counteract gravity and centrifugal force 
  #C = float(str(C))   #convert C back to floating point number  
  C = (mu/(alt*alt) - vt*vt/alt) / acc 
  fr = A + C                                    #sin(pitch) at current time 
   
  #estimation
  #CT = (mu/(Decimal(tgt)*Decimal(tgt)) - Decimal(v_tgt)*Decimal(v_tgt)/Decimal(tgt)) / (Decimal(acc) / (1 - Decimal(oldT)/Decimal(tau)))
  #CT = float(str(CT)) 
  CT = (mu/(tgt*tgt) - v_tgt*v_tgt/tgt) / (acc / (1 - oldT/tau))
  frT = A + B*oldT + CT                         #sin(pitch) at burnout 
  frdot = (frT - fr)/oldT
  ftheta = 1 - fr*fr/2
  fthetadot = -(fr*frdot)
  fthetadotdot = -frdot*frdot/2

  #Ideal velocity-to-gain in order to get target angular momentum, based on current pitch guidance
  dv = (dh/rbar + ve*(oldT-cycletime)*(fthetadot + fthetadotdot*tau) + fthetadotdot*ve*(oldT-cycletime)*(oldT-cycletime)/2) / (ftheta + fthetadot*tau + fthetadotdot*tau*tau) 
  
  # Estimate updated burnout time
  T = tau * (1 - np.exp(-dv/ve))

  if T >= 7.5:
    b0 = -ve*np.log(1-T/tau)
    b1 = b0*tau - ve*T
    c0 = b0*T - b1
    c1 = c0*tau - ve*T*T/2
    MB = np.array([-vr, tgt - alt - vr*T])  #NOTE for now r_dot target is zero for circular target orbit  
    MA = np.array([[b0, b1],[c0, c1]])
    MX = np.linalg.solve(MA,MB)
    A = MX[0]
    B = MX[1] 
  else:
    A = oldA
    B = oldB
    
  return A, B, C, T


