# -*- coding: utf-8 -*-
"""
powered_explicit_guidance.py
Created on Mon Feb 22 13:47:41 2016

Copy and recreation of matlab code from https://github.com/Noiredd/PEGAS

Implementation of Powered Explicit Guidance major loop.
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660006073.pdf
http://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance

Implementation of 2 Stage PEG algorithm 
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660006073.pdf
http://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance

Dependencies: lib_physicalconstants.py
              poweredExplicitGuidance.py

@author: sig
"""

import numpy as np
from decimal import *  

def poweredExplicitGuidance3Stage(dT, alt, vt, vr, tgt, acc, a2_0, ve1, ve2, omega_T1,  oldA1, oldB1, oldT1, oldT2):
  #dT     - time since last itration [s] can be zero
  #oldT1  - burntime of first stage left [s], this value is always known, since we assume that second stage (stage 1) burns out completely
  
  #alt    - current altitude as distance from body center [m]
  #vt     - tangential velocity (horizontal - orbital) [m/s]
  #vr     - radial velocity (veritcal - away from the Earth) [m/s]
  #tgt    - target altitude (measured as current) [m]
  #acc    - current vehicle acceleration [m/s^2]
  #ve     - effective exhaust velocity (isp*g0) [m/s] 
  
  global mu     #requires to run lib_physicalconstants.py first
  
  #NOTE staging not in here
  tau = ve1 / acc      #A sort of normalized mass - time to burn the vehicle completely as if it were all propellant  [s]
    
  #update time first
  B1 = oldB1
  A1 = oldA1 + oldB1*dT
  T1 = oldT1     #time update of T1 is done outside this function
  T2 = oldT2
  T = oldT1 + oldT2 #total burn time left
    
  #42a and 42b 
  rdot_0 = vr             #radial velocity [m/s] right now
  r_0 = alt               #radial distance [m] right now
  b01_T1 = -ve1*np.log(1-T1/tau)
  b11_T1 = b01_T1*tau - ve1*T1
  b21_T1 = b11_T1*tau - ve1*T1*T1/2
  c01_T1 = b01_T1*T1 - b11_T1 
  c11_T1 = c01_T1*tau - ve1*T1*T1/2
    
  #eq. 6 for first stage for calculating eq 34. and ultimately for eq. 40 (dA and dB)
  rdot_T1 = rdot_0 + b01_T1*A1 + b11_T1*B1
  r_T1 = r_0 + rdot_0*T1 + c01_T1*A1 + c11_T1*B1
  
  #current required pitch 
  a1_0 = acc        #current acceleration or accceleration of stage 1 at time zero (now)  
  #C1_0 =  (mu/(Decimal(alt)*Decimal(alt)) - Decimal(vt*vt)/Decimal(alt)) / Decimal(a1_0)     #Portion of vehicle acceleration used to counteract gravity and centrifugal force 
  #C1_0 = float(str(C1_0))  
  C1_0 = (mu/(alt*alt) - vt*vt/alt)  / a1_0
  fr1_0 = A1 + C1_0 
    
  #Horizontal state at staging (delta time T1 in the future)
  #for loop to test convergence of omega_T1, only one step per PEG step needed I think
  a1_T1 = acc / (1 - T1/tau)
  #C1_T1 = (mu/(Decimal(r_T1)*Decimal(r_T1)) - Decimal(omega_T1*omega_T1)*Decimal(r_T1)) / Decimal(a1_T1)
  #C1_T1 = float(str(C1_T1))
  C1_T1 = (mu/(r_T1*r_T1) - omega_T1*omega_T1*r_T1) / a1_T1
  fr1_T1 = A1 + B1*T1 +  C1_T1
  fr1dot_T1 = (fr1_T1 - fr1_0)/T1 
  ftheta1_T1 = 1 - fr1_T1*fr1_T1/2      
  fthetadot1_T1 = - (fr1_T1*fr1dot_T1)
  fthetadotdot1_T1 = - fr1dot_T1*fr1dot_T1/2
  
  #42c equation (not sure if I have to use (r(0)- r(T) or r(0) - r(T1))
  dh_T1 = (r_0 + tgt)/2 * (ftheta1_T1*b01_T1 + fthetadot1_T1*b11_T1 + fthetadotdot1_T1*b21_T1)
  h_0 = np.linalg.norm( np.cross([0, 0, alt],[vt, 0, vr]) )    #current angular momentum [unit? m2/s]
  h_T1 = h_0 + dh_T1    
  vtheta_T1 = h_T1/r_T1       #tangential velocitiy [m/s] at time T1 (staging)
  omega_T1 = vtheta_T1 / r_T1 #estimate of the angular velocity at staging (should iteratively converge)
 
  #guidance staging discontinuities
  #a2_0 acceleration at staging event?? not sure what this is hypotheticaly if we would stage now what would be the acceleration?
  dA = C1_T1 * a1_T1 * (1/a1_T1 - 1/a2_0)
  #dB = (Decimal(-C1_T1) * Decimal(a1_T1) * (1/Decimal(ve1) - 1/Decimal(ve2))) + (3*Decimal(omega_T1*omega_T1) - 2*mu/Decimal(r_T1*r_T1*r_T1))*Decimal(rdot_T1)*(1/Decimal(a1_T1) - 1/Decimal(a2_0))
  #dB = float(str(dB))
  dB = -C1_T1*a1_T1*(1/ve1 - 1/ve2) + (3*omega_T1*omega_T1 - 2*mu/(r_T1*r_T1*r_T1))*rdot_T1*(1/a1_T1 - 1/a2_0)  
  
  #Solve Explicit Guidance Equations. An estimated stage 2 burn time T_2 is needed 
  tau2 = ve2/a2_0

  b02_T2 = -ve2*np.log(1-T2/tau2)
  b12_T2 = b02_T2*tau2 - ve2*T2
  c02_T2 = b02_T2*T2 - b12_T2
  c12_T2 = c02_T2*tau2 - ve2*T2*T2/2

  MA11 = b01_T1 + b02_T2  
  MA12 = b11_T1 + b12_T2 + b02_T2*T1
  MA21 = c01_T1 + c02_T2 + b01_T1*T2
  MA22 = c11_T1 + b11_T1*T2 + c02_T2*T1 + c12_T2
  MA = np.array([[MA11, MA12],[MA21, MA22]])
  
  rdot_T = 0    #target radial velocity is zero for circular orbit   [m/s]
  
  r_T = tgt     #target altitude (radial distance) [m]  
  MB1 = rdot_T - rdot_0 - b02_T2*dA - b12_T2*dB
  MB2 = r_T - r_0 - rdot_0*T - c02_T2*dA - c12_T2*dB
  MB = np.array([MB1, MB2])
  MX = np.linalg.solve(MA,MB)
  A1 = MX[0]
  B1 = MX[1] 
  
  #Estimate stage 2 burn time T2, using the same method as for the single stage, but from staging state to target state 
  A2 = dA + A1 + B1*T1
  B2 = dB + B1
   
  v_tgt = np.sqrt(mu/tgt) 
  h_T =  v_tgt*tgt    #ht_vec = np.cross([0, 0, tgt],[v_tgt, 0, 0])  and ht = np.linalg.norm(ht_vec)
  dh = h_T- h_T1
  rbar = (r_T1 + tgt)/2 
  #C2 = (mu/(Decimal(r_T1)*Decimal(r_T1)) - Decimal(omega_T1*omega_T1)*Decimal(r_T1)) / Decimal(a2_0)
  #C2 = float(str(C2))
  C2 = (mu/(r_T1*r_T1) - omega_T1*omega_T1*r_T1) / a2_0 
  fr2_0 = A2 + C2
  omega_T = v_tgt/tgt
  a2_T2 = a2_0 / (1 - T2/tau2)
  #C2_T = (mu/(Decimal(r_T)*Decimal(r_T)) - Decimal(omega_T*omega_T)*Decimal(r_T)) / Decimal(a2_T2)
  #C2_T = float(str(C2_T))
  C2_T = (mu/(r_T*r_T) - omega_T*omega_T*r_T ) / a2_T2 
  fr2_T = A2 + B2*T2 + C2_T
  fr2dot =  (fr2_T - fr2_0)/T2
  ftheta2 = 1 - fr2_T*fr2_T/2
  ftheta2dot = - (fr2_T*fr2dot)
  ftheta2dotdot = fr2dot*fr2dot/2
  
  dv = (dh/rbar + ve2*T2*(ftheta2dot + ftheta2dotdot*tau2) + ftheta2dotdot*ve2*T2*T2/2) / (ftheta2 + ftheta2dot*tau2 + ftheta2dotdot*tau2*tau2 )
  T2 = tau2 * (1 - np.exp(-dv/ve2))    
    
  return A1, B1, C1_0, T1, omega_T1, T2


