# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 15:09:20 2016

@author: sig

Set up a "standard" two stage rocket

Recreation of matlab code from https://github.com/Noiredd/PEGAS
https://smartech.gatech.edu/bitstream/handle/1853/6820/dukeman_greg_a_200505_phd.pdf
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660006073.pdf
http://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance


"""

import numpy as np

class VESSEL:
  def __init__(self,
               s1_m0,           #stage 1: lauchmass [kg]
               s1_mp,           #stage 1: propellant mass [kg]                
               s1_thrust_asl,   #stage 1: thrust at sea level ASL [N]
               s1_isp_asl,      #stage 1: specific impulse isp at sea level [s]
               s1_isp_vac,      #stage 1: specific impulse isp in vaccum [s]
               s1_A,            #stage 1: reference cross sectional area [m2]
               s1_engT,        #stage 1: waittime until launchtower release engine preingition time [s]     
               s2_m0,           #stage 2: mass at separation [kg]
               s2_mp,           #stage 2: propellant mass [kg]
               s2_thrust,       #stage 2: thrust in vaccum [N]
               s2_isp,          #stage 2: specific impulse in vaccum [s]
               s2_A,            #stage 2: reference cross sectional area [m2]
               Cw,              #drag curve [#] as function of velocity [m/s]
               ):                
    global g0 
    global R
    global mu
           
    #vessel info stage 1
    self.s1_m0 = s1_m0
    self.s1_mp = s1_mp
    self.s1_thrust_asl = s1_thrust_asl 
    self.s1_isp_asl = s1_isp_asl
    self.s1_isp_vac = s1_isp_vac
    self.s1_A = s1_A
    self.s1_mfr = s1_thrust_asl/(s1_isp_asl*g0)                 #stage 1: mass flow rate (at sea level) [kg/s]
    self.s1_tb = s1_mp/self.s1_mfr                              #stage 1: burn time [s]
    self.s1_twr_asl = s1_thrust_asl/(s1_m0*g0)                  #stage 1: thrust weight ratio twr (asl) [#]
    self.s1_twr_final = self.s1_mfr*s1_isp_vac/(s1_m0-s1_mp)    #stage 1: final thrust weight ratio [#]
    self.Cw = Cw
    self.s1_engT = s1_engT
    
    #vessel info stage 2
    self.s2_m0 = s2_m0
    self.s2_mp = s2_mp
    self.s2_thrust = s2_thrust
    self.s2_isp = s2_isp
    self.s2_A = s2_A
    self.s2_mfr = s2_thrust/(s2_isp*g0)                         #stage 2: mass flow rate (at sea level) [kg/s]
    self.s2_tb = s2_mp/self.s2_mfr                              #stage 2: burn time [s]
    self.s2_twr_0 = s2_thrust/(s2_m0*g0)                        #stage 2: thrust weight ratio twr [#]
    self.s2_twr_final = self.s2_mfr*s2_isp/(s2_m0-s2_mp)        #stage 2: final thrust weight ratio [#]
 
    #initialize class
    self._init()
    self._setstate('s1')
    
  
  def _init(self):
    #vessel geometry paramters
    self.alt = 0                #initial altitude above sea level at launch site [m]
    self.v = [0, 0, 0]          #initial velocity through air [m/s]
    self.pitch = 90             #initial pitch = 90 [deg] (UP)
    self.yaw = 0
    self.roll = 0   
    self.lat = 0                #initial latitude = 0  [deg] (equator)
    self.R = R                  #initial radial distance to earth center [m]
    self.vx_gain = (2*np.pi*self.R/24/3600)  #horizontal velocitiy gain from earths rotation at launchsite [m/s]
    self.q = 0
    self.t = 0                  #total time elapsed (vessel lifetime) [s]

    #losses in delta-v
    self.vgloss = 0
    self.vdloss = 0
    self.maxq = 0
    
    #logging 
    self.trajectory =  np.zeros((1,7))    #log launch trajectory
    self.PEGparameters = np.zeros((1,5)) 
    
  def _setstate(self, state):
    if state == ('s1'): 
      alt = self.alt                                          #vessel altitude 
      self.m = self.s1_m0
      self.isp = self.s1_isp_asl*atmp(alt)/atmp(0)   + self.s1_isp_vac*(1-atmp(alt)/atmp(0) )           #calculate actual isp from interpolation between sea level isp and vacuum isp from air pressure
      self.mfr = self.s1_mfr
      self.thrust = self.isp*g0*self.mfr
      self.twr = self.thrust/(self.m*g0)
      self.A = self.s1_A
      self.state = 's1'
    elif state == 's2':
      alt = self.alt                                          #vessel altitude 
      self.m = self.s2_m0
      self.isp = self.s2_isp           
      self.mfr = self.s2_mfr
      self.thrust = self.isp*g0*self.mfr
      self.twr = self.thrust/(self.m*g0)
      self.A = self.s2_A
      self.state = 's2'
    return

  def add_launchsite(self, altitude, latitude):
    self.alt = altitude
    self.lat = latitude
    self.R = R + altitude
    self.vx_gain = (2*np.pi*self.R/24/3600)*np.cos(latitude*2*np.pi/360)   #horizontal velocity gained from Earth's rotation [m/s]
    
  def logtrajectory(self):
    self.q = 0.5*(self.v[0]*self.v[0]+self.v[1]*self.v[1]+self.v[2]*self.v[2])*rho(self.alt)                   #dynamic pressure  [kg/ms2]
    self.trajectory = np.append(self.trajectory,[[self.t, self.v[0], self.v[2], self.alt, self.pitch, self.m, self.q]],axis=0)
      
  def logPEGparameters(self):
    self.PEGparameters = np.append(self.PEGparameters,[[self.t, self.PEG_A, self.PEG_B, self.PEG_C, self.PEG_T]],axis=0)

  #function calculates apoapsis and periapsis for a coasting vessel (no acceleration) on earth
  def _getOrbital(self):
    #function from https://github.com/Noiredd/PEGAS
    #at any given point on the ellipse: v = sqrt( mu*( 2/r-1/a ) )    
    v_orb = self.v + [self.vx_gain, 0, 0]   #orbital velocity [m/s]
    #a - semimajor axis: 1/a = 2/r - v^2/mu
    sma = 1/( 2/self.R - np.power(np.linalg.norm(v_orb),2)/mu );
    #specific orbital energy
    soe = -mu/(2*sma);
    #from orbital eccentricity equations: e = sqrt( 1 + 2*soe*h^2/mu^2 )
    #srh - specific relative angular momentum: h = r*v
    srh = np.linalg.norm(np.cross([self.R, 0, 0], [v_orb[2], v_orb[0], 0]));
    Ecc = np.sqrt(1 + 2*soe*np.power(srh,2) / (mu*mu));
    #apsi distances are given: r_a/p = a*( 1 +/- e )
    Ap = sma*(1+Ecc)    #raw value
    Ap = (Ap-R)         #return in [m]
    Pe = sma*(1-Ecc)
    Pe = (Pe-R)
    
    return  [Ap, Pe, Ecc]
  
  def _propagate_lock(self,dt, pitch, yaw, roll)     :
    self.pitch = pitch  
    self.yaw = yaw
    self.roll = roll
    self._propagate_step(dt)
    self.pitch = pitch  #correct that rocket is steering 
    
    return 
      
  def _propagate_step(self,dt):
    #vessel mass changes
    m = self.m
    dm = dt*self.mfr
    self.m = m - dm    
    
    #update isp at altitude 
    alt = self.alt                                        #current vessel altitude over ground z = 0                                            
    if self.state == 's1':    
      self.isp = self.s1_isp_asl*atmp(alt)/atmp(0)   + self.s1_isp_vac*(1-atmp(alt)/atmp(0) )    
        
    #calculate drag force on vessel (angle of attack zero)
    q = 0.5*(self.v[0]*self.v[0]+self.v[1]*self.v[1]+self.v[2]*self.v[2])*rho(alt)                   #dynamic pressure  [kg/ms2]
    Cd =  np.interp(np.linalg.norm(self.v), self.Cw[:,0], self.Cw[:,1])   #current vessel drag coefficient     
    Fd = -Cd*self.A*q                                                 #drag force along vessel propagation
                               
    
    #calculate thrust along vessel propagation
    Ft = self.isp*g0*self.mfr
    
    #calculate gravity force along z-direction
    Fg = -self.m*mu/(self.R*self.R)
    
    #calculate centrifugal force 
    Fc = self.m*self.v[0]*self.v[0]/self.R
   
    #total force vector  //YAW = 0 FOR NOW! -> only propagation in x-z plane
    F = np.array([(Ft+Fd)*np.cos(self.pitch*2*np.pi/360), 0 , (Ft+Fd)*np.sin(self.pitch*2*np.pi/360) + Fg + Fc])
        
    #acceleration, velocity and position
    acc = F/m
    vel = self.v + acc*dt
    alt = self.alt + vel[2]*dt
    self.v = vel
    self.alt = alt   
    self.R = R + alt
    self.pitch = np.arctan2(vel[2],vel[0])*360/(2*np.pi)   
     
    #calculate losses in delta-v
    self.vgloss = self.vgloss + Fg/m*dt  
    self.vdloss = self.vdloss + Fd/m*dt
    
    #update time
    self.t = self.t + dt
     
    return
    
  def _propagatePEG_step(self, dt, cycletime, target_altitude, target_vz):
  
    self.PEG_alt = R + self.alt                                                         #set PEG altitude to current radial distance from center of earth [m]
    self.PEG_vt = self.v[0] + self.vx_gain      #current tangential velocity [m/s]
    self.PEG_vr = self.v[2]                                                             #current radial velocity (in z-direction) [m/s]
    self.PEG_ve = self.isp*g0                   #effective exhaust velocitz [m/s]
    self.PEG_a = self.PEG_ve*self.mfr/self.m    #current  acceleration [m/s2]
    
    #to safe computation time run major loop only every cyclretime seconds
    #to know the time passed, this function stores the burnout time T to self.PEG_lastcall [s]
    if  (self.PEG_lastcall - self.PEG_T) >= cycletime:
      [A, B, C, T] = poweredExplicitGuidance( cycletime, self.PEG_alt, self.PEG_vt, self.PEG_vr, target_altitude, self.PEG_a, self.PEG_ve, self.PEG_A, self.PEG_B, self.PEG_T);
      self.PEG_lastcall = T
      self.PEG_A = A
      self.PEG_B = B
      self.PEG_T = T  
      self.PEG_C = C
    else:
      self.PEG_T = self.PEG_T - dt
    
    #calculate pitch
    sinpitch = self.PEG_A - (self.PEG_lastcall - self.PEG_T)*self.PEG_B + self.PEG_C     
    pitch = np.arcsin(min(1, max(-1, sinpitch)))/(2*np.pi)*360
    self.pitch = pitch    #point vessel along pitch direction    
    
    #perform rocket physics step 
    #vessel mass changes
    m = self.m
    dm = dt*self.mfr
    self.m = m - dm    
    
    #calculate gravity and centrifugal force
    ca = (self.v[0] + self.vx_gain)*(self.v[0] + self.vx_gain) / self.R
    ga = mu / (self.R*self.R)
    vg = (ca - ga)*dt
    
    self.v[0] = self.v[0] + self.PEG_a*np.cos(self.pitch*2*np.pi/360)*dt
    self.v[2] = self.v[2] + vg + self.PEG_a*np.sin(self.pitch*2*np.pi/360)*dt
    self.alt = self.alt + self.v[2]*dt
    self.R = self.alt + R    
    
    #update time
    self.t = self.t + dt
    
    return
    
    
  #function propagates vessel to simulate launch with initial steering of pitch_init [deg] when velocitiy is greater than 
  #v_start [m/s]. After that the vessel launches along its velocity vector (zero angle of attack)
  #returns trajectory [time (s), velocity vector (m/s), altitude (m), pitch (deg), ridial position (m), status (1,2,3), dynamic pressure q (kg/ms2)]
  def AscentSimulation(self, v_start, pitch_init, target_altitude, target_vz, stage):
 
    #start simulation
    dt = 0.1                    #time step [s]     
    self.t = 0                 #simulation time [s]
    
    #first wait 1.5 s to launch after ignition of engines
    while self.t < self.s1_engT:
      self.m = self.m - self.mfr * dt
      self.t = self.t + dt
      
      #log launch trajectory data      
      self.logtrajectory()
    
    #perform simulation    
    while self.t < self.s1_tb:
      #check status of flight      
      if (np.linalg.norm(self.v) < v_start) :
        #lock steering to UP until velocity > v_init
        self._propagate_lock(dt, 90, 0 , 0)
      elif (np.linalg.norm(self.v) >= v_start and np.arctan2(self.v[2],self.v[0])*360/(2*np.pi) > pitch_init):
        #pitchover maneuver
        self._propagate_lock(dt, max(self.pitch-dt, pitch_init), 0 , 0)    # 1 deg/s change in pitch 
      elif (np.linalg.norm(self.v) >= v_start and np.arctan2(self.v[2],self.v[0])*360/(2*np.pi) <= pitch_init):
        #free gravity turn along prograde vector
        self._propagate_step(dt)
      else:
        break
      
      #break ascent simulation if altitude is negative or pitch is negative
      if (self.pitch < 0) or (self.alt < 0):
        return 
      
      #log launch trajectory data      
      self.logtrajectory()
      
      
    #update time elapsed
    self.s1_totalburntime = self.t
 
    #jettison stage 1 
    self._setstate('s2')
    if stage == 1:        
      return 
      
    
    #get to circular orbit with PEG
    #initialize PEG to orbit
    self.PEG_state = 1
    self.PEG_T = self.s2_tb
    self.PEG_lastcall = self.PEG_T
    self.PEG_A = 0
    self.PEG_B = 0

    #first PEG step with dt = 0    
    self.PEG_alt = R + self.alt                                                         #set PEG altitude to current radial distance from center of earth [m]
    self.PEG_vt = self.v[0] + (2*np.pi*R/24/3600)*np.cos(self.lat*2*np.pi/360)     #current tangential velocity [m/s]
    self.PEG_vr = self.v[2]                                                             #current radial velocity (in z-direction) [m/s]
    self.PEG_ve = self.isp*g0                   #effective exhaust velocitz [m/s]
    self.PEG_a = self.PEG_ve*self.mfr/self.m    #current  acceleration [m/s2]
      
    #print self.PEG_alt, self.PEG_vt, self.PEG_vr,  self.PEG_a, self.PEG_ve, self.PEG_T    
    [A, B, C, T] = poweredExplicitGuidance( 0, self.PEG_alt, self.PEG_vt, self.PEG_vr, target_altitude, self.PEG_a, self.PEG_ve, self.PEG_A, self.PEG_B, self.PEG_T);
    self.PEG_A = A
    self.PEG_B = B
    self.PEG_C = C
    self.PEG_T = T
    self.PEG_lastcall = T
 
    #new pitch angle
    sinpitch = A + C      #-dt*B term is always zero, since PEG is just freshly calculated
    pitch = np.arcsin(min(1, max(-1, sinpitch)))/(2*np.pi)*360
     
    #coast for ullage and pitch change, requires RCS capable second stage    
    pitchdot = 1     #0.5 deg change per second this is only an estimate
    
    maxpitchovertime = self.t + 15
    self.mfr = 0    
    while self.pitch != pitch:
      if np.abs(self.pitch - pitch) < pitchdot*dt:
        break
      
      if (self.t - maxpitchovertime) > 15:
        break
      
      self.pitch = self.pitch + pitchdot*np.sign(pitch-self.pitch)*dt
      self._propagate_lock(dt, self.pitch, 0 , 0)
           
      self.t = self.t + dt
      
      #log launch trajectory data      
      self.logtrajectory()
      
    #update time and ignite state 2
    self.mfr = self.s2_mfr    
    
    #ignition stage 2
    #continue with ascent guided PEG
    timeats2ignition = self.t
    timeats2burnout = self.t + self.s2_tb    
    while self.t < (timeats2burnout):   #burnout of stage 2
          
      if self.PEG_T <= 0.0: 
        break
      
      self._propagatePEG_step(dt, 2, target_altitude, target_vz)
     
      #log launch trajectory data      
      self.logtrajectory()
      
      #log PEG parameters
      self.logPEGparameters()
   
   
    return 
