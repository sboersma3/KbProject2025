
from sys import path
path.append(r"lib")

import lib as lib 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pdb #pdb.set_trace()
    

# model parameters
p          = lib.DefineParameters()

# model characteristics
ops        = {}
ops['np']  = len(p)                                 # #of model parameters
ops['nx']  = 6                                      # #of states
ops['ny']  = 13                                     # #of outputs
ops['nu']  = 3                                      # #controllable inputs
ops['nd']  = 4                                      # #disturbance inputs

# simulation parameters
c         = 86400
nDays     = 1                                       # #days in simulation
ops['h']  = 0.01                                   # sample period in hours
ops['t']  = np.arange(0,24*nDays+ops['h'],ops['h']) # initial time vector
ops['N']  = len(ops['t'])                           # number of samples in initial time vector
ops['N0'] = 0                                       # initial sample of the weather data

# signals
x         = np.zeros((ops['nx'],ops['N']+1))        # state
y         = np.zeros((ops['ny'],ops['N']))          # measurement
x[:,0]    = np.array([0,-5,0,0,0,0])                # initial state
u         = np.zeros((ops['nu'],ops['N']+1))        # control signal
u[:,0]    = np.array([0.24,0,0])                    # initial control signal
d         = lib.LoadDisturbances(ops)               # disturbance

        
for kk in range(ops['N']): # loop over time kk

           
    # propagate model one time step ahead
    x[:,kk+1]  = lib.fd(x[:,kk], u[:,kk], d[:,kk], p, ops['h'])
    
    # measure output from model
    y[:,kk]    = lib.g(x[:,kk], u[:,kk], d[:,kk], p)
       
    # control for the next time step
    u[:,kk+1]  = lib.controller(x[:,kk], u[:,kk], d[:,kk], p, ops)


###############################################################################
## plot
###############################################################################

# plot disturbance and controls
fig, ax = plt.subplots(ops['nu']+ops['nd'])
fig.set_size_inches(10.5, 10.5)

ax[0].plot(ops['t'][0:kk],u[0,0:kk], drawstyle='steps')
ax[0].grid(True)
ax[0].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[0].set_ylabel(r'$u_1$ (-)', fontsize=16)
#ax[0].set_ylim(bottom=0.9*min(u[0,1:kk]), top=1.1*max(u[0,1:kk]))

ax[1].plot(ops['t'][0:kk],u[1,0:kk], drawstyle='steps')
ax[1].grid(True)
ax[1].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[1].set_ylabel(r'$u_2$ (-)', fontsize=16)
#ax[1].set_ylim(bottom=0.9*min(u[1,1:kk]), top=1.1*max(u[1,1:kk]))

ax[2].plot(ops['t'][0:kk],u[2,0:kk], drawstyle='steps')
ax[2].grid(True)
ax[2].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[2].set_ylabel(r'$u_3$ (-)', fontsize=16)
#ax[2].set_ylim(bottom=0.9*min(u[2,1:kk]), top=1.1*max(u[2,1:kk]))

ax[3].plot(ops['t'][0:kk],d[0,0:kk], drawstyle='steps')
ax[3].grid(True)
ax[3].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[3].set_ylabel(r'$d_1$ (gC)', fontsize=16)
#ax[3].set_ylim(bottom=0.9*min(d[0,1:kk]), top=1.1*max(d[0,1:kk]))

ax[4].plot(ops['t'][0:kk],d[1,0:kk], drawstyle='steps')
ax[4].grid(True)
ax[4].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[4].set_ylabel(r'$d_2$ (kW)', fontsize=16)
#ax[4].set_ylim(bottom=0.9*min(d[1,1:kk]), top=1.1*max(d[1,1:kk]))

ax[5].plot(ops['t'][0:kk],d[2,0:kk], drawstyle='steps')
ax[5].grid(True)
ax[5].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[5].set_ylabel(r'$d_3$ (€/kWh)', fontsize=16)
#ax[5].set_ylim(bottom=0.9*min(d[2,1:kk]), top=1.1*max(d[2,1:kk]))

ax[6].plot(ops['t'][0:kk],d[3,0:kk], drawstyle='steps')
ax[6].grid(True)
ax[6].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[6].set_ylabel(r'$d_4$ (€/kWh)', fontsize=16)
#ax[6].set_ylim(bottom=0.9*min(d[3,1:kk]), top=1.1*max(d[3,1:kk]))

ax[6].set_xlabel(r'Time (hours)',fontsize=16)


# plot states
fig, ax = plt.subplots(ops['nx'])
fig.set_size_inches(10.5, 10.5)
ax[0].plot(ops['t'][0:kk],x[0,0:kk], drawstyle='steps')
ax[0].grid(True)
ax[0].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[0].set_ylabel(r'$x_1$ (kJ)', fontsize=16)
#ax[0].set_ylim(bottom=0.9*min(x[0,1:kk]), top=1.1*max(x[0,1:kk]))

ax[1].plot(ops['t'][0:kk],x[1,0:kk], drawstyle='steps')
ax[1].grid(True)
ax[1].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[1].set_ylabel(r'$x_2$ (gC)', fontsize=16)
#ax[1].set_ylim(bottom=0.9*min(x[1,1:kk]), top=1.1*max(x[1,1:kk]))

ax[2].plot(ops['t'][0:kk],x[2,0:kk], drawstyle='steps')
ax[2].grid(True)
ax[2].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[2].set_ylabel(r'$x_3$ (€)', fontsize=16)
#ax[2].set_ylim(bottom=0.9*min(x[2,1:kk]), top=1.1*max(x[2,1:kk]))

ax[3].plot(ops['t'][0:kk],x[3,0:kk], drawstyle='steps')
ax[3].grid(True)
ax[3].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[3].set_ylabel(r'$x_4$ (kW^pj(5).h)', fontsize=16)
#ax[3].set_ylim(bottom=0.9*min(x[3,1:kk]), top=1.1*max(x[3,1:kk]))

ax[4].plot(ops['t'][0:kk],x[4,0:kk], drawstyle='steps')
ax[4].grid(True)
ax[4].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[4].set_ylabel(r'$x_5$ (kW^pj(6).h)', fontsize=16)
#ax[4].set_ylim(bottom=0.9*min(x[4,1:kk]), top=1.1*max(x[4,1:kk]))

ax[5].plot(ops['t'][0:kk],x[5,0:kk], drawstyle='steps')
ax[5].grid(True)
ax[5].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[5].set_ylabel(r'$x_6$ (kW^pj(7))', fontsize=16)
#ax[5].set_ylim(bottom=0.9*min(x[5,1:kk]), top=1.1*max(x[5,1:kk]))

ax[5].set_xlabel(r'Time (hours)',fontsize=16)


# plot other outputs 
fig, ax = plt.subplots(ops['ny']-ops['nx'])
fig.set_size_inches(10.5, 10.5)
ax[0].plot(ops['t'][0:kk],y[6,0:kk], drawstyle='steps')
ax[0].grid(True)
ax[0].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[0].set_ylabel(r'$y_7$ (gC)', fontsize=16)

ax[1].plot(ops['t'][0:kk],y[7,0:kk], drawstyle='steps')
ax[1].grid(True)
ax[1].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[1].set_ylabel(r'$y_8$ (m)', fontsize=16)

ax[2].plot(ops['t'][0:kk],y[8,0:kk], drawstyle='steps')
ax[2].grid(True)
ax[2].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[2].set_ylabel(r'$y_9$ (-)', fontsize=16)

ax[3].plot(ops['t'][0:kk],y[9,0:kk], drawstyle='steps')
ax[3].grid(True)
ax[3].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[3].set_ylabel(r'$y_{10}$ (kW)', fontsize=16)

ax[4].plot(ops['t'][0:kk],y[10,0:kk], drawstyle='steps')
ax[4].grid(True)
ax[4].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[4].set_ylabel(r'$y_{11}$ (kW)', fontsize=16)

ax[5].plot(ops['t'][0:kk],y[11,0:kk], drawstyle='steps')
ax[5].grid(True)
ax[5].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[5].set_ylabel(r'$y_{12}$ (kW)', fontsize=16)

ax[6].plot(ops['t'][0:kk],y[12,0:kk], drawstyle='steps')
ax[6].grid(True)
ax[6].set_xlim(left=ops['t'][0], right=ops['t'][-1])
ax[6].set_ylabel(r'$y_{13}$ (kW)', fontsize=16)

ax[6].set_xlabel(r'Time (hours)',fontsize=16)
