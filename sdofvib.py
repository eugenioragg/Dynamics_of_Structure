# -*- coding: utf-8 -*-
"""
Sdofvib
-------

Driver for sdof free and forced response using:      
    timestep.m

    Department of Mechanical Engineering
    Technical University of Denmark

Created: June 2020

@author: Oscar Bondo Ellekvist, s163774

"""

# Import packages
import numpy as np
import matplotlib.pyplot as plt

# Import functions
from timestep import timestep


### DATA SECTION -----------------------------------------------

# Structural parameters
m = 1000                # mass [kg]
omega0 = 8              # natural (angular) frequency [rad/s]
zeta = 0.02             # damping ratio [-]

k = np.square(omega0)*m # stiffness from frequency [N/m]
c = 2*zeta*omega0*m     # viscous coefficient from damping ratio [N/(m/s)]

T0 = 2*np.pi/omega0     # vibration period

# Time record parameters
t0 = 0                  # initial time
dt = 0.01*T0            # time increment (here 100 points pr. period T0)
n = 1000                # number of increments (here 10 periods)

# forcing parameters
omega = omega0                          # forcing (angular) frequency [rad/s]
f0 = 100                                # forcing amplitude [N]
t = np.arange(0,n+1,1,dtype=float)*dt   # time record
#f = f0*np.sin(omega*t)                  # harmonic load record
f = np.zeros(n+1,dtype=float)           # zero load record (free vibrations)

# Initial conditions
x0 = 0.1
v0 = 1.0


### RESPONSE CALCULATIONS --------------------------------------

# Call timestep function - leave out f for free vibrations
(x,v,a,f,t) = timestep(m,c,k,x0,v0,t0,dt,n,f)


### POST-PROCESSING

# Save data to .csv file
np.savetxt('sdofvib_response.csv',np.stack((t,x,v,a),axis=1),delimiter=',')


# Plot response history using matplotlib
plt.close('all')

fig1 = plt.figure(1,figsize=[12,4])      # plot of disp and vel

fig1.add_subplot(121)           # subplot 1 = disp
plt.plot(t/T0,x,color='red')
plt.title('Displacement')
plt.xlabel(r'$t/T_0$')
plt.ylabel(r'$x(t)$')
plt.grid()
plt.xlim([0, 10])
plt.ylim([-0.2, 0.2])

fig1.add_subplot(122)           # subplot 2 = vel
plt.plot(t/T0,v/omega0,color='blue')
plt.title('Velocity')
plt.xlabel(r'$t/T_0$')
plt.ylabel(r'$\dot{x}(t)$')
plt.grid()
plt.xlim([0, 10])
#plt.ylim([-0.2, 0.2])
plt.tight_layout(pad=3)   # spacing between subplots
plt.show()


fig2 = plt.figure(2,figsize=[6,4],)       # plot of disp and acc
plt.plot(t/T0,x,color='blue')
plt.plot(t/T0,a/np.square(omega0),color='red')
plt.title('Displacement and Acceleration')
plt.xlabel(r'$t/T_0$')
plt.ylabel(r'$x(t)$, $\ddot{x}(t)/\omega_0^2$')
plt.grid()
plt.xlim([0, 10])
plt.legend([r'$x(t)$', r'$a$'])
#plt.ylim([-0.2, 0.2])
plt.tight_layout()


