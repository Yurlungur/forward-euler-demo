#!/usr/bin/env python2

# planet_forward_euler.py
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2016-05-09 15:54:04 (jmiller)>

# This is a simple program which demonstrates how to use the Forward
# Euler method to solve for the motion of a planet near a star.

# Because it's easier to visualize, we assume two dimensions, not
# three.

# imports
# ----------------------------------------------------------------------
import numpy as np # for handling arrays
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation
# ----------------------------------------------------------------------

# constants
# ----------------------------------------------------------------------
G = 6.6741E-11 # Newton's constant
Msun = 2E30 # mass of star
Rearth = 150E6 # Distance from planet to star
t0 = 0 # initial time
tf = 3*(365) # final time
nt = 30*(tf+1) # number of Euler steps
animation_name = 'forward-euler-ypert-{}.mp4'
# ----------------------------------------------------------------------

def get_r(x,y):
    """Given x and y coordinates of the planet,
    returns distance to star."""
    return np.sqrt(x**2 + y**2)

def get_acceleration(x,y):
    """Given x and y coordinates of the planet, returns the magnitude of
    the acceleration due to gravity from the star.
    """
    r = get_r(x,y)
    mag_a = -1.0*(G*Msun)/(r**3)
    ax = mag_a*x
    ay = mag_a*y
    return ax,ay

def rhs(p):
    """Given a time t and a state vector p, returns
    the right-hand-side of the ODE system.
    p = [x,y,vx,vy]
    where (x,y) are the coordinates of the planet
    and (vx,vy) are the velocity.
    The ODE is
    dp/dt = [vx,vy,ax,ay]
    where (ax,ay) are the x- and y-components of the
    acceleration due to gravity.
    """
    x,y,vx,vy = p # unwrap
    ax,ay = get_acceleration(x,y)
    return np.array([vx,vy,ax,ay])

def forward_euler_step(p,dt):
    dp = dt*rhs(p)
    return p + dp

def forward_euler_integration(x0,y0,vx0,vy0,t0,tf,nt):
    """Given initial conditions
    p = [x0,y0,vx0,vy0]
    at time t0,
    solves the system up to time tf with nt euler steps.
    Returns 2 arrays:
    t, [x,y,vx,vy]
    where the second array, called the state,
    is a (nt X 4) matrix."""
    t = np.linspace(t0,tf,nt)
    dt = t[1] - t[0]
    state = np.empty((len(t),4),dtype=float)
    state[0,...] = np.array([x0,y0,vx0,vy0])
    print "Starting the integration."
    print "dt = {}".format(dt)
    for i in range(1,len(t)):
        state[i,...] = forward_euler_step(state[i-1,...],dt)
    return t,state

def get_initial_data(vy_pert):
    """Produces initial position and velocity data
    based on the orbit of the Earth around the sun.
    When theta_pert = 0
    we have initial data for a circular orbit.
    When theta_pert != 0 we get something else.
    theta_pert is just an angle. Radius and velocity are assumed
    to be earth-like."""
    r = Rearth
    theta = -np.pi/4
    x,y = r*np.cos(theta),r*np.sin(theta)
    mag_v = np.sqrt(G*Msun/r) # Kepplerian velocity
    vx,vy = -mag_v*np.sin(theta),mag_v*np.cos(theta) # for circular orbit
    # rotate based on vy_pert
    vy_new = vy + vy_pert*mag_v
    vx_new = vx 
    return np.array([x,y,vx_new,vy_new])

def animate_trajectory(vy_pert,t0,tf,nt,ypertname):
    print "-------------------------------------"
    FRAME_FACTOR=180
    p0 = get_initial_data(vy_pert)
    x0,y0,vx0,vy0 = p0
    print "Initial conditions:"
    print "\tx = {}\n\ty = {}\n\tvx = {}\n\t vy = {}".format(x0,y0,
                                                             vx0,vy0)
    if ypertname == 'ellipse':
        tf *= 2.75

    t,state = forward_euler_integration(x0,y0,vx0,vy0,t0,tf,nt)

    print "Simulation finished. Preparing visualization."

    xmin = np.min(state[...,0])
    ymin = np.min(state[...,1])
    xmax = np.max(state[...,0])
    ymax = np.max(state[...,1])

    xmin *= 1.1
    ymin *= 1.1
    xmax *= 1.1
    ymax *= 1.1
    print "xlim = [{}, {}]".format(xmin,xmax)
    print "ylim = [{}, {}]".format(ymin,ymax)

    fig = plt.figure()
    ax = plt.axes(xlim=(xmin,xmax),
                  ylim=(ymin,ymax))
    line, = ax.plot([],[],'b-',lw=2)
    sun, = ax.plot(0,0,'yo',ms=21)
    planet, = ax.plot([],[],'go',ms=7)

    def init():
        line.set_data([],[])
        planet.set_data([],[])
        return line,planet

    def animate(mi):
        i = FRAME_FACTOR*mi
        line.set_data(state[:i+1,0],state[:i+1,1])
        planet.set_data(state[i,0],state[i,1])
        return line,planet
    
    ani = animation.FuncAnimation(fig,animate,
                                  frames=nt/FRAME_FACTOR,
                                  interval=1./30.,
                                  blit=True,init_func=init)
    
    print "saving animation"
    print "num frames = {}".format(nt/FRAME_FACTOR)
    ani.save(animation_name.format(ypertname),fps=30,
             extra_args=['-vcodec','libx264'])
    print "done!"
    print "-------------------------------------"

if __name__ == "__main__":
    animate_trajectory(0,t0,tf,nt,'circle')
    animate_trajectory(-0.9,t0,tf,nt,'spiral')
    animate_trajectory(0.3,t0,tf,nt,'ellipse')
    animate_trajectory(0.7,t0,tf,nt,'hyperbola')
