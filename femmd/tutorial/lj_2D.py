import os
import sys
#sys.path.append("../src/pysrc")
sys.path.append("../femmdlib")
import femmd as md
import femmd_module as fm
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import mdmesh
import Plotter
import pdb

#This script runs a 2D Lennard Jones molecular dynamics simulation
#the purpose is to show how the code works.

def run0():


    #Create a data directory for all data produced by the program
    print(os.getcwd())
    data_dir = 'data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    os.chdir(data_dir)
    print(os.getcwd())
    print(dir(md))

    # define some variables
    # delta time, time increment
    dt = 0.005
    # number of cells, will be used to assign the initial position  of
    # particles and for size of system
    ucells = 10
    # diameter of particles
    sigma = 1
    # Number of particles, we have one particle per cell
    N = ucells*ucells
    # Density, N/V, V volume
    dens = 0.1
    # packing fraction N*sigme*3/V
    nu = 0.3
    # size of of square
    L = math.sqrt(N / dens) 
    print(L)
    # Initial temperature
    T = 2.0
    # Initial average velocity of particles
    v0 = math.sqrt(2.0 * T )
    # lennard jones potential epsilon, strenght of potential
    # it has to be declared using a list in case we have more type of particles
    eps=np.array([1.0])
    # Cut for the potetial, 
    rc = sigma*math.pow(2,1.0/6.0)
    rc = 2.5
    rcut=np.array([rc])
    # potential zero reference
    shift=np.array([0.0])

    #Now we create the simulation system

    # initialize the system with a box size of [L,L]
    md.system_init(np.array([L,L]))
    # set the particles on a square lattice
    fm.init_particles_simple_cubic(md,ucells,L)
    lbox = md.system_get_box()
    print("box",lbox)
    #set boundary condiations 
    # we use the following construct in C++ enum PBCTYPE {NOPBC, XPBC, XYPBC, XYZPBC};
    # 0, is no periadic boundary conditions, all hard walls
    # 1 is boundary conditions in the x direction
    # 2 is boundary conditions in both x and y direction
    md.system_set_boundary_conditions(2)
    #set the initial veloctities random with average v0
    md.system_set_velocities(v0)
    #substract the center of mass velocity
    md.system_set_zero_center_of_mass_vel()
    # set delta_t for the system
    md.system_set_dt(dt)
    # number of steps to average
    md.system_set_avg_steps(100)
    # set the potential using the varialbles eps, rcut and shift defined before
    md.system_set_potential(eps,rcut,shift)
    # initialize verlet list
    md.system_init_neighbor_list()
    # Get a pointer to array with position and velocities for later work
    pos = md.system_get_positions1()
    vel = md.system_get_velocities()

    # we can save the position using any of these formats, uncomment the wa
    #xyzfile = fm.Saver(0,"test_test",pos)
    #vtkfile = fm.Saver(1,"test_test",pos,vel=vel,cells=tet)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[lbox[0],lbox[1]])
    #vtffile.append(0);
    
    # If you want to plot the particles along with the simulation,
    #uncomment these lines
    #plot = Plotter.Plotter(pos,vel,None,lbox[0],lbox[1],'pts',sigma) # 2 is the radius
    #plot.update(0)

    for i in range(1000):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        #xyzfile.append(avgs[0])
        #vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #md.system_unfold_positions()
        #plot.update(avgs[0])

    print("finished")


run0()
