import os
import sys
sys.path.append("../src/pysrc")
import femmd as md
import femmd_module as fm
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import mdmesh
import Plotter
import pdb

# In this script, we set up 4 elastomers, or bodies in a box a let it evolve
# over time.

def run0():

    # Switch to data directory
    print(os.getcwd())
    data_dir = 'data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    os.chdir(data_dir)
    print(os.getcwd())
    print(dir(md))

    # set delta t for time evolution
    dt = 0.002
    # number of cell
    ucells = 2
    # size of particles
    sigma = 1
    # total number of elastormes 4, in this case
    N = ucells*ucells
    # set the density and packing fraction, will be modified later
    # density
    dens = 0.1
    # packing fraction
    nu = 0.5
    # size of the box
    L = math.sqrt(N / dens) 
    print(L)
    # temperature and initial average velocity
    T = 2.0
    v0 = math.sqrt(2.0 * T )
    # potetial parameters
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    rcut=np.array([rc])
    shift=np.array([1.0])

    # young's modulus
    E=200

    # set the box size and initialize the system
    box = np.array([L,L])
    md.system_init(box)

    # The elastomers and structure as a mesh with nodes and cells (triangle)
    # you can read the nodes and triangles from a file or you can 
    # use the function mdmesh to produce the shape, in this example we 
    # use mdmesh, it returs the nodes and cells as an array of vector positions and
    # another array of index of vertices for triangles, the .3 controls how 
    # fine we want the mesh, lower numbers will give you more triangles
    # in this case the mesh is a disk of diameter 1
    
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.3)

    # Now that we have the basic mesh for one elastomer we will replicate it 4 times
    # in the box, we will also rescale the size of of the box and meshes so that ieach 
    # vertices contains a particle of diameter 1. Refore to function in  add cluster
    # to see how it is done
    L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc')

    # get the new size of the box
    box = md.system_get_box()
    print("box",box)

    # set boundary conditions, 0 for hard walls in all directions
    md.system_set_boundary_conditions(0)
    # set initial velocities
    md.system_set_velocities(v0)
    # substract center of mass
    md.system_set_zero_center_of_mass_vel()
    # set delta t for the sytem
    md.system_set_dt(dt)
    # set average steps
    md.system_set_avg_steps(100)
    # set the potential, interaction between particles
    md.system_set_potential(eps,rcut,shift)
    # set the youngs modulus for elactic interaction
    md.system_set_youngs_modulus(E)

    # if we wnat to initially swell the particles , uncomment this line
    #md.system_set_swelling_energy_constants(100,50,0)
    # if we want to get rid of thermal fluctiation, uncomment this line
    #md.system_set_friction(.1)
    #if we want to simulate a system with spring interactions instead, we use
    # the line below, but we would need to set k, and set E = 0
    #md.system_set_bond_spring_constant(k)

    # initialize neighbor list
    md.system_init_neighbor_list()
    # get pointers to positions, velocities and unfolded positions
    # also get a pointer ot the array of triangles
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    
    # get pointers to arrays fo initial area (refVol), and 
    # current area, (currVol) of triangles
    # also, array with the invariants of the triangles

    refVol, currVol = md.system_get_elements_volumes()
    I1,I2 = md.system_get_elements_invariants()
    print("Volume check: ",refVol.shape[0], currVol.shape[0], np.sum(refVol), np.sum(currVol))
    #pdb.set_trace()

    # we rescale the particles at the vertices 
    sigmas = md.system_get_particle_sizes()
    sigmas *= sigma
    R = R - 0.5 + sigma/2

    #if we want to save the configuration, we can use the formats below
    #xyzfile = fm.Saver(0,"test_test",pos)
    #vtkfile = fm.Saver(1,"test_test",upos,vel=vel,cells=tet)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[box[0],box[1]])
    #vtffile.append(0)
    
    # lets create a plot for animation
    plot = Plotter.Plotter(pos,vel,tet,box[0],box[1],'tri',I1,sigma) # 2 is the radius
    plot.update(0)
    
    # run the simulation
    for i in range(1000):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        #xyzfile.append(avgs[0])
        #vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #md.system_unfold_positions()
        if i % 1 == 0: plot.update(avgs[0])

    print("finished")
    plt.show()


run0()

