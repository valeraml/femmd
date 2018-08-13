import os
import femmd as md
import femmd_module as fm
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import mdmesh
import Plotter
import pdb

import mpl_toolkits.mplot3d as a3


def run3D():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    #md.system_c_main()

    dt = 0.005
    ucells = 1
    sigma = 1
    N = ucells**3
    dens = 0.1
    nu = 0.1
    L = math.pow(N / dens,1/3.0) 
    print(L)
    T = 2.0
    v0 = math.sqrt(2.0 * T )
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    #rc = 2.5
    rcut=np.array([rc])
    shift=np.array([1.0])
    E=200
    k=50

    pdb.set_trace()

    md.system_init(np.array([L,L,L]))
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_ball(.2)
    L,R = fm.add_3D_clusters(md,nodes,cells,1.0,ucells,nu,'sc1')
    #fm.init_particles_simple_cubic3D(md,ucells,L)
    lbox = md.system_get_box()
    md.system_set_boundary_conditions(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_set_youngs_modulus(E)
    md.system_set_swelling_energy_constants(100,10,0)
    md.system_set_friction(1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    print(pos)
    vel = md.system_get_velocities()
    #vel += 1
    #upos = md.system_get_unfolded_positions()
    #tet = md.system_get_tetras()
    #bonds = md.system_get_bonds()

    #sigmas = md.system_get_particle_sizes()
    #sigmas *= sigma
    #R = R - 0.5 + sigma/2

    xyzfile = fm.Saver(0,"test_test",pos,dims=3)
    vtkfile = fm.Saver(1,"test_test",pos,dims=3,vel=vel,cells=tet)
    vtkfile.append(0)
    '''
    lx = lbox[0]
    fig = plt.figure()
    ax = a3.Axes3D(fig)
    ax.set_aspect('equal')
    ax.set_xlim(0,lx)
    ax.set_ylim(0,lx)
    ax.set_zlim(0,lx)
    fig.set_size_inches(15, 15)
    ax.scatter(pos[:,0],pos[:,1],pos[:,2])
    ax.add_collection3d(tri)
    plt.show()
    '''

    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[lbox[0],lbox[1]])
    #vtffile.append(0);
    #plot = Plotter.Plotter(pos,vel,tet,lbox[0],lbox[1],'pts',sigma) # 2 is the radius
    #plot.update(0)
    #r = 0.0
    #md.system_set_walls_moving_rate(r,r)
    #pdb.set_trace()
    for i in range(1000):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        xyzfile.append(avgs[0])
        vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #md.system_unfold_positions()
        #plot.update(avgs[0])

    print("finished")
    plt.show()

run3D()