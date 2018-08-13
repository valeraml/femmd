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


def run2():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))

    dt = 0.002
    ucells = 2
    sigma = 1
    N = ucells*ucells
    dens = 0.1
    nu = 0.5
    L = math.sqrt(N / dens) 
    print(L)
    T = 2.0
    v0 = math.sqrt(2.0 * T )
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    rcut=np.array([rc])
    shift=np.array([1.0])
    E=200
    k=50

    pdb.set_trace()
    box = np.array([L,L])
    md.system_init(box)
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.3)
    L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc')
    box = md.system_get_box()
    print("box",box)
    md.system_set_boundary_conditions(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_set_youngs_modulus(E)
    #md.system_set_swelling_energy_constants(100,50,0)
    #md.system_set_friction(.1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    bonds = md.system_get_bonds()
    #plt.scatter(pos[:,0],pos[:,1])
    #plt.show()
    #plt.triplot(pos[:,0], pos[:,1], tet, 'go-', lw=1.0)
    bonds = md.system_get_bonds()
    refVol, currVol = md.system_get_elements_volumes()
    I1,I2 = md.system_get_elements_invariants()
    print("Volume check: ",refVol.shape[0], currVol.shape[0], np.sum(refVol), np.sum(currVol))
    pdb.set_trace()

    sigmas = md.system_get_particle_sizes()
    sigmas *= sigma
    R = R - 0.5 + sigma/2

    #xyzfile = fm.Saver(0,"test_test",pos)
    #vtkfile = fm.Saver(1,"test_test",upos,vel=vel,cells=tet)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[box[0],box[1]])
    #vtffile.append(0);
    plot = Plotter.Plotter(pos,vel,tet,box[0],box[1],'tri',I1,sigma) # 2 is the radius
    plot.update(0)
    #plt.show()

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


run2()

