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


def run0():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))

    dt = 0.005
    ucells = 10
    sigma = 1
    N = ucells*ucells
    dens = 0.1
    nu = 0.3
    L = math.sqrt(N / dens) 
    print(L)
    T = 2.0
    v0 = math.sqrt(2.0 * T )
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    rc = 2.5
    rcut=np.array([rc])
    shift=np.array([0.0])
    E=200
    k=50

    pdb.set_trace()

    md.system_init(np.array([L,L]))
    fm.init_particles_simple_cubic(md,ucells,L)
    lbox = md.system_get_box()
    print("box",lbox)
    md.system_set_boundary_conditions(2)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_init_neighbor_list()
    pos = md.system_get_positions1()
    vel = md.system_get_velocities()

    #upos = md.system_get_unfolded_positions()
    #tet = md.system_get_tetras()
    #bonds = md.system_get_bonds()

    #sigmas = md.system_get_particle_sizes()
    #sigmas *= sigma
    #R = R - 0.5 + sigma/2

    #xyzfile = fm.Saver(0,"test_test",pos)
    #vtkfile = fm.Saver(1,"test_test",pos,vel=vel,cells=tet)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[lbox[0],lbox[1]])
    #vtffile.append(0);
    
    plot = Plotter.Plotter(pos,vel,None,lbox[0],lbox[1],'pts',sigma) # 2 is the radius
    plot.update(0)
    #r = 0.0
    #md.system_set_walls_moving_rate(r,r)
    #pdb.set_trace()
    for i in range(1000):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        #xyzfile.append(avgs[0])
        #vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #md.system_unfold_positions()
        plot.update(avgs[0])

    print("finished")
    #plt.show()

run0()