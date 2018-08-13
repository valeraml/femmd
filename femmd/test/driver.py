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

import timeit

#Simple leenard jones for testing
def run0():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    #md.system_c_main()

    dt = 0.005
    ucells = 20
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

    md.system_init(L,L)
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    #nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.2)
    #L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc11')
    fm.init_particles_simple_cubic(md,ucells,L)
    Lx,Ly = md.system_get_box()
    md.system_set_boundary_conditions(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    #md.system_set_youngs_modulus(E)
    #md.system_set_swelling_energy_constants(100,100,0)
    #md.system_set_friction(.1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    #vel += 1
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    bonds = md.system_get_bonds()

    sigmas = md.system_get_particle_sizes()
    #sigmas *= sigma
    #R = R - 0.5 + sigma/2

    xyzfile = fm.Saver(0,"test_test",pos)
    vtkfile = fm.Saver(1,"test_test",upos,vel=vel,cells=tet)
    vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[Lx,Lx])
    vtffile.append(0);
    plot = Plotter.Plotter(pos,vel,tet,Lx,Lx,'pts',sigma) # 2 is the radius
    plot.update(0)
    r = 0.0
    md.system_set_walls_moving_rate(r,r)
    #pdb.set_trace()
    for i in range(100):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        xyzfile.append(avgs[0])
        vtkfile.append(avgs[0])
        vtffile.append(avgs[0])
        md.system_unfold_positions()
        plot.update(avgs[0])

    print("finished")
    plt.show()

def run2():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    #md.system_c_main()

    dt = 0.005
    ucells = 2
    sigma = 1
    N = ucells*ucells
    dens = 0.1
    nu = 0.2
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

    md.system_init(L,L)
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.2)
    L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc')
    #fm.init_particles_simple_cubic(md,ucells,L)
    Lx,Ly = md.system_get_box()
    md.system_set_boundary_conditions(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_set_youngs_modulus(E)
    #md.system_set_swelling_energy_constants(100,50,0)
    md.system_set_friction(.1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    #vel += 1
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    bonds = md.system_get_bonds()

    sigmas = md.system_get_particle_sizes()
    sigmas *= sigma
    R = R - 0.5 + sigma/2

    xyzfile = fm.Saver(0,"test_test",pos)
    vtkfile = fm.Saver(1,"test_test",upos,vel=vel,cells=tet)
    vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[Lx,Lx])
    vtffile.append(0);
    plot = Plotter.Plotter(pos,vel,tet,Lx,Lx,'tripts1',sigma) # 2 is the radius
    plot.update(0)
    r = 0.01
    md.system_set_walls_moving_rate(r,r)
    #pdb.set_trace()
    for i in range(1000):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print(avgs)
        xyzfile.append(avgs[0])
        vtkfile.append(avgs[0])
        vtffile.append(avgs[0])
        md.system_unfold_positions()
        if i % 10 == 0: plot.update(avgs[0])

    print("finished")
    plt.show()



def run1():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    #md.system_c_main()

    dt = 0.005
    ucells = 2
    sigma = 2
    N = ucells*ucells
    dens = 0.1
    nu = 0.2
    dens = nu/np.pi
    L = math.sqrt(N / dens) 
    print(L)
    T = 1.1
    #v0 = math.sqrt(2.0 * T * (1.0-1.0/N))
    v0 = math.sqrt(2.0 * T )
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    rcut=np.array([rc])
    shift=np.array([1.0])
    E=200
    k=200

    md.system_init(L,L)
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.1)
    L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc1')
    #init_particles_simple_cubic(ucells,L)
    Lx,Ly = md.system_get_box()
    md.system_set_boundary_conditions(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_set_youngs_modulus(E)
    md.system_set_swelling_energy_constants(100,50,0)
    md.system_set_friction(1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    #vel += 1
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    bonds = md.system_get_bonds()
    refVol, currVol = md.system_get_elements_volumes()

    sigmas = md.system_get_particle_sizes()
    sigmas *= sigma
    R = R - 0.5 + sigma/2
    n1 = 0.79
    L1 = np.sqrt(ucells**2*np.pi*R**2/n1)
    nf = 0.95
    Lf = np.sqrt(ucells**2*np.pi*R**2/nf)
    print(L,Lf,R)

    #xyzfile = fm.Saver(0,"test_test",pos)
    vtkfile = fm.Saver(1,"test_test",upos,vel=vel,cells=tet)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[Lx,Lx])
    #vtffile.append(0);
    plot = Plotter.Plotter(pos,vel,tet,Lx,Lx,'tripts1',sigma) # 2 is the radius
    plot.update(0)
    r = Lf/800
    md.system_set_walls_moving_rate(0,0)
    data = []
    totRefVol = np.sum(refVol)
    md.system_set_walls_moving_rate(r,r)    
    for i in range(10):
        md.system_run(100);
        avgs = md.system_get_avgs()
        walls = md.system_get_walls_pos()
        #print(walls)
        Lx = walls[1] - walls[0]
        Ly = walls[3] - walls[2]
        #if i == 20:
        #    md.system_set_walls_moving_rate(r,r)
        if Ly < L1 and Ly > Lf:
            md.system_set_walls_moving_rate(-r,r)
        if Lx > L:
            md.system_set_walls_moving_rate(0.0,0.0)
        md.system_unfold_positions()
        #xyzfile.append(avgs[0])
        #vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #if i % 1 == 0:
        print(avgs)
        wf = md.system_get_walls_forces()
        print("force on walls :", wf)
        data.append([dt*avgs[0],wf[0],wf[1],wf[2],wf[3]])
        totCurrVol = np.sum(currVol)
        print("Volumes: ",totRefVol, totCurrVol, totCurrVol/totRefVol*100)
    
        title = "steps = %d Lini = %f Lx = %f Ly = %f nu = %f" % (avgs[0],L, Lx, Ly, ucells**2*np.pi*R**2/(Lx*Ly))
        if i % 1 == 0: plot.update(avgs[0],walls,title)

    wf = np.array(data)
    np.savetxt("data.dat",wf)
    print("finished")
    #plt.show()


#run2()
#md.system_c_main()
a = timeit.timeit(run0,number=1)
print(a)


