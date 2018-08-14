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

# just saving this, T1 transition with 8 particles succesful
def runT1_8():

    print(os.getcwd())
    os.chdir('c:\\tmp\\test1')
    print(os.getcwd())
    print(dir(md))
    #md.system_c_main()

    dt = 0.002
    ucells = 3
    sigma = 2
    N = ucells*ucells
    dens = 0.1
    #nu = 0.5
    nu = 0.15
    dens = nu/np.pi
    L = math.sqrt(N / dens) 
    print(L)
    T = 1.0
    #v0 = math.sqrt(2.0 * T * (1.0-1.0/N))
    v0 = math.sqrt(2.0 * T )
    eps=np.array([1.0])
    rc = sigma*math.pow(2,1.0/6.0)
    rcut=np.array([rc])
    shift=np.array([1.0])
    E=200
    k=200

    box = np.array([L,L])
    md.system_init(box)
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    nodes, cells = mdmesh.uniform_mesh_on_unit_circle(.25)
    L,R = fm.add_clusters(md,nodes,cells,1.0,ucells,nu,'sc1')
    #L,R = fm.add_binary_clusters(md,nodes,cells,1.0,1.0,.15,ucells,nu,'fcc')
    box = md.system_get_box()
    md.system_set_boundary_conditions(0)
    md.system_set_moving_walls(0)
    md.system_set_velocities(v0)
    md.system_set_zero_center_of_mass_vel()
    md.system_set_dt(dt)
    md.system_set_avg_steps(100)
    md.system_set_potential(eps,rcut,shift)
    md.system_set_youngs_modulus(E)
    #md.system_set_swelling_energy_constants(100,250,0.0)
    md.system_set_friction(1)
    #md.system_set_bond_spring_constant(k)
    md.system_init_neighbor_list()
    pos = md.system_get_positions()
    vel = md.system_get_velocities()
    upos = md.system_get_unfolded_positions()
    tet = md.system_get_tetras()
    #plt.scatter(pos[:,0],pos[:,1])
    #plt.show()
    #plt.triplot(pos[:,0], pos[:,1], tet, 'go-', lw=1.0)
    bonds = md.system_get_bonds()
    refVol, currVol = md.system_get_elements_volumes()
    I1,I2 = md.system_get_elements_invariants()
    print("Volume check: ",refVol.shape[0], currVol.shape[0], np.sum(refVol), np.sum(currVol))

    sigmas = md.system_get_particle_sizes()
    sigmas *= sigma
    R = R - 0.5 + sigma/2
    n1 = 0.79
    L1 = np.sqrt(ucells**2*np.pi*R**2/n1)
    nf = 0.95
    Lf = np.sqrt(ucells**2*np.pi*R**2/nf)
    print(L,Lf,R)

    #xyzfile = fm.Saver(0,"test_test",pos)
    vtkfile = fm.Saver(1,"test_test",pos,vel=vel,cells=tet, I2=I2)
    #vtffile = fm.Saver(2,"test_test",upos,bonds=bonds,box=[Lx,Lx])
    vtkfile.append(0);
    walls = md.system_get_walls_pos()
    print("walls",walls)
    walls = [walls[0][0],walls[1][0],walls[0][1],walls[1][1]]
    print("walls",walls)
    plot = Plotter.Plotter(pos,vel,tet,box[0],box[1],'tripts1',I1,sigma) # 2 is the radius
    plot.update(0,walls,"")
    #plt.show()
    r = Lf/800
    r00 = np.array([1.01*r,r])
    r10 = np.array([-1.01*r,-r])
    r0 = np.array([r, r])
    r1 = np.array([-r,-r])
    r02 = np.array([-r,r])
    r12 = np.array([r,-r])
    r_stop = np.array([0,0])
    data = []
    totRefVol = np.sum(refVol)
    totCurrVol = np.sum(currVol)
    md.system_set_moving_walls(0)
    md.system_set_walls_moving_rate(r0,r1)  
    pdb.set_trace()  
    shear = 0

    Ws = 0
    for Ws in np.arange(0,250,10):
        md.system_run(100);
        avgs = md.system_get_avgs()
        print("\nstep %d" % avgs[0])
        print(avgs)
        title = "steps = %d" % (avgs[0])
        #if Ws % 10 == 0:
            #Ws += 10 
        md.system_set_swelling_energy_constants(100,Ws,0.0)
        plot.update(avgs[0],walls,title)
        vtkfile.append(avgs[0])


    md.system_set_moving_walls(1)
    md.system_set_walls_moving_rate(r0,r1)
    for i in range(10000):
        #pdb.set_trace()
        md.system_run(100);
        avgs = md.system_get_avgs()
        walls = md.system_get_walls_pos()
        print("walls",walls)
        walls = [walls[0][0],walls[1][0],walls[0][1],walls[1][1]]
        print("\nstep %d" % avgs[0])
        print("walls",walls)
        Lx = walls[1] - walls[0]
        Ly = walls[3] - walls[2]
        print("Box Vol = ",Lx*Ly)
        
        if i == 20:
            md.system_set_walls_moving_rate(r00,r10)
        if Ly/Lx > 1.1 and shear == 0:
            md.system_set_walls_moving_rate(r0,r1)
        if totCurrVol/totRefVol < 1.5 and i > 20:
            shear = 1;
        if shear == 1:
            c = Lx/Ly
            r02 = np.array([-c*r,r])
            r12 = np.array([c*r,-r])
            md.system_set_walls_moving_rate(r02,r12) 
            #md.system_set_walls_shear(1)
        if Lx > L and i > 1000:
            md.system_set_walls_moving_rate(r_stop,r_stop)
        
        md.system_unfold_positions()
        #xyzfile.append(avgs[0])
        vtkfile.append(avgs[0])
        #vtffile.append(avgs[0])
        #if i % 1 == 0:
        print(avgs)
        wf = md.system_get_walls_forces()
        #return Py_BuildValue("((dd)(dd))", f[0][0], f[0][1], f[1][0], f[1][1]);
        wf = [wf[0][0],wf[0][1],wf[1][0],wf[1][1]]
        print("force on walls :", wf)
        totCurrVol = np.sum(currVol)
        data.append([dt*avgs[0],wf[0],wf[1],wf[2],wf[3],totCurrVol/totRefVol,avgs[5]])
        print("Volumes: ",totRefVol, totCurrVol, totCurrVol/totRefVol*100)
        print("I1 max", I1.max());
    
        title = "steps = %d Lini = %f Lx = %f Ly = %f nu = %f V/Vref= %f" % (avgs[0],L, Lx, Ly, ucells**2*np.pi*R**2/(Lx*Ly), totCurrVol/totRefVol)
        if i % 10 == 0: 
            plot.update(avgs[0],walls,title)
            wf = np.array(data)
            np.savetxt("data.dat",wf)
            vtkfile.append(avgs[0])

    wf = np.array(data)
    np.savetxt("data.dat",wf)
    print("finished")
    plt.show()

runT1_8()