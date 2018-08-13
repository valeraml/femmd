import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import time

def init_particles_simple_cubic(md,unitCells,L):
    
    dx = L/unitCells;
    dy = L/unitCells;
    n = 0
    for ny in np.arange(unitCells):
        for nx in np.arange(unitCells):
            x = dx*(nx + 0.5)
            y = dy*(ny + 0.5)
            vx = (2*np.random.rand()-1)
            vy = (2*np.random.rand()-1)
            #print(x,y,vx,vy)
            md.system_add_particle(x,y,vx,vy,1.0,1.0,0)


#Given and number of cluster nc = ucells**2, and a density dens
#sets nc clusters in a single cubic.
#calculate a scale factor assuming nodes have particles of diameter 1
#area of cluster has to be 1
# 
def add_clusters(md,p,t,R0,ucells,nu_of_system,type='sc'):
    
    #nodes = np.loadtxt("diskR1Vertices1.txt");
    #cells = np.loadtxt("diskR1Triangles1.txt") - 1;
    #nodes, cells = fm.uniform_mesh_on_unit_circle(.25)
    nodes = p
    cells = t
    cells = cells.astype(int)
    print(nodes.shape, cells.shape)
    print(nodes.dtype, cells.dtype)

    # Need to calculate a scale factor for the box or new box size for given density
    # we want disks with particles of diameter 1 and nu = 0.91 ad hoc
    # Initial Radius of disk is R0 = 1
    # Initial Area of disk is A0 = PI
    # nu = n_vertices*pi*r0^2/A0
    # r0 radius of inner particles when disk has Radius 1
    # r0 -> r = 0.5
    # r0 = sqrt(nu*R0/n_vertices)
    #R0 = 1
    n_vertices = nodes.shape[0]
    nu = 0.91  # ad hoc packing fraction, max packing fraction for disks is 0.9069
    #nu = Nv*pi*r0^2/(pi*R0^2)
    r0 = R0 * np.sqrt(nu/ n_vertices); # all calculation without outer radius of particle included
    r0 = R0 * np.sqrt(nu/n_vertices)/(1.0-np.sqrt(nu/n_vertices))
    r0 = R0 / (np.sqrt(n_vertices/nu) - 1)
    r = 0.5
    # new disk radius from keeping same packing fraction inside disk r0/R0 = r/R use r=.5
    R = r*R0/r0 #new radious of disk with n_vertices of diameter 1 and nu = 0.91
    # L0 old box size, L new box size, keeping same packing fraction
    # R0/L0 = R/L
    N = ucells**2
    if type == 'fcc': N = 2*N
    #L0 = np.sqrt(N*np.pi*R0**2/nu_of_system)
    L0 = np.sqrt(N*np.pi*(R0+r0)**2/nu_of_system)
    L = R*L0/R0
    scale = L/L0
    md.system_set_scale_factor(scale);
    print("r0 %f scale %f r %f R %f A0 %f A %f expected A %f\n" % ( r0, scale, scale*r0, R, np.pi*r0*r0, np.pi*r0*r0*scale*scale, np.pi));
    print("L0 = %f, L = %f\n" % (L0, L))
	#scale simulation box
    md.system_set_box(L,L);
	#scale particles pos
    nodes *= scale
    dx = L/ucells
    dy = dx
    mass = 1.0*n_vertices
    exclude = 1
    if ucells == 1:
        x,y = L/2,L/2
        print("adding cluster at %f %f %f %f\n" % (x,y,dx,dy))
        md.system_add_cluster(nodes,cells,x,y,0,mass,exclude)
    elif type == 'sc':
        for ny in np.arange(ucells):
            for nx in np.arange(ucells):
                x = dx*(nx + 0.5)
                y = dy*(ny + 0.5)
                #PyArg_ParseTuple(pArgs, "OOddidi", &vert, &cells, &x, &y, &typ, &mass, &exc)
                print("adding cluster at %f %f %f %f\n" % (x,y,dx,dy))
                md.system_add_cluster(nodes,cells,x,y,0,mass,exclude)
    elif type == 'fcc':
        a1 = np.array([0.0,0.0])
        a2 = np.array([0.5,0.5])
        dr = np.array([dx,dy])
        for ny in np.arange(ucells):
            for nx in np.arange(ucells):
                nr = np.array([nx,ny])
                r1 = dr*(a1 + nr + 0.25)
                r2 = dr*(a2 + nr + 0.25)
                print(R,r1,r2)
                #PyArg_ParseTuple(pArgs, "OOddidi", &vert, &cells, &x, &y, &typ, &mass, &exc)
                print("adding  2 clusters at cell %f %f %f %f\n" % (r1[0],r1[1],dx,dy))
                md.system_add_cluster(nodes,cells,r1[0],r1[1],0,mass,exclude)
                md.system_add_cluster(nodes,cells,r2[0],r2[1],0,mass,exclude)

    else:
        xy = [[L/4,L/2], [L/2,L/4],[3*L/4,L/2],[L/2,3*L/4]]
        for x,y in xy:
            #PyArg_ParseTuple(pArgs, "OOddidi", &vert, &cells, &x, &y, &typ, &mass, &exc)
            print("adding cluster at %f %f %f %f\n" % (x,y,dx,dy))
            md.system_add_cluster(nodes,cells,x,y,0,mass,exclude)

    sigma = 1
    return L, R + sigma/2 


class Saver:
    #def __init__(self, type, filename, pos, vel=None, tetras=None, bonds=None):
    def __init__(self, type, filename, pos, **kwargs):
        '''
        self.n_particles = pos.shape[0]
        self.pos = pos
        self.type = type
        self.vel = vel
        self.tetras = tetras
        self.bonds = bonds
        self.ntetras = 0
        self.nbonds = 0
        if self.tetras is None:
            pass
        else:
            self.ntetras = tetras.shape[0]
        if self.bonds is None:
            pass
        else:
            self.nbonds = bonds.shape[0]
        '''
        self.n_particles = pos.shape[0]
        self.pos = pos
        self.type = type
        self.ncells = 0
        self.nbonds = 0
        for key in kwargs:
            if key == 'vel':
                self.vel = kwargs[key]
            elif key == 'cells':
                self.cells = kwargs[key]
                self.ncells = self.cells.shape[0]
            elif key == 'bonds':
                self.bonds = kwargs[key]   
                self.nbonds = self.bonds.shape[0]
            elif key == 'box':
                self.box = kwargs[key]

        if(type==0):
            self.filename = filename + '.xyz'
            self.file = open(self.filename,'w')
            self.file.close()
        if(type==1):
            self.filename = filename
        if(type==2):
            self.filename = filename + '.vtf'
            typName = ['A', 'B', 'C', 'D', 'E', 'F']
            self.file = open(self.filename,'w')
            for i,p in enumerate(self.pos):
                self.file.write('atom %d type %s radius %f\n' % (i, 'A', 0.5))
            if self.nbonds > 0:
                for b in self.bonds:
                    self.file.write("bond %d:%d\n" % (b[0], b[1]))
            print("pbc %f %f %f\n"% (self.box[0], self.box[1], 1.0));
            self.file.close()


    def append(self,steps):
        if self.type == 0:
            self.savexyz()
        if self.type == 1:
            self.savevtk(steps)
        if self.type == 2:
            self.savevtf()

    def savexyz(self):
        self.file = open(self.filename,'a')
        self.file.write("%d\n\n" % self.n_particles)
        for p in self.pos:
            self.file.write('A %f %f %f\n' % (p[0], p[1], 0.0))
            #print(p)
        self.file.close()

    def savevtk(self,steps):
        #if n_clusters == 0: return
        filename = "%s_%d%s" % (self.filename,steps,'.vtk')
        self.file = open(filename,'w')
        f = self.file
        f.write("# vtk DataFile Version 3.1\n")
        f.write("Really cool data\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n\n")
        f.write("POINTS %d FLOAT\n" % self.n_particles)
        for p in self.pos: f.write("%f %f %f\n" % (p[0], p[1], 0.0)) 
        f.write("\n")
        if self.ncells > 0:
            f.write("CELLS %d %d\n" % (self.ncells, 4*self.ncells))
            for tri in self.cells:
                f.write("%d %d %d %d\n" % (3, tri[0], tri[1], tri[2]))
            f.write("\n")
            f.write("CELL_TYPES %d\n" % self.ncells)
            VTK_TRIANGLE = 5
            VTK_TETRA = 10
            for nt in range(self.ncells): f.write("%d " % VTK_TRIANGLE)
            f.write("\n\n");
            f.write("POINT_DATA %d\n" % self.n_particles)
            f.write("SCALARS VelocitySquare FLOAT\n")
            f.write("LOOKUP_TABLE default\n")
            #vv = (v*v).sum(1)
            for v in self.vel:
                #f.write("%f\n" % v)
                f.write("%f\n" % np.dot(v,v))
            f.write("\n\n")
            f.write("VECTORS Velocity FLOAT\n")
            for v in self.vel:
                f.write("%f %f %f\n" % (v[0],v[1],0.0))
            f.write("\n\n")
            #for v in self.vel:
            #    f.write("%f % %f\n" % (v[0],v[1],0.0))
            #f.write("\n\n")
            #f.write("POINT_DATA %d\n", self.n_particles);
            #f.write("\n");
            #f.write("box size  %f\t%f\t%f\n", box.x, box.y,1.0);
        
        f.close()

    def savevtf(self):
        self.file = open(self.filename,'a')
        self.file.write("timestep\n")
        for p in self.pos:
            self.file.write('%f %f %f\n' % (p[0], p[1], 0.0))
            #print(p)
        self.file.close()

'''
oid MDSystem::saveVTF(const char* fileName, int s) {
	FILE *dataFile;
	char typName[] = { 'A', 'B', 'C', 'D', 'E', 'F' };
	int i;
	if (s == 0){
		dataFile = fopen(fileName, "w");
		for (int i = 0; i < N; i++){
			fprintf(dataFile, "atom %d type %c radius %f\n", i, typName[particles.type[i]], particles.radius[i]); //FIXME change %c to %s
		}
		if (particles.bonds.bondList.size() > 0){
			for (int b = 0; b < particles.bonds.numBonds; b++){
				fprintf(dataFile, "bond %d:%d\n", particles.bonds.bondList[2*b], particles.bonds.bondList[2*b+1]);
			}
		}
		fprintf(dataFile, "pbc %f %f %f\n", box.x, box.y, 1.0);
	}
	else
		dataFile = fopen(fileName, "a");
	fprintf(dataFile, "timestep\n");
	for (i = 0; i < N; i++)
		fprintf(dataFile, "%f %f %f\n", particles.pos[i].x, particles.pos[i].y, 0.0f);
	fclose(dataFile);
}
'''

