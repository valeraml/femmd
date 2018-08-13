#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include "def.h"
#include "groups.h"

class MDSystem;

class Clusters {
public:
	MDSystem *sys;
	int initial_offset;
	int n_clusters;
	hvector<int> particles_offset_index;
	hvector<int> cells_offset_index;
	hvector<int> cluster_id_of_particle;
	hvector<int> cluster_num_of_particles;
	hvector<int> cluster_num_of_cells;
	hvector<real> cluster_mass;
	hvector<real> cluster_density;
	hvector<VectorR> center_of_mass;
	hvector<VectorR> center_of_mass_vel;
	Groups cluster_groups;

	void init(MDSystem *s) { 
		sys = s; 
		n_clusters = 0;
		initial_offset = -1; 
	}
	void reset() {
		n_clusters = 0;
		initial_offset = -1;
	}

	void addCluster(real *vert, int nvert, int *tet, int ntet, VectorR r, int typ, real M, bool exclude);
	int createCluster();
	void addParticlesToCluster(int cluster_index, real *vert, int nvert, VectorR r, int typ, bool exclude);
	void addCellsToCluster(int ind, int *cells, int ncells, real M);
};
#endif