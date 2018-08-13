#include "md3dsystem.h"
#include "clusters.h"
#include <iostream>



void Clusters::addCluster(real *vert, int nvert, int *tet, int ntet, VectorR r, int typ, real mass, bool exclude) {
	int cluster_index = createCluster();
	addParticlesToCluster(cluster_index, vert, nvert, r, typ, exclude);
	addCellsToCluster(cluster_index, tet, ntet, mass);
}

int Clusters::createCluster() {
	int ind = n_clusters;
	n_clusters += 1;
	cluster_mass.push_back(-1);
	cluster_density.push_back(-1);
	VectorR t;
	t.setZero();
	center_of_mass.push_back(t);
	center_of_mass_vel.push_back(t);
	cluster_groups.addGroup();
	return ind;
}

void Clusters::addParticlesToCluster(int cluster_index, real *vert, int nvert, VectorR r, int typ, bool exclude) {
	int pid;
	int particles_offset = sys->particles.n_particles;
	particles_offset_index.push_back(particles_offset);
	cluster_num_of_particles.push_back(nvert);
	int exclusion_group_id = sys->exclusion_groups.addGroup();
	hvector<VectorR> nodes;
	for (int i = 0; i < nvert; i++) {
		VectorR node;
		for (int j = 0; j < NDIM; j++) {
			node[j] = vert[NDIM*i + j];
		}
		nodes.push_back(node);
	}
	for (int i = 0; i < nvert; i++) {
		VectorR rp = nodes[i] + r;
		pid = sys->particles.addParticle(rp, 0.0, 1.0, typ); //Create particle with 0 mass and sigma = 1
		applyBoundaryCondition(sys->particles.pos[pid], sys->box, sys->pbcType);
		if (exclude) sys->exclusion_groups.insertParticleInGroup(pid, exclusion_group_id); 
		cluster_id_of_particle[pid]=cluster_index;
	}
}

void Clusters::addCellsToCluster(int cluster_index, int *cells, int ncells, real mass) {

	int nvert = cluster_num_of_particles[cluster_index];
	int particles_offset = particles_offset_index[cluster_index];
	int current_number_of_particle = sys->particles.n_particles;
	//sys->elements.refPos.resize((cluster_index + 1)*nvert);
	//sys->elements.unfoldedPos.resize((cluster_index + 1)*nvert); 
	sys->elements.refPos.resize(current_number_of_particle);
	sys->elements.unfoldedPos.resize(current_number_of_particle);
	for (int i = 0; i < nvert; i++) {
		sys->elements.refPos[i + particles_offset] = sys->particles.pos[i + particles_offset];
		sys->elements.unfoldedPos[i + particles_offset] = sys->particles.pos[i + particles_offset];
	}
	//Add elements
	int nb = 0;
	//int cells_offset = cluster_index*ncells;
	int current_number_of_cells = sys->elements.tetras.size();
	int cells_offset = current_number_of_cells;
	int new_number_of_cells = current_number_of_cells + ncells;
	cells_offset_index.push_back(cells_offset);
	cluster_num_of_cells.push_back(ncells);
	//sys->elements.tetras.resize((cluster_index + 1)*ncells);
	sys->elements.tetras.resize(new_number_of_cells);
	sys->elements.Bm.resize(new_number_of_cells);
	sys->elements.refVol.resize(new_number_of_cells);
	sys->elements.currVol.resize(new_number_of_cells);
	sys->elements.I1.resize(new_number_of_cells);
	sys->elements.I2.resize(new_number_of_cells);
	for (int icell = 0; icell < ncells; icell++) {
		//sys->elements.numTetras += 1;
		elementIndexes particle_index;
		VectorR r0[NDIM + 1];
		for (int j = 0; j < NDIM + 1; j++){
			particle_index[j] = cells[(NDIM + 1)*icell + j] + particles_offset;
			sys->elements.tetras[icell + cells_offset][j] = particle_index[j];
			r0[j] = sys->elements.refPos[particle_index[j]];
		}

		for (int j = 0; j < NDIM; j++) {
			for (int k = j + 1; k < NDIM + 1; k++) {
				sys->bonds.createBond(particle_index[j], particle_index[k]);
				nb++;
			}
		}

/*
		real Dm[2][2];
		real Bm[2][2];

		//Agorithm from fendefo.org page 30

		// calculate Dm from ref pos (Xi - X0)
		Dm[0][0] = r0[0].x() - r0[2].x();
		Dm[1][0] = r0[0].y() - r0[2].y();
		Dm[0][1] = r0[1].x() - r0[2].x();
		Dm[1][1] = r0[1].y() - r0[2].y();

		//det and volume of undeformed triangle;
		real detDm = Dm[0][0] * Dm[1][1] - Dm[0][1] * Dm[1][0];
		real vol0 = detDm / 2;

		// Bm inverse of Dm Bm = Inverse(Dm)
		Bm[0][0] = Dm[1][1] / detDm;
		Bm[1][0] = -Dm[1][0] / detDm;
		Bm[0][1] = -Dm[0][1] / detDm;
		Bm[1][1] = Dm[0][0] / detDm;
*/

		VectorR rij[NDIM];
		MatrixDxD Dm, Bm;
		for (int j = 0; j < NDIM; j++) { //loop through vertices 0, 1 .. NDIM in cell i
			//rij[j] = r0[j] - r0[NDIM];
			//Dm.col(j) = rij[j];
			Dm.col(j) = r0[j] - r0[NDIM];
		}
		//real vol = 0.5*((r0[2].y() - r0[0].y())*(r0[1].x() - r0[0].x()) - (r0[1].y() - r0[0].y())*(r0[2].x() - r0[0].x()));
		real detDm = Dm.determinant();
		real vol = detDm * one_over_ndim_factorial;
		sys->elements.refVol[icell + cells_offset] = vol;
		sys->elements.totRefVol += abs(vol);
		Bm = Dm.inverse();
		sys->elements.Bm[icell + cells_offset] = Bm;
		sys->elements.I1[icell + cells_offset] = 0;
		sys->elements.I2[icell + cells_offset] = 0;
		if (vol <= 0) {
			//printf("Error: neg vol: %f\n", vol);
			//std::cout << Bm << "\n\n";
		}

	}
	cluster_mass[cluster_index] = mass;
	cluster_density[cluster_index] = mass /sys->elements.totRefVol;
	sys->elements.numTetras = sys->elements.tetras.size();
	for (int itet = 0; itet < ncells; itet++) {
		for (int j = 0; j < (NDIM+1); j++) {
			int ip = sys->elements.tetras[itet + cells_offset][j];
			sys->particles.mass[ip] += abs( cluster_density[cluster_index]*sys->elements.refVol[itet] )/ (NDIM+1);
		}
	}
	printf("Total Volume: %f\n", sys->elements.totRefVol);

}

