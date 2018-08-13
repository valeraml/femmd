#include <vector>
#include <algorithm>
#include "md3dsystem.h"
#include "Elements.h"



void Elements::init(MDSystem *s){
	sys = s;
	sys->elements.totRefVol = 0;
	numTetras = 0;
	//offset = 0; //FIXME
	//offset = sys->clusters.offset;
	//printf("Num of triangles: %d\n", numTetras);
}

void Elements::reset(){

	refPos.clear();
	unfoldedPos.clear();
	tetras.clear();
	//xm.clear();
	refVol.clear();
	currVol.clear();

	numTetras = 0;
	offset = 0;
	tetraListLen = 0;
	tetraNN.clear();
	offset = 0;
}

void Elements::checkAreas(VectorR box) {
	VectorR r[3], r0[3];

	for (int it = 0; it < tetras.size(); it++){
		for (int i = 0; i < 3; i++) {
			int pi = tetras[it][i];
			r[i] = sys->particles.pos[pi];
			//r[i].y() = sys->particles.pos[pi].y();
			r0[i] = refPos[pi];
			//r0[i].y() = refPos[pi].y();
		}
		VectorR r01 = r[0] - r[1];
		VectorR r02 = r[0] - r[2];

		real dx01 = r[0].x() - r[1].x();
		if (dx01 > box.x() / 2.0) r[1].x() += box.x();
		if (dx01 < -box.x() / 2.0) r[1].x() -= box.x();

		real dx02 = r[0].x() - r[2].x();
		if (dx02 > box.x() / 2.0) r[2].x() += box.x();
		if (dx02 < -box.x() / 2.0) r[2].x() -= box.x();

		real dy01 = r[0].y() - r[1].y();
		if (dy01 > box.y() / 2.0) r[1].y() += box.y();
		if (dy01 < -box.y() / 2.0) r[1].y() -= box.y();

		real dy02 = r[0].y() - r[2].y();
		if (dy02 > box.y() / 2.0) r[2].y() += box.y();
		if (dy02 < -box.y() / 2.0) r[2].y() -= box.y();

		real vol = 0.5*((r0[2].y() - r0[0].y())*(r0[1].x() - r0[0].x()) - (r0[1].y() - r0[0].y())*(r0[2].x() - r0[0].x()));
		real newVol = 0.5*((r[2].y() - r[0].y())*(r[1].x() - r[0].x()) - (r[1].y() - r[0].y())*(r[2].x() - r[0].x()));
		if (vol*newVol < 0)printf("Error in triangle %d %f %f \n", it, newVol, vol);
	}
}

void Elements::unfoldPos(VectorR box){

	VectorR rmin; 
	rmin.fill(100000);
	VectorR rmax;
	rmax.fill(-100000);
	VectorR rc;
	rc.setZero();

	for (int ei = 0; ei < sys->clusters.n_clusters; ei++) {
		int eoffset = sys->clusters.particles_offset_index[ei];
		int esize = sys->clusters.cluster_num_of_particles[ei];
		for (int pi = 0; pi < esize; pi++){
			unfoldedPos[pi + eoffset] = sys->particles.pos[pi + eoffset];
			rc += unfoldedPos[pi + eoffset];
			rmin = rmin.cwiseMin(unfoldedPos[pi + eoffset]);
			rmax = rmax.cwiseMax(unfoldedPos[pi + eoffset]);
		}
		if (rmax[0] - rmin[0] > box.x() / 2){
			rc.x() = 0;
			for (int pi = 0; pi < esize; pi++){
				if (unfoldedPos[pi + eoffset].x() > box.x() / 2){
					unfoldedPos[pi + eoffset].x() -= box.x();
				}
				rc.x() += unfoldedPos[pi + eoffset].x();
			}
			if (rc.x() / esize < 0){
				for (int pi = 0; pi < esize; pi++){
					unfoldedPos[pi + eoffset].x() += box.x();
				}
			}
		}
		if (rmax[1] - rmin[1] > box.y() / 2){
			rc.y() = 0;
			for (int pi = 0; pi < esize; pi++){
				if (unfoldedPos[pi + eoffset].y() > box.y() / 2){
					unfoldedPos[pi + eoffset].y() -= box.y();
				}
				rc.y() += unfoldedPos[pi + eoffset].y();
			}
			if (rc.y() / esize < 0){
				for (int pi = 0; pi < esize; pi++){
					unfoldedPos[pi + eoffset].y() += box.y();
				}
			}
		}
		if (NDIM == 3) {
			if (rmax[2] - rmin[2] > box[2] / 2) {
				rc[2] = 0;
				for (int pi = 0; pi < esize; pi++) {
					if (unfoldedPos[pi + eoffset][2] > box[2] / 2) {
						unfoldedPos[pi + eoffset][2] -= box[2];
					}
					rc[2] += unfoldedPos[pi + eoffset][2];
				}
				if (rc[2] / esize < 0) {
					for (int pi = 0; pi < esize; pi++) {
						unfoldedPos[pi + eoffset][2] += box[2];
					}
				}
			}
		}
		
		VectorR vc;
		vc.setZero();
		rc.setZero();
		for (int pi = 0; pi < esize; pi++){
			real mpi = sys->particles.mass[pi + eoffset];
			rc += mpi*unfoldedPos[pi + eoffset];
			vc += mpi*sys->particles.vel[pi + eoffset];
		}
		//FIXME
		sys->clusters.center_of_mass[ei] = rc / sys->clusters.cluster_mass[ei];
		sys->clusters.center_of_mass_vel[ei] = vc / sys->clusters.cluster_mass[ei];
	}
}

// calculate cm using
// https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
// or  http://lammps.sandia.gov/threads/msg44589.html
// compute 1 all property / atom xs
// variable xi atom lx / 2.0 / PI*cos(c_1*2.0*PI)
// variable zeta atom lx / 2.0 / PI*sin(c_1*2.0*PI)
// compute 2 all reduce ave v_xi v_zeta
// variable xb equal lx / 2.0 / PI*(atan2(-c_2[2], -c_2[1]) + PI)

/*
void Elements::updateClustersCentroid(VectorR &box){
	real twopi = 2.0*PI;
	for (int ei = 0; ei < sys->clusters.data.size(); ei++){
		real cxavg = 0;
		real sxavg = 0;
		real cyavg = 0;
		real syavg = 0;
		real x, y, tx, ty;
		int eoffset = sys->clusters.data[ei].offset;
		int esize = sys->clusters.data[ei].nvertices;
		for (int pi = 0; pi < esize; pi++){
			x = particles->pos[pi + eoffset].x;
			tx = twopi * x / box.x;
			real ci = cos(tx);
			real si = sin(tx);
			cxavg += ci;
			sxavg += si;

			y = particles->pos[pi + eoffset].y;
			ty = twopi * y / box.y;
			ci = cos(ty);
			si = sin(ty);
			cyavg += ci;
			syavg += si;

		}
		cxavg = cxavg / esize;
		sxavg = sxavg / esize;
		real txbar = atan2(-sxavg, -cxavg) + PI;
		sys->clusters.data[ei].centroid.x = box.x * txbar / twopi;
		cyavg = cxavg / esize;
		syavg = sxavg / esize;
		real tybar = atan2(-syavg, -cyavg) + PI;
		sys->clusters.data[ei].centroid.y = box.y * tybar / twopi;
	}
}
*/
/*
void Elements::calcClusterProps(){
	if (sys->clusters.data.size() > 0){
		//unfoldPos(particles->sys->box);
		updateClustersCentroid(particles->sys->box);
	}
}
*/
