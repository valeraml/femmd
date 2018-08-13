
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "particles.h"
#include "md3dsystem.h"


void Particles::init(MDSystem *s){
	sys = s;
	n_particles = 0;
	numTypes = 1;
	maxTypes = 10;
}

void Particles::reset(){

	pos.clear();
	vel.clear();
	acc.clear();
	force.clear();
	vv.clear();

	radius.clear();
	sigma.clear();
	mass.clear();
	type.clear();
	charge.clear();

	enei.clear();
	viri.clear();

	n_particles = 0;
	numTypes = 1;
	maxTypes = 10;
}

int Particles::addParticle(VectorR r0, real m, real sig, int typ) {

	int pid = pos.size();
	n_particles += 1;
	sys->density = sys->particles.n_particles / (sys->box.prod());
	VectorR t;
	t = r0;
	pos.push_back(t);
	t.fill(0.0);
	vel.push_back(t);
	acc.push_back(t);
	force.push_back(t);
	vv.push_back(0);
	mass.push_back(m);
	type.push_back(typ);
	sigma.push_back(sig);
	radius.push_back(sig / 2.0);
	charge.push_back(1);
	sys->exclusion_groups.id_of_particle.push_back(-1);
	sys->clusters.cluster_id_of_particle.push_back(-1);
	return pid;
}

int Particles::addParticle(VectorR r0, VectorR v0, real m, real sig, int typ) {
	int pid = addParticle(r0, m, sig, typ);
	vel[pid] = v0;
	return pid;
}

void Particles::setParticlesVel(real v0){
	for (int p = 0; p < n_particles; p++) {
		vel[p].setRandom()*v0 / 2;
		//vel[p].x = (RANDOM01 - .5f)*v0;
		//vel[p].y = (RANDOM01 - .5f)*v0;
		//real s = RANDOM01;
		//s = 2.0*PI*s;
		//vel[p].x = cos(s)*v0;
		//vel[p].y = sin(s)*v0;
		//printf("%f %f %f\n", v0, vel[p].x, vel[p].y);
		//vel[p].x = 1.0;
		//vel[p].y = 1.0;

	}
}

void Particles::zeroCM(){
	int p;
	VectorR vSum;
	real totMass = 0;;
	vSum.fill(0.0);
	for (p = 0; p < n_particles; p++){
		vSum += mass[p] * vel[p];
		totMass += mass[p];
	}
	// with zero total momentum
	for (p = 0; p < n_particles; p++){
		vel[p] -= vSum / totMass;
	}	
}

void Particles::setSigMax(){
	sigMax = sys->particles.sigma[0];
}


void Particles::checkOverlap() {
	for (int i = 0; i<n_particles; i++) {
		for (int j = 0; j<n_particles; j++) {
			//console.log(this.N,i,j);
			//int gi = sys->particles.exclusionGroupIdOfPart[i];
			//int gj = sys->particles.exclusionGroupIdOfPart[j];
			//if ((i != j) && (gi == -1 || gj == -1 || (gi != gj))){
			if (i != j) {
				real sigmaij = (sigma[i] + sigma[j]) / 2;

				Box L = sys->box;
				//real dx = pos[i].x - pos[j].x;
				//real dy = pos[i].y - pos[j].y;
				VectorR dr = pos[i] - pos[j];
				//if (sys->pbc) {
				//if (0) {
				//	if (dx > 0.5*L.x) dx -= L.x;
				//	else if (dx < -0.5*L.x) dx += L.x;
				//	if (dy > 0.5*L.y) dy -= L.y;
				//	else if (dy < -0.5*dy) dy += L.y;
				//};
				real dr2 = dr.dot(dr);
				real rCut2 = (0.95*sigmaij)*(0.95*sigmaij);
				if (dr2 < rCut2) {
					printf("overlap: %d %d %f\n", i, j, dr2);
				}
			}
		}
	}
}

