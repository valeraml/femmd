#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include "def.h"
#include "bonds.h"
#include "Elements.h"

class MDSystem;

class Particles{
public:
	MDSystem *sys;
	int n_particles;

	hvector<VectorR> pos;
	hvector<VectorR> vel;
	hvector<VectorR> acc;
	hvector<VectorR> force;
	hvector<real> vv;

	hvector<real> radius;
	hvector<real> sigma;
	hvector<real> mass;
	hvector<real> charge;
	hvector<int> type;

	int numTypes;
	int maxTypes;
	real sigMax;

	hvector<real> enei;
	hvector<real> viri;

	void init(MDSystem *sys);
	void reset();

	int addParticle(VectorR r0, real m, real sig, int typ);
	int addParticle(VectorR r0, VectorR v0, real m, real sig, int typ);
	void setParticlesVel(real v0);
	void zeroCM();
	void setSigMax();
	void checkOverlap();
};

#endif
