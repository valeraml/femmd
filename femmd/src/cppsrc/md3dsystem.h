#ifndef __MD3DSYSTEM_H__
#define __MD3DSYSTEM_H__

#include <vector>
#include "def.h"
#include "particles.h"
#include "bonds.h"
#include "nn.h"
#include "interactions.h"
#include "integrator.h"
#include "properties.h"
#include "groups.h"

#define PYBIND11_INTERFACE
//#define CPYTHON_INTERFACE

// Properties

class MDSystem;

class Walls {
public:
	MDSystem *sys;
	hvector<VectorR> forces;
	hvector<VectorR> pos;
	hvector<VectorR> motion_rate;
	real dt;
	int moveWalls;
	int shearWalls;

	void init(MDSystem *s, VectorR l0) {
		sys = s;
		VectorR zero;
		zero.setZero();
		forces.push_back(zero); // left, bottom and back forces
		forces.push_back(zero); // right, top, and front forces
		motion_rate.push_back(zero);
		motion_rate.push_back(zero);
		pos.push_back(zero);
		pos.push_back(l0);
		moveWalls = 0;
		shearWalls = 0;
	}

//if (pos > sys->box.x()) right_pos = sys->box.x();
	void move(real dt);
	void shear(real dt);
	void linear_move(real dt);
	void set_motion_rate(VectorR v0, VectorR v1) {
		motion_rate[0] = v0;
		motion_rate[1] = v1;
	}
	void set_walls_postions(VectorR b0, VectorR b1) {
		pos[0] = b0;
		pos[1] = b1;
	}

};


class MDSystem {

public:
	bool start;

	Particles particles;
	Bonds<real> bonds;
	Elements elements;
	Clusters clusters;
	NeighborList neighborlist;
	Interactions interactions;
	Integrator integrator;
	Properties props;
	Groups exclusion_groups;
	Walls walls;

	int steps;
	real dt;
	real simTime;
	real initialTemperature;

	real density;
	VectorR box;
	bool useNN;
	PBCTYPE pbcType;
	real scale;

	real vvMax, clusterKinEneSum, kinEneSum, potEnergy, bondedEnergy, pairEnergy,
		feaEnergy, virial, feaVirial, clusterVirial, temperature;
	real wallPressure;
	int averageSteps;
	void zeroCurrVals(){
		vvMax = kinEneSum = potEnergy = bondedEnergy = 0.0;
		feaEnergy = virial = temperature = wallPressure = 0.0;
		feaVirial = clusterVirial = clusterKinEneSum = 0;
		pairEnergy = 0.0;
	}

	bool adjustTemperature;
	int pairCount;

	MDSystem() {}
	MDSystem(int nn) {
		//n_particles = 0;
		//n_clusters = 0;
	}

	void init(VectorR l0);
	void set_defaults();
	void reset();

	void setBox(VectorR lbox) {
		box = lbox;
		VectorR b0;
		b0.setZero();
		walls.set_walls_postions(b0, box);
	}

	void setTemperature(real temp);
	void rescaleVelocities(real vFac);
	void evalProps();
	void saveXYZ(const char* fileName, int s);

	void printSystemInfo(){
		printf("Number of Particles: %d\n", particles.n_particles);
		printf("delta t: %f\n", dt);
		printf("Number of Elastomers: %d\n", clusters.n_clusters);
		printf("Number of bonds %d\n", (int)bonds.bondList.size());
		printf("Number of Triangles %d\n", (int)elements.tetras.size());
	}

};



#endif