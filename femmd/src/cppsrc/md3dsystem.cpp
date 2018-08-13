// JavaScript Document

/*
	A 2D molecular dynamics simulation in c/c++
	
	Copyright 2015, Manuel Valera
	
	Permission is hereby granted, free of charge, to any person obtaining a copy of 
	this software and associated data and documentation (the "Software"), to deal in 
	the Software without restriction, including without limitation the rights to 
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
	of the Software, and to permit persons to whom the Software is furnished to do 
	so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all 
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR 
	ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
	OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
	OTHER DEALINGS IN THE SOFTWARE.

	Except as contained in this notice, the name of the author shall not be used in 
	advertising or otherwise to promote the sale, use or other dealings in this 
	Software without prior written authorization.
	
	
*/


#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "md3dsystem.h"

real PI;
real one_over_ndim_factorial;

//////////// Class system functions ///////////////////////////

//////////// Inititalization functions///////////////


void MDSystem::init(VectorR lbox){
	std::cout << "size of VectorR " << sizeof(VectorR) << std::endl;
	std::cout << "size of real " << sizeof(real) << std::endl;
	std::cout << "size of MatrixDxD " << sizeof(MatrixDxD) << std::endl;

	PI = 4.0 * atan(1.0);
	if (NDIM == 2)
		one_over_ndim_factorial = 0.5;
	else if (NDIM == 3)
		one_over_ndim_factorial = 1.0 / 6.0;
	start = true;
	steps = 0;
	dt = 0.005;
	simTime = 0;
	averageSteps = 100;
	initialTemperature = 1.0;
	pbcType = XYPBC;

	density = 0;
	useNN = true;

	//Properties 
	bool adjustTemperature = true;

	box = lbox;
	density = particles.n_particles / (box.prod());
	walls.init(this, box);
	
	particles.init(this);
	bonds.init(this);
	elements.init(this);
	clusters.init(this);

	neighborlist.init(this);
	interactions.init(this);
	integrator.init(this);
	props.sys = this;

	particles.numTypes = 1;
}

void MDSystem::set_defaults() {
	useNN = false;
	dt = 0.001;
	props.step_avg = 100;
	props.step_equi = 0;
	particles.numTypes = 1;
	pbcType = XYPBC;
	interactions.gravity = 0;
	interactions.gravityForce = false;
	interactions.E = 5;
	interactions.kBond = 10;
	interactions.bondForces = false;
	interactions.feaForces = false;
	interactions.areaForces = false;
	interactions.setElasticConstants();

	double v0 = 0.5;
	//sys->particles.setParticlesVel(v0);
	//if (sys->useNN)
	//	sys->neighborlist.init();

	start = true;
}

void MDSystem::reset(){
	particles.reset();
	bonds.reset();
	elements.reset();
	neighborlist.reset();
	interactions.reset();
	props.reset();
}



void MDSystem::setTemperature(real temp){

	real vvSum = 0;
	real kinEne = 0;
	real vFac;

	for (int i = 0; i<particles.n_particles; i++){
		int pi = i;
		particles.vv[pi] = particles.vel[pi].squaredNorm();
		vvSum += particles.vv[pi];
		kinEne += 0.5f*this->particles.mass[pi] * particles.vv[pi];
	}
	//var vMag = Math.sqrt(2*temp/mass); 
	//vFac = vMag/Math.sqrt(vvSum/numAtoms);
	real currTemp = kinEne / particles.n_particles;
	vFac = sqrt(temp / currTemp);
	for (int pi = 0; pi<particles.n_particles; pi++){
		particles.vel[pi] *= vFac;
	}
	//particles.zeroCM();
}

void MDSystem::evalProps(){

		for (int pi = 0; pi<particles.n_particles; pi++){
			//particles.vv[pi] = particles.vel[pi].x*particles.vel[pi].x + particles.vel[pi].y*particles.vel[pi].y;
			particles.vv[pi] = particles.vel[pi].squaredNorm();
			vvMax = std::max(vvMax, particles.vv[pi]);
			real kin = particles.mass[pi] * particles.vv[pi];
			//vvSum += particles.vv[pi];
			kinEneSum += kin;
		}
		//if (clusters.data.size()>0){
		//	for (int ci = 0; ci < clusters.data.size(); ci++){
		//		real cvv = clusters.data[ci].cmVel.x*clusters.data[ci].cmVel.x +
		//			clusters.data[ci].cmVel.y*clusters.data[ci].cmVel.y;
		//		clusterKinEneSum += clusters.data[ci].mass* cvv;
		//	}
		//}
	neighborlist.dispHi += sqrt(vvMax) * dt;

}

/*
void MDSystem::saveXYZ(const char* fileName, int s) {
	//FILE *dataFile;
	ofstream dataFile;
	int i;
	if (s == 0)
		//errno_t err = fopen(&dataFile, fileName, "w");
		//dataFile = fopen(fileName, "w");
	else
		//errno_t err = fopen(&dataFile, fileName, "w");
		//dataFile = fopen(fileName, "a");
	fprintf(dataFile, "%d\n", particles.n_particles);
	fprintf(dataFile, "\n");
	//fprintf(dataFile, "box size  %f\t%f\t%f\n", box.x, box.y,1.0);
	for (i = 0; i < particles.n_particles; i++)
		fprintf(dataFile, "A\t%f\t%f\t%f\n", particles.pos[i].x(), particles.pos[i].y(), 0.0f);
	//fprintf(dataFile,"\n");
	fclose(dataFile);
}
*/
void Walls::move(real dt) {
	if (shearWalls) shear(dt);
	else linear_move(dt);
}


void Walls::linear_move(real dt) {
	pos[0] += dt*motion_rate[0];
	for (int i = 0; i<NDIM; i++) if (pos[0][i] < 0) pos[0][i] = 0.0;
	pos[1] += dt*motion_rate[1];
	for (int i = 0; i<NDIM; i++) if (pos[1][i] > sys->box[i]) pos[1][i] = sys->box[i];
}

void Walls::shear(real dt) {
	VectorR L = pos[1] - pos[0];
	pos[0][0] += dt*motion_rate[0][0] * L[1] / L[0];
	pos[0][1] += dt*motion_rate[0][1];
	pos[1][0] += dt*motion_rate[1][0] * L[1] / L[0];
	pos[1][1] += dt*motion_rate[1][1];
	pos[1] += dt*motion_rate[1];
	for (int i = 0; i < NDIM; i++) {
		if (pos[0][i] < 0) pos[0][i] = 0.0;
		if (pos[1][i] > sys->box[i]) pos[1][i] = sys->box[i];
	}
}


