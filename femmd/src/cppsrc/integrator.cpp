///////// Integration functions //////////////////////////

#include <ctime>
#include "particles.h"
#include "md3dsystem.h"
#include "integrator.h"

void Integrator::doStep(){

	sys->steps += 1;
	sys->simTime = sys->steps*sys->dt;
	integratorStep1();
	if (sys->useNN)
		sys->neighborlist.update(0);
	sys->interactions.calcForces();
	integratorStep2();
	sys->props.compute();
	//sys->particles.checkOverlap();

}

void Integrator::run(int steps)
{
	if (sys->useNN)
		sys->neighborlist.update(1);
	zeroStuff();
	sys->interactions.calcForces();
	updateAcc();
	clock_t start1;
	double duration;
	start1 = clock();
	for (int i = 0; i < steps; i++){
		doStep();
	}
	duration = (clock() - start1) / (double)CLOCKS_PER_SEC;
	double tps = (double)steps / duration;
	printf("steps: %d time: %f; %f tps\n", steps, duration, tps);
	if(sys->useNN) printf("tnn updates %d %f\n", sys->neighborlist.totalUpdates, (real)steps / sys->neighborlist.totalUpdates);

}
/*
void setThermostat(real g){

	//langevin_pref1 = -langevin_gamma/time_step;
	//langevin_pref2 = sqrt(24.0*temperature*langevin_gamma/time_step);
	//p->f.f[j] = langevin_pref1*p->m.v[j]*PMASS(*p) + langevin_pref2*(d_random()-0.5)*massf;

	printf("setting thermostat, gamma=%f\n", g);
	gamma = g;

	c1 = -gamma;
	c2 = sqrt(24.0*sys->temperature*gamma / sys->dt);
}

void ApplyLangevin(){
	int n;
	real c1, c2;
	c1 = sys->thermostat.c1;
	c2 = sys->thermostat.c2;
	for (n = 0; n<sys->n_particles; n++){
		//cout << GaussianRandom() << endl;
		//mol[n].ra.x += 2*sysp.langevinNoise*(GaussianRandom()-0.5);
		//mol[n].ra.y += 2*sysp.langevinNoise*(GaussianRandom()-0.5);
		//mol[n].ra.z += 2*sysp.langevinNoise*(GaussianRandom()-0.5);

		//pa[n].x += 2*sys.langevinNoise*(genrand_real2()-0.5);
		//pa[n].y += 2*sys.langevinNoise*(genrand_real2()-0.5);
		//pa[n].z += 2*sys.langevinNoise*(genrand_real2()-0.5);

		//sys->pf[n].x += c2*(genrand_real1()-0.5) + c1 * sys->pv[n].x;
		//sys->pf[n].y += c2*(genrand_real1()-0.5) + c1 * sys->pv[n].y;
		//sys->pf[n].z += c2*(genrand_real1()-0.5) + c1 * sys->pv[n].z;

		sys->pf[n].x += c2*(d_random() - 0.5) + c1 * sys->pv[n].x;
		sys->pf[n].y += c2*(d_random() - 0.5) + c1 * sys->pv[n].y;
		sys->pf[n].z += c2*(d_random() - 0.5) + c1 * sys->pv[n].z;

	}
}
*/
/////////////////Leap Frog Steps////////////////

void Integrator::leapFrogStep1(){
	real dt = sys->dt;
	real hdt = 0.5*dt;
	//#pragma omp parallel for
	for (int i = 0; i < sys->particles.n_particles; i++){
		sys->particles.vel[i] += hdt*sys->particles.acc[i];
		sys->particles.pos[i] += dt*sys->particles.vel[i];
	}
}

void Integrator::leapFrogStep2(){
	real hdt = 0.5*sys->dt;
	//#pragma omp parallel for
	for (int i = 0; i < sys->particles.n_particles; i++){
		sys->particles.vel[i] += hdt*sys->particles.acc[i];
	}
}

void Integrator::applyPBC(){
	//#pragma omp parallel for
	for (int i = 0; i < sys->particles.n_particles; i++){
		applyBoundaryCondition(sys->particles.pos[i], sys->box, sys->pbcType);
	}
}

void Integrator::updateAcc(){
	//#pragma omp parallel for
	for (int i = 0; i < sys->particles.n_particles; i++){
		sys->particles.acc[i] = sys->particles.force[i] / sys->particles.mass[i];
	}
}

void Integrator::rescaleVelocities(real vFac){
	//#pragma omp parallel for
	for (int i = 0; i<sys->particles.n_particles; i++){
		sys->particles.vel[i] *= vFac;
	}
}


void Integrator::zeroStuff(){

	sys->zeroCurrVals();
	//#pragma omp parallel for
	for (int i = 0; i<sys->particles.n_particles; i++){
		sys->particles.acc[i].setZero();
		sys->particles.force[i].setZero();
	}
}

void Integrator::integratorStep1(){
		leapFrogStep1();
		applyPBC();
		zeroStuff();
}

void Integrator::integratorStep2(){
		updateAcc();
		leapFrogStep2();
}