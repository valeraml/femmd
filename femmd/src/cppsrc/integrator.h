#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__
#include "def.h"


class MDSystem;
class Particles;
struct curandState;

class Integrator{
public:
	MDSystem *sys;
	void init(MDSystem *s){ sys = s; }

	void doStep();
	void run(int steps);
	void leapFrogStep1();
	void leapFrogStep2();
	void applyPBC();
	void updateAcc();
	void rescaleVelocities(real vFac);
	void zeroStuff();

	void integratorStep1();
	void integratorStep2();

	curandState* devCurandStates;

};

void applyBoundaryCondition(VectorR &r, VectorR &box, PBCTYPE pbcType);


#endif
