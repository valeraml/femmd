#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "def.h"
#include "clusters.h"

class MDSystem;

class Elements{
public:
	MDSystem *sys;
	//Particles *particles;
	int numTetras;
	int offset;

	hvector< VectorR > refPos;
	hvector< VectorR > unfoldedPos; //FIXME unfolde pos should not be in elements, should in particles
	hvector< elementIndexes > tetras;
	hvector< MatrixDxD > Bm;
	hvector< real > refVol;
	hvector< real > currVol;
	hvector< real > I1, I2;
	real totRefVol;
	real totCurrVol;

	int tetraListLen;
	hvector<elementIndexes> tetraNN;

	void init(MDSystem *sys);
	void reset();
	void unfoldPos(VectorR box);
	void checkAreas(VectorR box);
	void updateClustersCentroid(VectorR &box);
	void calcClusterProps();

};


#endif