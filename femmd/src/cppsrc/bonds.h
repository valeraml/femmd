#ifndef BONDS_H
#define BONDS_H

#include "def.h"

class MDSystem;

template<typename NUMTYPE>
class Bonds{

public:
	MDSystem *sys;
	int numBonds;
	int bondListLen;
	hvector<int> bondList;
	hvector<NUMTYPE> bondLength;

	hvector<int>bondListNN;
	hvector<NUMTYPE>bondLengthNN;
	hvector<int>bondCountNN;
	int maxBondPartnersNN;

	void createBond(int atomi, int atomj);
	void initBonds();
	void init(MDSystem *s) { sys = s; }
	void reset();

};

/*
template <typename NUMTYPE>
void Bonds<NUMTYPE>::createBond(int atomi, int atomj) {
	bondList.push_back(atomi);
	bondList.push_back(atomj);
	VectorR r = sys->particles.pos[atomi] - sys->particles.pos[atomj];
	bondLength.push_back(r.norm());
	numBonds += 1;
	bondListLen += 2;
}

template <typename NUMTYPE>
void Bonds<NUMTYPE>::reset() {
	bondList.clear();
	bondLength.clear();

	bondListNN.clear();
	bondLengthNN.clear();
	bondCountNN.clear();

	numBonds = 0;
	bondListLen = 0;
}
*/

#endif