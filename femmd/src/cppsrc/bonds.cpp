#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "md3dsystem.h"
#include "particles.h"
#include "bonds.h"


template <typename NUMTYPE>
void Bonds<NUMTYPE>::initBonds(){
	int N = sys->particles.n_particles;
	maxBondPartnersNN = 300;

	std::vector<int> tempBondList(N*maxBondPartnersNN, 0);
	std::vector<NUMTYPE> tempBondLengthList(N*maxBondPartnersNN, 0);
	std::vector<int> tempBondCount(N, 0);

	for (int b = 0; b < numBonds; b++){
		int i = bondList[2 * b];
		int j = bondList[2 * b + 1];

		int jInd = tempBondCount[i];
		tempBondList[i*maxBondPartnersNN + jInd] = j;
		tempBondLengthList[i*maxBondPartnersNN + jInd] = bondLength[b];
		tempBondCount[i] += 1;
		//printf("x %d (%d %d) %d\n", b, i, j, tempBondCount[i]);
		int iInd = tempBondCount[j];
		tempBondList[j*maxBondPartnersNN + iInd] = i;
		tempBondLengthList[j*maxBondPartnersNN + iInd] = bondLength[b];
		tempBondCount[j] += 1;
		//printf("xx %d (%d %d) %f (%d %d)\n\n", b,i,j,bondLength[b],tempBondCount[i],tempBondCount[j]);
	}
	int oldMaxBondPartners = maxBondPartnersNN;
	maxBondPartnersNN = *max_element(tempBondCount.begin(), tempBondCount.end());
	int sum = std::accumulate(tempBondCount.begin(), tempBondCount.end(), 0);
	printf("total bonds: %d\n", sum);
	bondListNN.resize(N*maxBondPartnersNN);
	bondLengthNN.resize(N*maxBondPartnersNN);
	bondCountNN.resize(N);
	for (int i = 0; i < N; i++){
		bondCountNN[i] = tempBondCount[i];
		for (int j = 0; j < bondCountNN[i]; j++){
			bondListNN[i*maxBondPartnersNN + j] = tempBondList[i*oldMaxBondPartners + j];
			bondLengthNN[i*maxBondPartnersNN + j] = tempBondLengthList[i*oldMaxBondPartners + j];
			//printf("(%d %d) %d->%d %f\n",i,j,i,bondListNN[i*maxBondPartnersNN+j],bondLengthNN[i*maxBondPartnersNN+j]);
		}
	}

}


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




template class Bonds<float>;
template class Bonds<double>;