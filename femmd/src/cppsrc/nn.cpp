
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "nn.h"
#include "md3dsystem.h"

int cellId(VectorR &r, VectorR &box, VectorI &numCellsPerDim){
	int ix, iy, iz, cellId;
	//FIXME the ifs are a hack, if x = -0, then bc sends then to box size and cell is out of range
	ix = r.x()*numCellsPerDim.x() / box.x(); if (ix == numCellsPerDim.x()) ix -= 1;
	iy = r.y()*numCellsPerDim.y() / box.y(); if (iy == numCellsPerDim.y()) iy -= 1;
	cellId = ix + iy*numCellsPerDim.x();
	if (NDIM == 3) {
		iz = r[2]*numCellsPerDim[2] / box[2]; if (iz == numCellsPerDim[2]) iz -= 1;
		cellId += iz*numCellsPerDim[0] * numCellsPerDim[1];
	}
	return cellId;
}

void ind2sub(int cellId, VectorI &cellSubscripts, VectorI &numCellsPerDim){
//void cell2DIndex(int cellId, VectorI &cell2DId, VectorI &numCellsPerDim){
	if (NDIM == 3) {
		//cell2DId.z() = cellId / numCellsPerDim.y()*numCellsPerDim.x();
		//cellId -= cell2DId.z()*numCellsPerDim.y()*numCellsPerDim.x();
		cellSubscripts[2] = cellId / (numCellsPerDim[1]*numCellsPerDim[0]);
		cellId -= cellSubscripts[2]*numCellsPerDim[1]*numCellsPerDim[0];
	}
	cellSubscripts.y() = cellId / numCellsPerDim.x();
	cellId -= cellSubscripts.y() * numCellsPerDim.x();
	cellSubscripts.x() = cellId;
}

int sub2ind(VectorI &cellSubscripts, VectorI &numCellsPerDim) {
	int ind;
	//neighborCellId = neighborCell.x() + neighborCell.y()*numCellsPerDim.x();
	ind = cellSubscripts.x() + cellSubscripts.y()*numCellsPerDim.x();
	if (NDIM == 3) {
		//ind += cellSubscripts.z() * numCellsPerDim.x() * numCellsPerDim.y();
		ind += cellSubscripts[2] * numCellsPerDim[0] * numCellsPerDim[1];
	}
	return ind;
}


void cellBCDim(VectorI &curCell, VectorI &numCellsPerDim, int idim) {
	if (curCell[idim] < 0)
		curCell[idim] = numCellsPerDim[idim] - 1;
	else if (curCell[idim] == numCellsPerDim[idim])
		curCell[idim] = 0;
}

int applyCellBC(VectorI &curCell, VectorI &numCellsPerDim, PBCTYPE pbcType){

	int re = 0;
	switch (pbcType){
	case XYZPBC:
		cellBCDim(curCell, numCellsPerDim, 0);
		cellBCDim(curCell, numCellsPerDim, 1);
		cellBCDim(curCell, numCellsPerDim, 2);
		break;
	case XYPBC:
		cellBCDim(curCell, numCellsPerDim, 0);
		cellBCDim(curCell, numCellsPerDim, 1);
		if(NDIM == 3)
			if (curCell[2] < 0 || curCell[2] == numCellsPerDim[2]) re = 1;
	break;
	case XPBC:
		cellBCDim(curCell, numCellsPerDim, 0);
		if (curCell.y() < 0 || curCell.y() == numCellsPerDim.y()) re = 1;
	break;
	case NOPBC:
		if (curCell.x() < 0 || curCell.x() == numCellsPerDim.x()) re = 1;
		if (curCell.y() < 0 || curCell.y() == numCellsPerDim.y()) re = 1;
		if (NDIM == 3)
			if (curCell[2] < 0 || curCell[2] == numCellsPerDim[2]) re = 1;

	break;
	default:
	break;
	}
	//printf("cell bc %d %d\n", curCell.x(), curCell.y());
	return re;
}


void updateCellsCPU(hvector<VectorR> &r, VectorR &box, VectorI &numCellsPerDim, int n_particles, int maxParticlesPerCell,
				hvector<int> &cellCount, hvector<int> &cellContent){
	int cId, particlesInCellcId;
	//#pragma omp parallel for
	for (int pi = 0; pi< n_particles; pi++){
		cId = cellId(r[pi], box, numCellsPerDim);
		particlesInCellcId = cellCount[cId];
		cellCount[cId]++;
		int ind = cId*maxParticlesPerCell + particlesInCellcId;
		//#pragma omp critical
		cellContent[ind] = pi;
	}
}

void buildNN(hvector<VectorR> &r, 
					  int n_particles,
	                  VectorR &box, 
	                  PBCTYPE pbcType, 
					  hvector<int> &cellCount, 
					  hvector<int> &cellContent,
					  VectorI &numCellsPerDim,
					  int maxParticlesPerCell,  
					  hvector<int> &neighborList, 
					  hvector<int> &neighborsCount,
					  int maxNeighborsPerParticle, 
					  real nnRadiusSq,
					  hvector<int> &exclusionGroupIdOfPart){
	//#pragma omp parallel for
	for (int pi = 0; pi<n_particles; pi++){
		int neighborCellId, npCell, pj;
		VectorR dr;
		real dr2;
		VectorI piCell;
		VectorI neighborCell;
		int piCellId = cellId(r[pi], box, numCellsPerDim);
		ind2sub(piCellId, piCell, numCellsPerDim);
		int neighbors = 0;
		for (int ci = 0; ci < NeighborList::nn_iter.size(); ci++) {
			neighborCell = piCell + NeighborList::nn_iter[ci];
			if (applyCellBC(neighborCell, numCellsPerDim, pbcType) == 1)
				continue;
			neighborCellId = sub2ind(neighborCell, numCellsPerDim);
			npCell = cellCount[neighborCellId];
			for (int j = 0; j < npCell; j++) {
				pj = cellContent[neighborCellId*maxParticlesPerCell + j];
				int gi = exclusionGroupIdOfPart[pi];
				int gj = exclusionGroupIdOfPart[pj];
				if ((pi != pj) && (gi == -1 || gj == -1 || (gi != gj))) {
					dr = r[pi] - r[pj];
					nearestImage(dr, box, pbcType);
					dr2 = dr.squaredNorm();
					if (dr2 <= nnRadiusSq) {
						neighborList[pi*maxNeighborsPerParticle + neighbors] = pj;
						neighbors++;
					}
				}
			}
		}
		neighborsCount[pi] = neighbors;
	}
}



void buildNNAll(hvector<VectorR> &r,
	int n_particles,
	VectorR &box,
	PBCTYPE pbcType,
	hvector<int> &cellCount,
	hvector<int> &cellContent,
	VectorI &numCellsPerDim,
	int maxParticlesPerCell,
	hvector<int> &neighborList,
	hvector<int> &neighborsCount,
	int maxNeighborsPerParticle,
	real nnRadiusSq,
	hvector<int> &exclusionGroupIdOfPart) {

	VectorR dr;
	real dr2;
	int neighbors;
	for (int pi = 0; pi<n_particles; pi++){
		neighbors = 0;
		for (int pj = 0; pj<n_particles; pj++){
			int gi = exclusionGroupIdOfPart[pi];
			int gj = exclusionGroupIdOfPart[pj];
			if ((pi != pj) && (gi == -1 || gj == -1 || (gi != gj))){
				dr = r[pi] - r[pj];
				nearestImage(dr, box, pbcType);
				dr2 = dr.squaredNorm();
				if (dr2 < nnRadiusSq){
					neighborList[pi*maxNeighborsPerParticle + neighbors] = pj;
					neighbors++;
				}
			}
		}
		neighborsCount[pi] = neighbors;
	}
}


void NeighborList::buildNN(){
	
	memset(&cellCount[0], 0, cellCount.size() * sizeof(cellCount[0]));
	memset(&cellContent[0], 0, cellContent.size() * sizeof(cellCount[0]));
	//memset(&neighborList[0], 0, neighborList.size() * sizeof(neighborList[0]));
	//memset(&neighborCount[0], 0, neighborCount.size() * sizeof(neighborCount[0]));

	
	updateCellsCPU(sys->particles.pos, 
		           sys->box, 
				   numCellsPerDim,
		           sys->particles.n_particles,
				   maxParticlesPerCell, 
		           cellCount, 
		           cellContent);

	
	::buildNN(sys->particles.pos, 
					 sys->particles.n_particles,
		             sys->box, 
		             sys->pbcType, 
		             cellCount, 
		             cellContent, 
		             numCellsPerDim,
		             maxParticlesPerCell,  
					 neighborList,
					 neighborCount,
					 maxNeighborsPerParticle, 
					 nnRadiusSq,
					 sys->exclusion_groups.id_of_particle);
	
    /*
	buildNNAll(sys->particles.pos,
		sys->particles.n_particles,
		sys->box,
		sys->pbcType,
		cellCount,
		cellContent,
		numCellsPerDim,
		maxParticlesPerCell,
		neighborList,
		neighborCount,
		maxNeighborsPerParticle,
		nnRadiusSq,
		sys->exclusion_groups.id_of_particle);
	*/
}


void NeighborList::update(int restart){
	//printf("xxx list %d %d %f %f\n", totalUpdates, totalSteps, dispHi, 0.5*skin);
	//buildNN();
	if (restart)initialUpdate = 1;
	if (dispHi > 0.5*skin) {
		totalUpdates++;
		totalSteps = sys->steps - initialSteps;
		//printf("updating list %d %d %d %f %f\n", sys->steps, totalUpdates, totalSteps, dispHi, (real) totalSteps / totalUpdates);
		dispHi = 0.0;
		buildNN();
	}
	if (initialUpdate){
		totalUpdates = 1;
		totalSteps = 0;
		initialSteps = sys->steps;
		initialUpdate = 0;
		dispHi = 0;
		buildNN();
		//printf("initial update list at %d\n",sys->steps);
	}
}


//void NeighborList::buildNN(){
//	buildNN();
//}


void NeighborList::allocate(){

	cellCount.resize(numCells);
	cellContent.resize(maxParticlesPerCell*numCells);
	neighborList.resize(sys->particles.n_particles*maxNeighborsPerParticle);
	neighborCount.resize(sys->particles.n_particles);

	memset(&cellCount[0], 0, cellCount.size() * sizeof(cellCount[0]));
	memset(&cellContent[0], 0, cellContent.size() * sizeof(cellCount[0]));
	memset(&neighborList[0], 0, neighborList.size() * sizeof(neighborList[0]));
	memset(&neighborCount[0], 0, neighborCount.size() * sizeof(neighborCount[0]));

	//cellList.resize(numCells + sys->particles.n_particles);
	//std::fill(cellList.begin(), cellList.end(), -1);

}

void NeighborList::init(){
	//normalize to max sigma;
	//real s = sys->particles.sigMax;
	//printf("sigmaMax %f\n", s);
	real s = 1.0;
	init(s*sys->neighborlist.skin, s*sys->interactions.rCutMax);
}

hvector<VectorI> NeighborList::nn_iter;

void NeighborList::init(real s, real rc){
	//dimGrid = sys->dimGrid;
	//dimBlock = sys->dimBlock;
	N = sys->particles.n_particles;
	skin = s;
	rCut = rc; //FIXME incorporate sigma
	rCutSq = rc*rc;
	nnRadius = (rc + s);
	dispHi = 0;
	totalUpdates = 0;
	initialUpdate = 1;
	nnRadiusSq = nnRadius*nnRadius;
	for(int i=0; i<NDIM; i++)
		numCellsPerDim[i] = int(floor(sys->box[i] / nnRadius));
	real vCell = sys->box.prod() / numCellsPerDim.prod();
	real factor = 2;
	real dens = sys->density;
	if (dens < 1.0) dens = 1.0;
	numCells = numCellsPerDim.prod();
	real temp1 = (factor*dens*vCell);
	maxParticlesPerCell = temp1;
	//(approx 4/3*pi ~ 5)
	//maxNeighborsPerParticle = (int)(factor * 5.0*nnRadius*nnRadius*nnRadius*dens);
	real temp = (factor * 5.0*nnRadius*nnRadius*nnRadius*dens);
	maxNeighborsPerParticle = temp;


	int nn_iter_size = (int)std::pow(3, NDIM);
	nn_iter.resize(nn_iter_size);
	int nn_ind = 0;
	if (NDIM == 2) {
		for (int cx = -1; cx <= 1; cx++) {
			for (int cy = -1; cy <= 1; cy++) {
				nn_iter[nn_ind][0] = cx;
				nn_iter[nn_ind][1] = cy;
				nn_ind += 1;
			}
		}
	}
	else if (NDIM == 3) {
		for (int cx = -1; cx <= 1; cx++) {
			for (int cy = -1; cy <= 1; cy++) {
				for (int cz = -1; cz <= 1; cz++) {
					nn_iter[nn_ind][0] = cx;
					nn_iter[nn_ind][1] = cy;
					nn_iter[nn_ind][2] = cz;
					nn_ind += 1;
				}
			}
		}
	}
	
	printf("Neighbor List Initialized: %d %d %d\n", maxParticlesPerCell, maxNeighborsPerParticle, numCells);

	allocate();
};

void NeighborList::reset(){
	cellCount.clear();
	cellContent.clear();
	neighborList.clear();
	neighborCount.clear();

}
