#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include "particles.h"
#include "nn.h"
#include "md3dsystem.h"

void Interactions::init(MDSystem *s){
	sys = s;
	rCutMax = pow(2.0, 1.0 / 6.0); 
	kWall = 100;
	bondForces = false;
	feaForces = false;
	gravityForce = false;
	frictionForce = false;
	elasticForceType = -1;
	gravity = 0;
	frictionCoeff = 0.0;
	nu = 0.3;
	E = 0;
	kBond = 0;
	kWall = 10;
	kArea = 0;
	We0 = 0;
	Wmix0 = 0;
	Chi = 0;
	setElasticConstants();
	real eps[] = { 1.0 };
	real rCut[] = { pow(2.0, 1.0 / 6.0) };
	real uShift[] = { 1.0 };
	setPairForce(eps, rCut, uShift);

}

void Interactions::reset(){
	rCut.clear();
	uShift.clear();
	epsilon.clear();
}

void Interactions::setElasticConstants(){
	lambda = E*nu / ((1 - 2 * nu)*(1 + nu));
	mu = E / (2 * (1 + nu));
	//We0 = 0;
	//Wmix0 = 0;
	//Chi = 0;

}


void Interactions::setPairForce(real *epsA, real *rCutA, real *ushiftA){
	rCutMax = 0;
	int numTypes = sys->particles.numTypes;
	epsilon.init(numTypes);
	rCut.init(numTypes);
	uShift.init(numTypes);
	for (int i = 0; i < numTypes; i++){
		for (int j = 0; j<numTypes; j++){
			epsilon[i][j] = epsA[i*numTypes + j];
			rCut[i][j] = rCutA[i*numTypes + j];
			if (rCut[i][j] > rCutMax) rCutMax = rCut[i][j];
			uShift[i][j] = ushiftA[i*numTypes + j];
		}
	}
}

void Interactions::calcForces(){
	switch (sys->pbcType){
	case NOPBC:
		if (sys->walls.moveWalls) 
			sys->walls.move(sys->dt);
		calcWallForces();
		break;
	case XPBC:
		calcTopBottomWallForces();
	case XYPBC:
		break;
	default:
		break;
	}
	if (frictionForce) calcFrictionForces();
	if (gravityForce) calcGravityForces();
	if (bondForces) calcBondForces();
	if (feaForces) calcFeaForces();
	if (sys->useNN)
		calcPairForcesNN();
	else
		calcPairForces();
	//sys->particles.elements.checkAreas(sys->box);
	//printf("%f\n", sys->potEnergy);
}

void Interactions::calcGravityForces(){
	if (gravity > 0){
		#pragma omp parallel for
		for (int i = 0; i < sys->particles.n_particles; i++){
			sys->particles.force[i][1] -= sys->particles.mass[i] * gravity;
		}
	}
}

void Interactions::calcFrictionForces() {
	if (frictionCoeff > 0) {
//#pragma omp parallel for
		for (int i = 0; i < sys->particles.n_particles; i++) {
			sys->particles.force[i] -= sys->particles.vel[i] * frictionCoeff;
		}
	}
}

void Interactions::calcBondForces(){
	calcBondForcesCPU();
}

void Interactions::calcBondForcesCPU(){

		//#pragma omp parallel for
	for (int l = 0; l < sys->bonds.numBonds; l++){
		//for (int l = 0; l < 1; l++){
		int i1 = sys->bonds.bondList[2 * l];
		int i2 = sys->bonds.bondList[2 * l + 1];
		real r0 = sys->bonds.bondLength[l];
		int i = i1;
		int j = i2;
		VectorR dr;
		dr = sys->particles.pos[i] - sys->particles.pos[j];
		nearestImage(dr, sys->box, sys->pbcType);
		real r = dr.norm();
		real deltar = r - r0;
		//kBond = E / r0;
		//kBond = 10;
		real kdr = -kBond*deltar / r;
		VectorR f = kdr*dr;
		sys->particles.force[i] += f;
		sys->particles.force[j] -= f;
		real vir = kBond*deltar*deltar;
		real u = 0.5*vir;
		sys->potEnergy += u;
		//sys->virial += vir;
	}
}


void Interactions::calcWallForces(){
	sys->walls.forces[0].setZero();
	sys->walls.forces[1].setZero();
	real rc = pow(2.0, 1.0 / 6.0);

	for (int i = 0; i<sys->particles.n_particles; i++){
		real s0 = sys->particles.sigma[i];
		//index 1, top right and front wall
		for(int idim=0; idim<NDIM;idim++){
			if ((sys->walls.pos[1][idim] - sys->particles.pos[i][idim]) < s0*rc) { // collision detection
				real di = sys->walls.pos[1][idim] - sys->particles.pos[i][idim];
				real f;
				real u = ljwallforce(di, s0, f);
				sys->particles.force[i][idim] += -f;
				sys->walls.forces[1][idim] += f;
				sys->potEnergy += u;
			}
		}
		//index 0 bottom, left and back wall
		for (int idim = 0; idim < NDIM; idim++) {
			if ((sys->particles.pos[i][idim] - sys->walls.pos[0][idim]) < s0*rc) { // collision detection
				real di = sys->particles.pos[i][idim] - sys->walls.pos[0][idim];
				real f;
				real u = ljwallforce(di, s0, f);
				sys->particles.force[i][idim] += f;
				sys->walls.forces[0][idim] += -f;
				sys->potEnergy += u;
			}
		}
	}

	sys->wallPressure = 0; //FIXME
}

void Interactions::calcTopBottomWallForces(){
	sys->walls.forces[0].setZero();
	sys->walls.forces[1].setZero();
	real rc = pow(2.0, 1.0 / 6.0);
	int idim = 2;
	for (int i = 0; i<sys->particles.n_particles; i++){
		real s0 = sys->particles.sigma[i];
		//index 1, top right and front wall
		if ((sys->walls.pos[1][idim] - sys->particles.pos[i][idim]) < s0*rc) { // collision detection
			real di = sys->walls.pos[1][idim] - sys->particles.pos[i][idim];
			real f;
			real u = ljwallforce(di, s0, f);
			sys->particles.force[i][idim] += -f;
			sys->walls.forces[1][idim] += f;
			sys->potEnergy += u;
		}
		//index 0 bottom, left and back wall
		if ((sys->particles.pos[i][idim] - sys->walls.pos[0][idim]) < s0*rc) { // collision detection
			real di = sys->particles.pos[i][idim] - sys->walls.pos[0][idim];
			real f;
			real u = ljwallforce(di, s0, f);
			sys->particles.force[i][idim] += f;
			sys->walls.forces[0][idim] += -f;
			sys->potEnergy += u;
		}
	}
	sys->wallPressure /= 4;
}


//if (m1 != m2 || j2 < j1) {
//if ((m1 != m2 || j2 < j1) && (mol[j1].inChain == -1 || mol[j1].inChain != mol[j2].inChain || abs(j1 - j2) > 1))

/*
Condition
i       j
-1 and -1  interaction
g  and -1  interaction
-1 and  g  interaction
g  and  g  no interaction

*/

real ljwallforce(real zi, real sigma, real &fz){
	zi = sigma / zi;
	real zi3 = zi*zi*zi;
	real zi6 = zi3*zi3*zi3;
	real u = real(4.0)*zi6 * (zi6 - real(1.0)) + real(1.0);
	real fcVal = real(48.0) * zi6 * (zi6 - real(0.5)) / zi;
	fz = fabs(fcVal);
	return u;

}


real ljforce(real dr2, real dr, VectorR &drvec, real eps, real sigma, real uShift, VectorR &f, real &vir){
	real c = real(48.0) * eps / (sigma*sigma);
	real dr2i = sigma*sigma / dr2;
	real dr6i = dr2i*dr2i*dr2i;
	real fcVal = c*dr6i*(dr6i - real(0.5))*dr2i;
	f = fcVal*drvec;
	vir = fcVal*dr2;

	real u = real(4.0) * eps*dr6i*(dr6i - real(1.0)) + uShift;
	return u;
}

real pairForce(
	int i,
	int j,
	VectorR* pos,
	int* type,
	real *sigma,
	real *epsilon,
	real *rCut,
	real *uShift,
	int numTypes,
	VectorR box,
	PBCTYPE pbcType,
	VectorR *force,
	real &vir,
	VectorR &fij
	){

	real u = 0;
	VectorR f;
	VectorR dr;
	int typei = type[i];
	int typej = type[j];
	real eij = epsilon[typei*numTypes + typej];
	real rCutij = rCut[typei*numTypes + typej];
	real uShiftij = uShift[typei*numTypes + typej];
	real sigmaij = (sigma[i] + sigma[j]) / 2;

	dr = pos[i] - pos[j];
	nearestImage(dr, box, pbcType);
	real dr2 = dr.squaredNorm();
	real rCut2 = rCutij*rCutij*sigmaij*sigmaij;

	if (dr2 < rCut2){
		//real dr1 = sqrt(dr2);
		u = ljforce(dr2, 0.0, dr, eij, sigmaij, uShiftij, f, vir);
		force[i] += f;
		fij = f;
	}
	return u;
}


void Interactions::calcPairForces(){

		//#pragma omp parallel for schedule(static)
		sys->virial = 0;
		VectorR fij;
		for (int i = 0; i < sys->particles.n_particles; i++){
			//real temp = 0;
			for (int j = 0; j < sys->particles.n_particles; j++){
				int gi = sys->exclusion_groups.id_of_particle[i];
				int gj = sys->exclusion_groups.id_of_particle[j];
				if ((i != j) && (gi == -1 || gj == -1 || (gi != gj))){
				//if (i != j){
					real vir = 0.0;
					real u = pairForce(i, j, sys->particles.pos.data(), sys->particles.type.data(), sys->particles.sigma.data(),
						epsilon.data, rCut.data, uShift.data, sys->particles.numTypes, sys->box, sys->pbcType, sys->particles.force.data(),
						vir,fij);
					//#pragma omp atomic
					sys->virial += 0.5*vir;
					sys->potEnergy += 0.5*u;
					sys->pairEnergy += 0.5*u;
					if (sys->clusters.n_clusters > 0){
						int ci = sys->clusters.cluster_id_of_particle[i];
						int cj = sys->clusters.cluster_id_of_particle[j];
						VectorR rci, rcj, drcij;
						rci = sys->clusters.center_of_mass[ci];
						rcj = sys->clusters.center_of_mass[cj];
						//drcij.x = rci.x - rcj.x;
						//drcij.y = rci.y - rcj.y;
						//real dr2 = drcij.x*drcij.x + drcij.y*drcij.y;
						//sys->clusterVirial += fij.x*drcij.x + fij.y*drcij.y;
					}
					//temp += 0.5*u;
				}
			}
			//printf("%f\n", temp);
		}
}

void Interactions::calcPairForcesNN(){

		//#pragma omp parallel for schedule(static)
		for (int pi = 0; pi < sys->particles.n_particles; pi++){
			VectorR fij;
			int pj;
			real u;
			for (int j = 0; j < sys->neighborlist.neighborCount[pi]; j++){
				pj = sys->neighborlist.neighborList[pi*sys->neighborlist.maxNeighborsPerParticle + j];
				real vir = 0;
				u = pairForce(pi, pj, sys->particles.pos.data(), sys->particles.type.data(), sys->particles.sigma.data(),
						epsilon.data, rCut.data, uShift.data, sys->particles.numTypes, sys->box, sys->pbcType, sys->particles.force.data(), 
						vir, fij);

				//#pragma omp atomic
				sys->potEnergy += u*real(0.5);
				sys->pairEnergy += u*real(0.5);
				//#pragma omp atomic
				sys->virial += vir*real(0.5);
				if (sys->clusters.n_clusters > 0){
					int ci = sys->clusters.cluster_id_of_particle[pi];
					int cj = sys->clusters.cluster_id_of_particle[pj];
					VectorR rci, rcj, drcij;
					rci = sys->clusters.center_of_mass[ci];
					rcj = sys->clusters.center_of_mass[cj];
					//drcij.x = rci.x - rcj.x;
					//drcij.y = rci.y - rcj.y;
					//nearestImage(drcij, sys->box, sys->pbcType);
					//real dr2 = drcij.x*drcij.x + drcij.y*drcij.y;
					//sys->clusterVirial += 0.5*(fij.x*drcij.x + fij.y*drcij.y);
					//printf("%d %d %f %f %f\n",ci, cj, dr2, fij, fij*dr2);
				}

			}
		}
		//printf("%f\n", sys->virial);
}

///// Tetra Forces //////////////////////

void Interactions::calcFeaForces() {
		//calcFeaForcesCPU();
		sys->elements.offset = 0;
		//sys->elements.calcClusterProps();
		calcFeaForcesCPU1();

}

void Interactions::calcFeaForcesCPU1() {

	VectorR box = sys->box;
	MatrixDxD F, P, H, E, I, invF;
	I.setIdentity();
	int offset = sys->elements.offset;
	sys->elements.totCurrVol = 0;
	//printf("In calcFeaForcesCPU1\n");

	//#pragma omp parallel num_threads(8)
	//#pragma omp for
	for (int it = 0; it < sys->elements.tetras.size(); it++){
		VectorR rij[NDIM], f[NDIM + 1];
		MatrixDxD Ds;
		elementIndexes pind;

		// Get indexes of particle in cell
		for (int i = 0; i < NDIM+1; i++) 
			pind[i] = sys->elements.tetras[it][i];

		// Get vertor edges ri - r_dim i.e. (r0 -r2, r1 -r2 in 2D)
		for (int edge = 0; edge < NDIM; edge++) {
			rij[edge] = sys->particles.pos[pind[edge]] - sys->particles.pos[pind[NDIM]]; //edges
		}

		// Check that particles are inside box
		//Calculate matrix Ds page 28 equ 4.4
		for (int edge = 0; edge < NDIM; edge++) { // for each edge
			for (int i = 0; i < NDIM; i++) {
				if (rij[edge][i] >  sys->box[i] / 2) rij[edge][i] -= sys->box[i];
				if (rij[edge][i] < -sys->box[i] / 2) rij[edge][i] += sys->box[i];
			}
			Ds.col(edge) = rij[edge];
		}
		real vol = sys->elements.refVol[it];		
		real detDs = Ds.determinant();
		real newVol = detDs*one_over_ndim_factorial;
		if (std::isnan(newVol)) 
			std::cout << sys->elements.Bm[it] << "\n###########\n" << Ds << std::endl << std::endl;
		sys->elements.currVol[it] = newVol;
		sys->elements.totCurrVol += abs(newVol);
		if (vol*newVol < 0)printf("Error in triangle %d %f %f \n", it, newVol, vol);
		vol = fabs(vol);

		// Agorithm from fendefo.org page 30
		// calculated Dm from ref pos (Xi - X0)
		// Bm inverse of Dm Bm = Inverse(Dm)
		// calculated Ds from world positions Ds = (xi-x0)
		// Deformation gradient F = Ds*Bm
		F = Ds * sys->elements.Bm[it];

		// St Venant-Kirchooff 
		// Green lagrange strain tensor E = 0.5(Transpose(F)*F - 1)
		// P(F) = F(2*mu*E + lambda*Tr(E)*I)
		real trE, ltrE, detF, pe;
		if (elasticForceType == 0){
			E = 0.5*(F.transpose()*F - I);
			trE = E.trace();
			ltrE = lambda*trE;
			P = ltrE * F + 2 * mu * F * E;
			pe = mu*E.squaredNorm() + 0.5*lambda*trE*trE;
		}
		//NeoHookean
		//P(F) = mu(F-Transpose(inv(F)) + lambda*log(J)*Transpose(inv(F))
		else if(elasticForceType == 1){
			detF = F.determinant();
			real I1 = F.squaredNorm();
			real logJ = log(detF);
			real llogJ = lambda*logJ;
			invF = F.inverse();
			P = mu*F + (llogJ - mu) * invF.transpose();
			//P.noalias() = mu*F + (llogJ - mu) * F.inverse().transpose();
			pe = 0.5*mu*(I1 - 2) - mu*logJ + 0.5*lambda*logJ*logJ;
		}
		//Flory-Huggins
		else if (elasticForceType == 2) {
			real We, Wm;
			real J = F.determinant(); //det(F)
			real logJ = log(J);
			real logJ_1overJ = log((J - 1) / J);
			//We0 = mu;
			real lambda1 = lambda;
			lambda1 = 0;
			real llogJ = lambda1*logJ;
			MatrixDxD C = F.transpose()*F;
			real trC = C.trace();
			real trC2 = (C*C).trace();
			real I2 = 0.5*(trC*trC - trC2);
			sys->elements.I2[it] = I2;
			MatrixDxD invF = F.inverse();
			
			//Wmix0 = 0;
			real dWmix_dJ = 0;
			if (J - 1 < 0.0001) {
				//dWmix_dJ = Wmix0*(Chi / J + 1.0 / J);
				dWmix_dJ = 0;
				Wm = 0;
			}
			else {
				dWmix_dJ = Wmix0*(Chi / J + 1.0 / J + logJ_1overJ);
				Wm = Wmix0*(J - 1)*(logJ_1overJ + Chi / J);
			}
			P = We0*F + (dWmix_dJ*J - We0 + llogJ)*invF.transpose();
			
			real I1 = F.squaredNorm();
			sys->elements.I1[it] = I1;
			//pe = 0.5*mu*(I1 - 2) - mu*logJ + 0.5*lambda*logJ*logJ;
			We = 0.5*We0*(I1 - NDIM) - mu*logJ + 0.5*lambda1*logJ*logJ ;
			pe = We + Wm;

		}
		real vol0 = vol;
		//Force calculation
		H = -vol0 * P * sys->elements.Bm[it].transpose();

		f[NDIM].setZero();
		for (int ivert = 0; ivert < NDIM; ivert++) {
			f[ivert] = H.col(ivert);
			f[NDIM] -= f[ivert];
		}
		
		for (int i = 0; i < NDIM + 1; i++){
			int j = pind[i];
			//#pragma omp atomic
			sys->particles.force[j + offset] += f[i];
			//#pragma omp atomic
			sys->feaVirial += f[i].dot(sys->elements.unfoldedPos[j + offset]);

		}
		//printf("E: %f\n", pe);
		sys->potEnergy += pe*vol0;
		sys->feaEnergy += pe*vol0;

	}
	//printf("currVol: %f\n", sys->elements.totCurrVol);
}



void Interactions::calcFeaForcesCPU2() {

	VectorR box = sys->box;
	//printf("In calcFeaForcesCPU2\n");
	int offset = sys->elements.offset;
	sys->elements.totRefVol = 0;
	sys->elements.totCurrVol = 0;

	//#pragma omp parallel num_threads(8)
	//#pragma omp for
	for (int it = 0; it < sys->elements.tetras.size(); it++) {
		VectorR u[3], r[3], r0[3], f[3];
		int pind[3];
		for (int i = 0; i < 3; i++) {
			int pi = sys->elements.tetras[it][i];
			pind[i] = pi;
			r[i] = sys->particles.pos[pi];
			r0[i] = sys->elements.refPos[pi];
		}

		real dx01 = r[0].x() - r[1].x();
		if (dx01 >  sys->box.x() / 2.0) r[1].x() += sys->box.x();
		if (dx01 < -sys->box.x() / 2.0) r[1].x() -= sys->box.x();

		real dx02 = r[0].x() - r[2].x();
		if (dx02 > sys->box.x() / 2.0) r[2].x() += sys->box.x();
		if (dx02 < -sys->box.x() / 2.0) r[2].x() -= sys->box.x();

		real dy01 = r[0].y() - r[1].y();
		if (dy01 > sys->box.y() / 2.0) r[1].y() += sys->box.y();
		if (dy01 < -sys->box.y() / 2.0) r[1].y() -= sys->box.y();

		real dy02 = r[0].y() - r[2].y();
		if (dy02 > sys->box.y() / 2.0) r[2].y() += sys->box.y();
		if (dy02 < -sys->box.y() / 2.0) r[2].y() -= sys->box.y();

		//real vol = 0.5*((r[1].x*r[2].y - r[2].x*r[1].y) - (r[2].x*r[0].y-r[0].x*r[2].y) - (r[1].x*r[2].y-r[2].x*r[1].y));
		//real vol = 0.5*((r0[1].x*r0[2].y - r0[2].x*r0[1].y) - (r0[2].x*r0[0].y - r0[0].x*r0[2].y) - (r0[1].x*r0[2].y - r0[2].x*r0[1].y));
		real vol = 0.5*((r0[2].y() - r0[0].y())*(r0[1].x() - r0[0].x()) - (r0[1].y() - r0[0].y())*(r0[2].x() - r0[0].x()));
		real newVol = 0.5*((r[2].y() - r[0].y())*(r[1].x() - r[0].x()) - (r[1].y() - r[0].y())*(r[2].x() - r[0].x()));
		sys->elements.refVol[it] = vol;
		sys->elements.currVol[it] = newVol;
		sys->elements.totRefVol += vol;
		sys->elements.totCurrVol += newVol;
		if (vol*newVol < 0)printf("Error in triangle %d %f %f \n", it, newVol, vol);
		vol = fabs(vol);
		//real vol = e->tetVol[itet];
		//xm = sys->particles.elements.xm[it];

		real Dm[2][2];
		real Bm[2][2];
		real Ds[2][2];
		real F[2][2];
		real P[2][2];
		real H[2][2];
		real E[2][2];
		//Agorithm from fendefo.org page 30

		// calculate Dm from ref pos (Xi - X0)
		Dm[0][0] = r0[0].x() - r0[2].x();
		Dm[1][0] = r0[0].y() - r0[2].y();
		Dm[0][1] = r0[1].x() - r0[2].x();
		Dm[1][1] = r0[1].y() - r0[2].y();

		//det and volume of undeformed triangle;
		real detDm = Dm[0][0] * Dm[1][1] - Dm[0][1] * Dm[1][0];
		real vol0 = detDm / 2;

		// Bm inverse of Dm Bm = Inverse(Dm)
		Bm[0][0] = Dm[1][1] / detDm;
		Bm[1][0] = -Dm[1][0] / detDm;
		Bm[0][1] = -Dm[0][1] / detDm;
		Bm[1][1] = Dm[0][0] / detDm;

		//calculate Ds from world positions Ds = (xi-x0)
		Ds[0][0] = r[0].x() - r[2].x();
		Ds[1][0] = r[0].y() - r[2].y();
		Ds[0][1] = r[1].x() - r[2].x();
		Ds[1][1] = r[1].y() - r[2].y();

		//Deformation gradient F = Ds*Bm
		F[0][0] = Ds[0][0] * Bm[0][0] + Ds[0][1] * Bm[1][0];
		F[1][0] = Ds[1][0] * Bm[0][0] + Ds[1][1] * Bm[1][0];
		F[0][1] = Ds[0][0] * Bm[0][1] + Ds[0][1] * Bm[1][1];
		F[1][1] = Ds[1][0] * Bm[0][1] + Ds[1][1] * Bm[1][1];

		//Green lagrange strain tensor E = 0.5(Transpose(F)*F - 1)
		E[0][0] = 0.5*(F[0][0] * F[0][0] + F[1][0] * F[1][0] - 1);
		E[1][0] = 0.5*(F[0][0] * F[0][1] + F[1][1] * F[1][0]);
		E[0][1] = 0.5*(F[0][0] * F[0][1] + F[1][0] * F[1][1]);
		E[1][1] = 0.5*(F[0][1] * F[0][1] + F[1][1] * F[1][1] - 1);

		//printf("\nE: %f %f %f %f\n", E[0][0], E[1][0], E[0][1], E[1][1]);
		//printf("\F: %f %f %f %f\n", F[0][0], F[1][0], F[0][1], F[1][1]);

		// Piola tensor
		//St Venant-Kirchooff model
		//P(F) = F(2*mu*E + lambda*Tr(E)*I)

		real trE, ltrE, detF, pe;
		if (elasticForceType == 0) {
			trE = (E[0][0] + E[1][1]);
			ltrE = lambda*trE;
			P[0][0] = ltrE * F[0][0] + 2 * mu* (F[0][0] * E[0][0] + F[0][1] * E[1][0]);
			P[1][0] = ltrE * F[1][0] + 2 * mu* (F[1][0] * E[0][0] + F[1][1] * E[1][0]);
			P[0][1] = ltrE * F[0][1] + 2 * mu* (F[0][0] * E[0][1] + F[0][1] * E[1][1]);
			P[1][1] = ltrE * F[1][1] + 2 * mu* (F[1][0] * E[0][1] + F[1][1] * E[1][1]);
			pe = mu*(E[0][0] * E[0][0] + E[1][0] * E[1][0] + E[0][1] * E[0][1] + E[1][1] * E[1][1]) + 0.5*lambda*trE*trE;
		}
		else if (elasticForceType == 1) {
			//NeoHookean
			//
			detF = F[0][0] * F[1][1] - F[0][1] * F[1][0];
			real logJ = log(detF);
			real llogJ = lambda*logJ;
			real invF[2][2];
			invF[0][0] = F[1][1] / detF;
			invF[1][0] = -F[1][0] / detF;
			invF[0][1] = -F[0][1] / detF;
			invF[1][1] = F[0][0] / detF;
			P[0][0] = mu*(F[0][0] - invF[0][0]) + llogJ*invF[0][0];
			P[1][0] = mu*(F[1][0] - invF[0][1]) + llogJ*invF[0][1];
			P[0][1] = mu*(F[0][1] - invF[1][0]) + llogJ*invF[1][0];
			P[1][1] = mu*(F[1][1] - invF[1][1]) + llogJ*invF[1][1];
			real I1 = F[0][0] * F[0][0] + F[1][0] * F[1][0] + F[0][1] * F[0][1] + F[1][1] * F[1][1];
			pe = 0.5*mu*(I1 - 2) - mu*logJ + 0.5*lambda*logJ*logJ;
		}
		else if (elasticForceType == 2) {
			real We, Wm;
			real J = F[0][0] * F[1][1] - F[0][1] * F[1][0]; //det(F)
			real logJ = log(J);
			real logJ_1overJ = log((J - 1) / J);
			//We0 = mu;
			real lambda1 = lambda;
			lambda1 = 0;
			real llogJ = lambda1*logJ;

			real invF[2][2];
			invF[0][0] = F[1][1] / J;
			invF[1][0] = -F[1][0] / J;
			invF[0][1] = -F[0][1] / J;
			invF[1][1] = F[0][0] / J;

			//Wmix0 = 0;
			real dWmix_dJ = 0;
			if (J - 1 < 0.0001) {
				//dWmix_dJ = Wmix0*(Chi / J + 1.0 / J);
				dWmix_dJ = 0;
				Wm = 0;
			}
			else {
				dWmix_dJ = Wmix0*(Chi / J + 1.0 / J + logJ_1overJ);
				Wm = Wmix0*(J - 1)*(logJ_1overJ + Chi / J);
			}

			P[0][0] = We0*(F[0][0] - invF[0][0]) + llogJ*invF[0][0] + dWmix_dJ*J*invF[0][0];
			P[1][0] = We0*(F[1][0] - invF[0][1]) + llogJ*invF[0][1] + dWmix_dJ*J*invF[0][1];
			P[0][1] = We0*(F[0][1] - invF[1][0]) + llogJ*invF[1][0] + dWmix_dJ*J*invF[1][0];
			P[1][1] = We0*(F[1][1] - invF[1][1]) + llogJ*invF[1][1] + dWmix_dJ*J*invF[1][1];

			real I1 = F[0][0] * F[0][0] + F[1][0] * F[1][0] + F[0][1] * F[0][1] + F[1][1] * F[1][1];
			//pe = 0.5*mu*(I1 - 2) - mu*logJ + 0.5*lambda*logJ*logJ;
			We = 0.5*We0*(I1 - 2) - mu*logJ + 0.5*lambda1*logJ*logJ;
			pe = We + Wm;

		}

		//Force calculation
		H[0][0] = -vol0*(P[0][0] * Bm[0][0] + P[0][1] * Bm[0][1]);
		H[1][0] = -vol0*(P[1][0] * Bm[0][0] + P[1][1] * Bm[0][1]);
		H[0][1] = -vol0*(P[0][0] * Bm[1][0] + P[0][1] * Bm[1][1]);
		H[1][1] = -vol0*(P[1][0] * Bm[1][0] + P[1][1] * Bm[1][1]);

		//        calculate the x,y forces on the 3 nodes
		//VectorR f[3];
		f[0].x() = H[0][0];
		f[0].y() = H[1][0];
		f[1].x() = H[0][1];
		f[1].y() = H[1][1];
		f[2].x() = -f[0].x() - f[1].x();
		f[2].y() = -f[0].y() - f[1].y();

		for (int i = 0; i < 3; i++) {
			int j = pind[i];
			//#pragma omp atomic
			sys->particles.force[j + offset].x() += f[i].x();
			//#pragma omp atomic
			sys->particles.force[j + offset].y() += f[i].y();
			//#pragma omp atomic
			sys->feaVirial += f[i].x()*sys->elements.unfoldedPos[j + offset].x() +
				f[i].y()*sys->elements.unfoldedPos[j + offset].y();

		}
		sys->potEnergy += pe*vol0;
		sys->feaEnergy += pe*vol0;
	}
	//printf("Vols: %f %f\n", sys->elements.totRefVol, sys->elements.totCurrVol);
}


