#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "md3dsystem.h"
#include "properties.h"


void Properties::evaluate()
{
	sys->evalProps();
	uKinetic.val = 0.5*sys->kinEneSum / sys->particles.n_particles;
	uPotential.val = sys->potEnergy /sys->particles.n_particles;
	uFea.val = sys->feaEnergy / sys->particles.n_particles;
	uPair.val = sys->pairEnergy / sys->particles.n_particles;
	uTotal.val = uKinetic.val + uPotential.val;
	temperature.val = uKinetic.val;

	clusterUKinetic.val = 0.5*sys->clusterKinEneSum / sys->clusters.n_clusters;

	real V = sys->box.prod();
	real dens = sys->particles.n_particles / V;
	//virial.val = sys->virial/(2.0*sys->N);
	virial.val = sys->virial;
	//clusterVirial.val = sys->clusterVirial / (2 * sys->NC);
	clusterVirial.val = sys->clusterVirial;
	//feaVirial.val = sys->feaVirial / (2 * sys->N);
	feaVirial.val = sys->feaVirial;
	//pressure.val = sys->density*(sys->kinEneSum + sys->virial + sys->feaVirial) / (2 * sys->N);
	pressure.val = 0.5*sys->kinEneSum/V + 0.5*(sys->virial + sys->feaVirial)/V;
	molPressure.val = 0.5*sys->clusterKinEneSum/V - 0.5* sys->clusterVirial / V;
}


void Properties::accumulate(int icode)
{
	//Prop uPotential, uKinetic, uBonded[10], virial, pressure;
	if (icode == 0) {
		count = 0;
		uPotential.zero();
		uKinetic.zero();
		clusterUKinetic.zero();
		uTotal.zero();
		pressure.zero();
		molPressure.zero();
		virial.zero();
		clusterVirial.zero();
		feaVirial.zero();
		uBonded.zero();
		uFea.zero();
		uPair.zero();
		temperature.zero();
	}
	else if (icode == 1) {
		count++;
		uPotential.accum();
		uKinetic.accum();
		clusterUKinetic.accum();
		uTotal.accum();
		pressure.accum();
		molPressure.accum();
		virial.accum();
		clusterVirial.accum();
		feaVirial.accum();
		uBonded.accum();
		uFea.accum();
		uPair.accum();
		temperature.accum();
	}
	else if (icode == 2) {
		uPotential.avg(count);
		uKinetic.avg(count);
		clusterUKinetic.avg(count);
		uTotal.avg(count);
		pressure.avg(count);
		molPressure.avg(count);
		virial.avg(count);
		clusterVirial.avg(count);
		feaVirial.avg(count);
		uBonded.avg(count);
		uFea.avg(count);
		uPair.avg(count);
		temperature.avg(count);
	}
}

void Properties::compute(){
	evaluate();
	accumulate(1);
	if (sys->steps > sys->props.step_equi){
		evalClusterRdf();
		if (sys->steps % sys->props.step_avg == 0) {
			sys->props.accumulate(2);
			sys->props.accumulate(0);
		}
	}
	/*
	printf("%d %f %f %f %f %f %f\n%f %f %f\n",
		sys->steps,
		sys->props.uKinetic.val,
		sys->props.uPotential.val,
		sys->props.uTotal.val,
		sys->props.virial.val,
		sys->props.feaVirial.val,
		sys->props.pressure.val,
		sys->props.clusterUKinetic.val,
		sys->props.clusterVirial.val,
		sys->props.molPressure.val
	);
	*/
}

void printAvgProps(){
	//printf()
}


void Properties::evalClusterRdf()
{
	/*
	if (sys->steps%skipStepsRdf != 0) return;
	VectorR dr;
	real deltaR, normFac, rr;
	int j1, j2, n;
	VectorR box = sys->box;
	//box.x = box.x / sys->scale;
	//box.y = box.y / sys->scale;
	rangeRdf = box.x / 2;

	if (countRdf == 0) {
		for (n = 0; n < sizeHistRdf; n++) histRdf[n] = 0.;
	}
	deltaR = rangeRdf / sizeHistRdf;
	for (j1 = 0; j1 < sys->n_clusters - 1; j1++) {
		for (j2 = j1 + 1; j2 < sys->n_clusters; j2++) {
			//VSub(dr, mol[j1].r, mol[j2].r);
			dr.x = sys->clusters.data[j1].cmPos.x - sys->clusters.data[j2].cmPos.x;
			dr.y = sys->clusters.data[j1].cmPos.y - sys->clusters.data[j2].cmPos.y;

			//VWrapAll(dr);
			//dr.x = dr.x / sys->scale;
			//dr.y = dr.y / sys->scale;
			nearestImage(dr, box, sys->pbcType);
			//rr = VLenSq(dr);
			rr = dr.x*dr.x + dr.y*dr.y;
			if (rr < rangeRdf*rangeRdf) {
				n = sqrt(rr) / deltaR;
				++histRdf[n];
				//printf("%d %d\n", n, histRdf[n]);
			}
		}
	}
	++countRdf;
	if (countRdf == avgCountRdf) {
		char name[255];
		sprintf_s(name, "%s%.9d.dat", "rdf", sys->steps);
		//FILE *fp = fopen_s(name, "w");
		FILE *fp;
		errno_t err = fopen_s(&fp, name, "w");
		normFac = sys->box.x*sys->box.y / (PI * deltaR*deltaR * sys->n_clusters*sys->n_clusters * countRdf);
		for (n = 0; n < sizeHistRdf; n++){
			histRdf[n] *= normFac / (n - 0.5);
			//printf("%f %f %f %d %d\n", normFac, rangeRdf, deltaR, n, histRdf[n]);
		}
		real rb;
		int n;
		//fprintf(fp, "rdf\n");
		for (n = 0; n < sizeHistRdf; n++) {
			rb = (n + 0.5) * rangeRdf / sizeHistRdf;
			fprintf(fp, "%8.4f %8.4f\n", rb, histRdf[n]);
		}
		fflush(fp);
		countRdf = 0;
	}
	*/
}


