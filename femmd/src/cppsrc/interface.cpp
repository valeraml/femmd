

#ifdef _WIN32
#include <direct.h>
#elif defined __linux__
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "md3dsystem.h"
#include <Python.h>
#include <numpy/arrayobject.h>


#include <cmath>   

MDSystem *sys;

int c_main();

PyObject* md_system_c_main(PyObject* pSelf, PyObject* pArgs) {
	c_main();
	Py_INCREF(Py_None);
	return Py_None;
}

void md_system_set_defaults() {
	sys->useNN = false;
	sys->dt = 0.001;
	sys->props.step_avg = 1000;
	sys->props.step_equi = 0;
	sys->particles.numTypes = 1;
	sys->pbcType = XYPBC;
	sys->interactions.gravity = 0;
	sys->interactions.gravityForce = false;
	sys->interactions.E = 5;
	sys->interactions.kBond = 10;
	sys->interactions.kArea = .0001;
	sys->interactions.bondForces = false;
	sys->interactions.feaForces = false;
	sys->interactions.areaForces = false;
	sys->interactions.setElasticConstants();

	double v0 = 0.5;
	//sys->particles.setParticlesVel(v0);
	//if (sys->useNN)
	//	sys->neighborlist.init();

	sys->start = true;
}

PyObject* md_system_init(PyObject* pSelf, PyObject* pArgs)
{
	PyArrayObject *box;
	if (!PyArg_ParseTuple(pArgs, "O", &box)) return NULL;
	real *box_data = (double*)PyArray_DATA(box);
	VectorR L;
	for (int i = 0; i < NDIM; i++)
		L[i] = box_data[i];
	MDSystem *system = new MDSystem(0);
	sys = system;
	sys->init(L);
	md_system_set_defaults();
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_reset(PyObject* pSelf, PyObject* pArgs)
{
	printf("system reset\n");
	sys->reset();
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_dt(PyObject *pSelf, PyObject *pArgs) {
	double dt;
	if (!PyArg_ParseTuple(pArgs, "d", &dt)) return NULL;
	sys->dt = dt;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_add_particle(PyObject* pSelf, PyObject* pArgs)
{
	PyArrayObject *r_, *v_;
	real m, sig;
	int typ;
	if (!PyArg_ParseTuple(pArgs, "OOddi", &r_, &v_, &m, &sig, &typ)) return NULL;
	real *r_data = (double*)PyArray_DATA(r_);
	real *v_data = (double*)PyArray_DATA(v_);
	VectorR r, v;
	for (int i = 0; i < NDIM; i++) {
		r[i] = r_data[i];
		v[i] = v_data[i];
	}
	sys->particles.addParticle(r, v, m, sig, typ);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_add_cluster(PyObject* pSelf, PyObject* pArgs)
{
	//TODO check ints and types
	real x, y, mass;
	int typ, exc;
	PyArrayObject *vert;
	PyArrayObject *cells;
	PyArrayObject *r0;
	if (!PyArg_ParseTuple(pArgs, "OOOidi", &vert, &cells, &r0, &typ, &mass, &exc)) return NULL;
	int nv = PyArray_SIZE(vert)/NDIM;
	int nc = PyArray_SIZE(cells)/(NDIM+1);
	real *vert_data = (real*)PyArray_DATA(vert);
	int *cells_data = (int*)PyArray_DATA(cells);
	real *r0_data = (real*)PyArray_DATA(r0);
	VectorR r;
	for (int i = 0; i < NDIM; i++) 
		r[i] = r0_data[i];
	//mass = (real)nv; //FIXME really we dont need the mass
	sys->clusters.addCluster(vert_data, nv, cells_data, nc, r, typ, mass, exc);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_set_velocities(PyObject* pSelf, PyObject* pArgs)
{
	real v0;
	if (!PyArg_ParseTuple(pArgs, "d", &v0)) return NULL;
	sys->particles.setParticlesVel(v0);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_zero_center_of_mass_vel(PyObject* pSelf, PyObject* pArgs)
{
	sys->particles.zeroCM();
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_set_potential(PyObject* pSelf, PyObject* pArgs)
{
	PyArrayObject *eps1;
	PyArrayObject *rcut;
	PyArrayObject *shift;
	if (!PyArg_ParseTuple(pArgs, "OOO", &eps1, &rcut, &shift)) return NULL;
	int size = PyArray_SIZE(eps1);
	if (size != (sys->particles.numTypes*sys->particles.numTypes)) {
		printf("wrong dimension from potential in set_potential\n");
	}
	real *eps_data = (double*) PyArray_DATA(eps1);
	real *rcut_data = (double*)PyArray_DATA(rcut);
	real *shift_data = (double*)PyArray_DATA(shift);
	for (int i = 0; i < size; i++) {
		printf("inter param %d %f %f %f\n", i, eps_data[i], rcut_data[i], shift_data[i]);
	}
	sys->interactions.setPairForce(eps_data, rcut_data, shift_data);

	Py_RETURN_NONE;
}

PyObject *md_system_set_gravity(PyObject *pSelf, PyObject *pArgs) {
	double g;
	if (!PyArg_ParseTuple(pArgs, "d", &g)) return NULL;
	if (g > 0) {
		sys->interactions.gravity = g;
		sys->interactions.gravityForce = true;
	}
	else
		sys->interactions.gravityForce = false;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_friction(PyObject *pSelf, PyObject *pArgs) {
	double mu;
	if (!PyArg_ParseTuple(pArgs, "d", &mu)) return NULL;
	if (mu > 0) {
		sys->interactions.frictionCoeff = mu;
		sys->interactions.frictionForce = true;
	}
	else
		sys->interactions.gravityForce = false;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_bond_spring_constant(PyObject *pSelf, PyObject *pArgs) {
	double k;
	if (!PyArg_ParseTuple(pArgs, "d", &k)) return NULL;
	if (k > 0) {
		sys->interactions.kBond = k;
		sys->interactions.bondForces = true;
	}
	else
		sys->interactions.bondForces = false;
	Py_RETURN_NONE;
}

PyObject *md_system_set_area_force_constant(PyObject *pSelf, PyObject *pArgs) {
	double k;
	if (!PyArg_ParseTuple(pArgs, "d", &k)) return NULL;
	if (k > 0) {
		sys->interactions.kArea = k;
		sys->interactions.areaForces = true;
	}
	else
		sys->interactions.areaForces = false;
	Py_RETURN_NONE;
}

PyObject *md_system_set_youngs_modulus(PyObject *pSelf, PyObject *pArgs) {
	double E;
	if (!PyArg_ParseTuple(pArgs, "d", &E)) return NULL;
	if (E > 0) {
		sys->interactions.E = E;
		sys->interactions.feaForces = true;
		sys->interactions.elasticForceType = 1;
	}
	else
		sys->interactions.feaForces = false;
	sys->interactions.setElasticConstants();
	printf("Elastic Constants: E=%f nu=%f mu=%f lambda=%f\n",
		sys->interactions.E, sys->interactions.nu, sys->interactions.mu, sys->interactions.lambda);
	Py_RETURN_NONE;
}

PyObject *md_system_set_swelling_energy_constants(PyObject *pSelf, PyObject *pArgs) {
	double we0, wmix0, chi;
	if (!PyArg_ParseTuple(pArgs, "ddd", &we0, &wmix0, &chi)) return NULL;
	if (we0 > 0) {
		sys->interactions.We0 = we0;
		sys->interactions.Wmix0 = wmix0;
		sys->interactions.Chi = chi;
		sys->interactions.feaForces = true;
		sys->interactions.elasticForceType = 2;
	}
	else
		sys->interactions.feaForces = false;
	sys->interactions.setElasticConstants();
	printf("Swelling Constants: We0=%f Wmix0=%f Chi=%f\n",
		sys->interactions.We0, sys->interactions.Wmix0, sys->interactions.Chi);
	Py_RETURN_NONE;
}

PyObject *md_system_set_equilibrium_steps(PyObject *pSelf, PyObject *pArgs) {
	int equi_steps;
	if (!PyArg_ParseTuple(pArgs, "i", &equi_steps)) return NULL;
	sys->props.step_equi = equi_steps;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_avg_steps(PyObject *pSelf, PyObject *pArgs) {
	int avg_steps;
	if (!PyArg_ParseTuple(pArgs, "i", &avg_steps)) return NULL;
	sys->props.step_avg = avg_steps;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_scale_factor(PyObject *pSelf, PyObject *pArgs) {
	double s;
	if (!PyArg_ParseTuple(pArgs, "d", &s)) return NULL;
	sys->scale = s;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_box(PyObject *pSelf, PyObject *pArgs) {

	PyArrayObject *lbox;
	if (!PyArg_ParseTuple(pArgs, "O", &lbox)) return NULL;
	real *lbox_data = (double*)PyArray_DATA(lbox);
	VectorR box;
	for (int i = 0; i < NDIM; i++) {
		box[i] = lbox_data[i];
	}
	sys->setBox(box);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_get_box(PyObject *pSelf, PyObject *pArgs) {
	VectorR lbox;
	lbox = sys->box;
	if(NDIM == 2)
		return Py_BuildValue("(dd)", lbox[0], lbox[1]);
	else
		return Py_BuildValue("(ddd)", lbox[0], lbox[1], lbox[2]);
}

PyObject *md_system_get_walls_pos(PyObject *pSelf, PyObject *pArgs) {
	double left, right, top, bottom;
	hvector<VectorR> p = sys->walls.pos;
	if(NDIM==2)
		return Py_BuildValue("[[dd][dd]]", p[0][0], p[0][1], p[1][0], p[1][1]);
	else
		return Py_BuildValue("((ddd)(ddd))", p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}

PyObject *md_system_get_walls_forces(PyObject *pSelf, PyObject *pArgs) {
	double left, right, top, bottom;
	hvector<VectorR> f = sys->walls.forces;
	if (NDIM == 2)
		return Py_BuildValue("((dd)(dd))", f[0][0], f[0][1], f[1][0], f[1][1]);
	else
		return Py_BuildValue("((ddd)(ddd))", f[0][0], f[0][1], f[0][2], f[1][0], f[1][1], f[1][2]);
}

PyObject *md_system_set_walls_pos(PyObject *pSelf, PyObject *pArgs) {
	PyArrayObject *p0;
	PyArrayObject *p1;
	if (!PyArg_ParseTuple(pArgs, "OO", &p0, &p1)) return NULL;
	real *p0_data = (double*)PyArray_DATA(p0);
	real *p1_data = (double*)PyArray_DATA(p1);
	for (int i = 0; i < NDIM; i++) {
		sys->walls.pos[0][i] = p0_data[i];
		sys->walls.pos[1][i] = p1_data[i];
	}
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_walls_moving_rate(PyObject *pSelf, PyObject *pArgs) {
	PyArrayObject *r0;
	PyArrayObject *r1;
	if (!PyArg_ParseTuple(pArgs, "OO", &r0, &r1)) return NULL;
	real *r0_data = (double*)PyArray_DATA(r0);
	real *r1_data = (double*)PyArray_DATA(r1);
	for (int i = 0; i < NDIM; i++) {
		sys->walls.motion_rate[0][i] = r0_data[i];
		sys->walls.motion_rate[1][i] = r1_data[i];
	}
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_moving_walls(PyObject *pSelf, PyObject *pArgs) {
	int flag;
	if (!PyArg_ParseTuple(pArgs, "i", &flag)) return NULL;
	sys->walls.moveWalls = flag;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_walls_shear(PyObject *pSelf, PyObject *pArgs) {
	int flag;
	if (!PyArg_ParseTuple(pArgs, "i", &flag)) return NULL;
	sys->walls.shearWalls = flag;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *md_system_set_boundary_conditions(PyObject *pSelf, PyObject *pArgs) {
	int pbc;
	if (!PyArg_ParseTuple(pArgs, "i", &pbc)) return NULL;
	if (pbc == 0)
		sys->pbcType = NOPBC;
	else if (pbc == 1)
		sys->pbcType = XPBC;
	else if (pbc == 2)
		sys->pbcType = XYPBC;
	else
		sys->pbcType = XYZPBC;
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_get_positions(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.pos.data());
	if (p == NULL) printf("error creating array\n");
	return p;
}

PyObject* md_system_get_velocities(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.vel.data());
	if (p == NULL) printf("error creating array\n");
	return p;

}

PyObject* md_system_get_accelerations(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.acc.data());
	if (p == NULL) printf("error creating array\n");
	return p;

}

PyObject* md_system_get_forces(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (p == NULL) printf("error creating array\n");
	return p;

}

PyObject* md_system_get_unfolded_positions(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->elements.unfoldedPos.data());
	if (p == NULL) printf("error creating array\n");
	return p;
}

PyObject* md_system_unfold_positions(PyObject* pSelf, PyObject* pArgs)
{
	sys->elements.unfoldPos(sys->box);
	Py_RETURN_NONE;
}

PyObject* md_system_get_tetras(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->elements.numTetras;
	dims[1] = NDIM + 1;
	PyObject *tet = NULL;
	tet = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT, sys->elements.tetras.data());
	if (tet == NULL) printf("error creating array\n");
	return tet;
}

PyObject* md_system_get_bonds(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->bonds.bondListLen/2;
	dims[1] = 2;
	PyObject *bonds = NULL;
	bonds = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT, sys->bonds.bondList.data());
	if (bonds == NULL) printf("error creating array\n");
	return bonds;
}

PyObject* md_system_get_particle_sizes(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	//dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *p = NULL;
	p = (PyObject *)PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sys->particles.sigma.data());
	if (p == NULL) printf("error creating array\n");
	return p;

}

PyObject* md_system_get_particle_data(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->particles.n_particles;
	dims[1] = sizeof(VectorR) / sizeof(real);
	PyObject *ppos = NULL, *pupos = NULL, *pvel = NULL, *pacc = NULL, *pforce = NULL;
	ppos = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (ppos == NULL) printf("error creating array\n");
	pupos = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (pupos == NULL) printf("error creating array\n");
	pvel = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (pvel == NULL) printf("error creating array\n");
	pacc = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (pacc == NULL) printf("error creating array\n");
	pforce = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, sys->particles.force.data());
	if (pforce == NULL) printf("error creating array\n");

	dims[0] = sys->elements.numTetras;
	dims[1] = 3;
	PyObject *tet = NULL;
	tet = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT, sys->elements.tetras.data());
	if (tet == NULL) printf("error creating array\n");

	dims[0] = sys->bonds.bondListLen / 2;
	dims[1] = 2;
	PyObject *bonds = NULL;
	bonds = (PyObject *)PyArray_SimpleNewFromData(2, dims, NPY_INT, sys->bonds.bondList.data());
	if (bonds == NULL) printf("error creating array\n");

	return Py_BuildValue("(OOOOOOO)", ppos, pupos, pvel, pacc, pforce, tet, bonds);

}

PyObject* md_system_get_elements_volumes(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->elements.numTetras;
	dims[1] = 1;
	PyObject *p1 = NULL;
	p1 = (PyObject *)PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sys->elements.refVol.data());
	if (p1 == NULL) printf("error creating array\n");
	PyObject *p2 = NULL;
	p2 = (PyObject *)PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sys->elements.currVol.data());
	if (p2 == NULL) printf("error creating array\n");
	return Py_BuildValue("(OO)",p1,p2);

}

PyObject* md_system_get_elements_invariants(PyObject* pSelf, PyObject* pArgs)
{
	npy_intp dims[2];
	dims[0] = sys->elements.numTetras;
	dims[1] = 1;
	PyObject *p1 = NULL;
	p1 = (PyObject *)PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sys->elements.I1.data());
	if (p1 == NULL) printf("error creating array\n");
	PyObject *p2 = NULL;
	p2 = (PyObject *)PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sys->elements.I2.data());
	if (p2 == NULL) printf("error creating array\n");
	return Py_BuildValue("(OO)", p1, p2);

}

PyObject* md_system_init_neighbor_list(PyObject* pSelf, PyObject* pArgs)
{
	sys->useNN = true;
	sys->neighborlist.init();
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_run(PyObject* pSelf, PyObject* pArgs)
{
	int steps;
	if (!PyArg_ParseTuple(pArgs, "i", &steps)) return NULL;
	//printf("in run(): running %d steps \n", steps);
	sys->integrator.run(steps);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject* md_system_print_avgs(PyObject* pSelf, PyObject* pArgs)
{
	std::cout << sys->particles.pos[1] << std::endl;
	printf("%d %f %f %f %f %f\n",
		sys->steps,
		sys->props.uKinetic.avgval,
		sys->props.uPotential.avgval,
		sys->props.uTotal.avgval,
		sys->props.uPair.avgval,
		sys->props.uFea.avgval
		//sys->props.pressure.avgval,
		//sys->props.clusterUKinetic.avgval,
		//sys->props.clusterVirial.avgval,
		//sys->props.molPressure.avgval
	);

	Py_INCREF(Py_None);
	return Py_None;
}


PyObject* md_system_get_properties(PyObject* pSelf, PyObject* pArgs)
{

	//if(sys->props.step_avg == 0) sys->props.accumulate(2);
	PyObject *props_dict;
	props_dict = PyDict_New();

	PyDict_SetItemString(props_dict, "uKinetic", Py_BuildValue("(ddd)",
		sys->props.uKinetic.val, sys->props.uKinetic.avgval, sys->props.uKinetic.stdv));
	PyDict_SetItemString(props_dict, "uPotential", Py_BuildValue("(ddd)",
		sys->props.uPotential.val, sys->props.uPotential.avgval, sys->props.uPotential.stdv));
	PyDict_SetItemString(props_dict, "uTotal", Py_BuildValue("(ddd)",
		sys->props.uTotal.val, sys->props.uTotal.avgval, sys->props.uTotal.stdv));
	PyDict_SetItemString(props_dict, "pressure", Py_BuildValue("(ddd)",
		sys->props.pressure.val, sys->props.pressure.avgval, sys->props.pressure.stdv));

	return props_dict;
}

PyObject* md_system_get_avgs(PyObject* pSelf, PyObject* pArgs)
{
	return Py_BuildValue("[iddddd]",
		sys->steps,
		sys->props.uKinetic.avgval,
		sys->props.uPotential.avgval,
		sys->props.uTotal.avgval,
		sys->props.uPair.avgval,
		sys->props.uFea.avgval
		//sys->props.feaVirial.avgval,
		//sys->props.pressure.avgval,
		//sys->props.clusterUKinetic.avgval,
		//sys->props.clusterVirial.avgval,
		//sys->props.molPressure.avgval
	);
}



static PyMethodDef md_methods[] = {
	{ "system_c_main", md_system_c_main, METH_NOARGS,"main function for testing" },
	{ "system_init", md_system_init, METH_VARARGS,"Set up MD initial syste" },
	{ "system_reset", md_system_reset, METH_NOARGS,"Set up MD initial syste" },
	{ "system_set_dt", md_system_set_dt, METH_VARARGS,"Set MD dt" },
	{ "system_add_particle", md_system_add_particle, METH_VARARGS,"Add one particle to system with pos, vel, mass, sigma and type" },
	{ "system_add_cluster", md_system_add_cluster, METH_VARARGS,"Add one cluster to system with pos, vel, mass, sigma and type" },
	{ "system_set_velocities", md_system_set_velocities, METH_VARARGS,"Set initial velocities" },
	{ "system_set_zero_center_of_mass_vel", md_system_zero_center_of_mass_vel, METH_NOARGS,"zero center of mass velecities" },
	{ "system_set_potential", md_system_set_potential, METH_VARARGS,"Set potential" },
	{ "system_set_equilibrium_steps", md_system_set_equilibrium_steps, METH_VARARGS,"Set up equilibrium steps" },
	{ "system_set_avg_steps", md_system_set_avg_steps, METH_VARARGS,"Set up avg steps" },
	{ "system_run", md_system_run, METH_VARARGS,"Run system a number of steps" },
	{ "system_print_avgs", md_system_print_avgs, METH_NOARGS,"Print system avgs" },
	{ "system_get_properties", md_system_get_properties, METH_NOARGS,"Return dict with properties" },
	{ "system_get_avgs", md_system_get_avgs, METH_NOARGS,"Returns list with avg vals" },
	{ "system_get_positions", md_system_get_positions, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_get_unfolded_positions", md_system_get_unfolded_positions, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_unfold_positions", md_system_unfold_positions, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_get_velocities", md_system_get_velocities, METH_NOARGS,"Returns list with avg vals" },
	{ "system_get_accelerations", md_system_get_accelerations, METH_NOARGS,"Returns list with avg vals" },
	{ "system_get_forces", md_system_get_forces, METH_NOARGS,"Returns list with avg vals" },
	{ "system_get_tetras", md_system_get_tetras, METH_NOARGS,"Returns list with triangles" },
	{ "system_get_bonds", md_system_get_bonds, METH_NOARGS,"Returns list with bonds" },
	{ "system_get_particle_data", md_system_get_particle_data, METH_NOARGS,"Returns list with bonds" },
	{ "system_get_particle_sizes", md_system_get_particle_sizes, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_get_elements_volumes", md_system_get_elements_volumes, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_get_elements_invariants", md_system_get_elements_invariants, METH_NOARGS,"Returns numpy array with pos" },
	{ "system_init_neighbor_list", md_system_init_neighbor_list, METH_NOARGS,"init neighbor list" },
	{ "system_set_scale_factor", md_system_set_scale_factor, METH_VARARGS,"Set up scale factor" },
	{ "system_set_boundary_conditions", md_system_set_boundary_conditions, METH_VARARGS,"Set up boundary conditions" },
	{ "system_set_box", md_system_set_box, METH_VARARGS,"Set up box for system, Lx, Ly" },
	{ "system_get_box", md_system_get_box, METH_NOARGS,"Get box for system, Lx, Ly" },
	{ "system_get_walls_pos", md_system_get_walls_pos, METH_NOARGS,"Get box for system, Lx, Ly" },
	{ "system_set_walls_pos", md_system_set_walls_pos, METH_VARARGS,"Get box for system, Lx, Ly" },
	{ "system_set_moving_walls", md_system_set_moving_walls, METH_VARARGS,"Get box for system, Lx, Ly" },
	{ "system_set_walls_shear", md_system_set_walls_shear, METH_VARARGS,"Get box for system, Lx, Ly" },
	{ "system_set_walls_moving_rate", md_system_set_walls_moving_rate, METH_VARARGS,"Get box for system, Lx, Ly" },
	{ "system_get_walls_forces", md_system_get_walls_forces, METH_NOARGS,"Get box for system, Lx, Ly" },
	{ "system_set_gravity", md_system_set_gravity, METH_VARARGS,"Set up box for system, gravity" },
	{ "system_set_friction", md_system_set_friction, METH_VARARGS,"Set up box for system, gravity" },
	{ "system_set_bond_spring_constant", md_system_set_bond_spring_constant, METH_VARARGS,"Set k for bonds in system" },
	{ "system_set_area_force_constant", md_system_set_area_force_constant, METH_VARARGS,"Set k for bonds in system" },
	{ "system_set_youngs_modulus", md_system_set_youngs_modulus, METH_VARARGS,"Set k for bonds in system" },
	{ "system_set_swelling_energy_constants", md_system_set_swelling_energy_constants, METH_VARARGS,"Set k for bonds in system" },

	// Terminate the array with an object containing nulls.
	{ nullptr, nullptr, 0, nullptr }
};


static PyModuleDef md_module = {
	PyModuleDef_HEAD_INIT,
	"femmd",                        // Module name
	"Provides some functions, but faster",  // Module description
	0,
	md_methods                   // Structure that defines the methods
};

#ifdef CPYTHON_INTERFACE

PyMODINIT_FUNC PyInit_femmd() {
	PyObject* m = NULL;
	m = PyModule_Create(&md_module);
	if (!m) return NULL;
	import_array();
	return m;
}

#endif





//{ "system_init", md_system_init, METH_VARARGS, "Set up MD initial syste" },
//{ "system_reset", md_system_reset, METH_NOARGS,"Set up MD initial syste" },

/*
PYBIND11_MODULE(femmd, m) {
	m.doc() = "femmd module, molecular dynamics and finite elements code"; // optional module docstring
	
	m.def("init", &md_system_init, "Set up MD initial system");
}
*/
