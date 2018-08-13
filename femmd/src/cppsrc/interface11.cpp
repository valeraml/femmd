

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
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"
//#include "pybind11/stl_bind.h"

// Modified this file for includes to work ... FIXME
#include "pybind11/eigen.h"

//#include <Python.h>
//#include <numpy/arrayobject.h>


#include <cmath>   

namespace py = pybind11;
MDSystem *sys11;
PYBIND11_MAKE_OPAQUE(hvector<VectorR>);


void sys_set_defaults() {
	sys11->useNN = false;
	sys11->dt = 0.001;
	sys11->props.step_avg = 1000;
	sys11->props.step_equi = 0;
	sys11->particles.numTypes = 1;
	sys11->pbcType = XYPBC;
	sys11->interactions.gravity = 0;
	sys11->interactions.gravityForce = false;
	sys11->interactions.E = 5;
	sys11->interactions.kBond = 10;
	sys11->interactions.kArea = .0001;
	sys11->interactions.bondForces = false;
	sys11->interactions.feaForces = false;
	sys11->interactions.areaForces = false;
	sys11->interactions.setElasticConstants();

	double v0 = 0.5;
	//sys->particles.setParticlesVel(v0);
	//if (sys->useNN)
	//	sys->neighborlist.init();

	sys11->start = true;
}

void sys_init(VectorR box) {
	printf("Using pybind11\n");
	printf("%f %f\n", box[0], box[1]);
	sys11 = new MDSystem(0);
	sys11->init(box);
	std::cout << "Creating System @ location: "<< sys11 << std::endl;
	sys_set_defaults();
}

void sys_add_particle(VectorR r, VectorR v, double m, double sig, int typ) {
	sys11->particles.addParticle(r, v, m, sig, typ);
}


void sys_add_cluster(py::array_t<real> vert, py::array_t<int> cells, VectorR r0, int typ, real mass, int exc) {
	py::buffer_info vert_info = vert.request();
	auto vert_ptr = static_cast<real *>(vert_info.ptr);
	int nv = vert_info.shape[0];
	py::buffer_info cells_info = cells.request();
	auto cells_ptr = static_cast<int *>(cells_info.ptr);
	int nc = cells_info.shape[0];
	sys11->clusters.addCluster(vert_ptr, nv, cells_ptr, nc, r0, typ, mass, exc);
}

void sys_set_potential(py::array_t<real> eps, py::array_t<real> rcut, py::array_t<real> shift) {

	py::buffer_info eps_info = eps.request();
	auto eps_data = static_cast<real *>(eps_info.ptr);
	auto size = eps_info.size;
	py::buffer_info rcut_info = rcut.request();
	auto rcut_data = static_cast<real *>(rcut_info.ptr);
	py::buffer_info shift_info = shift.request();
	auto shift_data = static_cast<real *>(shift_info.ptr);
	for (int i = 0; i < size; i++) {
		printf("inter param %d %f %f %f\n", i, eps_data[i], rcut_data[i], shift_data[i]);
	}
	sys11->interactions.setPairForce(eps_data, rcut_data, shift_data);
}

void sys_set_gravity(real g) {
	if (g > 0) {
		sys11->interactions.gravity = g;
		sys11->interactions.gravityForce = true;
	}
	else
		sys11->interactions.gravityForce = false;
}

void sys_set_friction(real mu) {
	if (mu > 0) {
		sys11->interactions.frictionCoeff = mu;
		sys11->interactions.frictionForce = true;
	}
	else
		sys11->interactions.gravityForce = false;
}

void sys_set_bond_spring_constant(real k) {
	if (k > 0) {
		sys11->interactions.kBond = k;
		sys11->interactions.bondForces = true;
	}
	else
		sys11->interactions.bondForces = false;
}

void sys_set_youngs_modulus(double E) {
	if (E > 0) {
		sys11->interactions.E = E;
		sys11->interactions.feaForces = true;
		sys11->interactions.elasticForceType = 1;
	}
	else
		sys11->interactions.feaForces = false;
	sys11->interactions.setElasticConstants();
	printf("Elastic Constants: E=%f nu=%f mu=%f lambda=%f\n",
		sys11->interactions.E, sys11->interactions.nu, sys11->interactions.mu, sys11->interactions.lambda);
}

void sys_set_swelling_energy_constants(real we0, real wmix0, real chi) {
	if (we0 > 0) {
		sys11->interactions.We0 = we0;
		sys11->interactions.Wmix0 = wmix0;
		sys11->interactions.Chi = chi;
		sys11->interactions.feaForces = true;
		sys11->interactions.elasticForceType = 2;
	}
	else
		sys11->interactions.feaForces = false;
	sys11->interactions.setElasticConstants();
	printf("Swelling Constants: We0=%f Wmix0=%f Chi=%f\n",
		sys11->interactions.We0, sys11->interactions.Wmix0, sys11->interactions.Chi);
}

auto sys_get_positions() {
	real *ptr = (real*)sys11->particles.pos.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles, NDIM }, // shape
	{ NDIM * sizeof(real), sizeof(real) }, // C-style contiguous strides for double
	ptr, 
	dummy); // the data pointer
}

auto pyarray_from_data(hvector<VectorR> &vec, int itemdim) {
	real* v = (real*)vec.data();
	auto shape = { (int)vec.size(), itemdim };
	auto strides = { itemdim * sizeof(real), sizeof(real) };
	auto capsule = py::capsule(v, [](void *v) { /*delete reinterpret_cast<std::vector<real>*>(v); */});
	return py::array_t<real>(shape, strides, v, capsule);
}

auto sys_get_positions1() {
	//py::buffer_info info;
	//info.ptr = sys11->particles.pos.data();
	//info.itemsize = sizeof(real);
	//info.ndim = 2;
	//info.shape = { (long long)sys11->particles.vel.size(), NDIM };
	//info.strides = { sizeof(real)*NDIM, sizeof(real) };
	//info.format = py::format_descriptor<real>::value;
	//return py::array_t<real>(info);

	real *v = (real*)sys11->particles.pos.data();
	auto shape = { sys11->particles.n_particles, NDIM };
	auto strides = { NDIM * sizeof(real), sizeof(real) };
	auto capsule = py::capsule(v, [](void *v) { delete reinterpret_cast<std::vector<real>*>(v); });
	//return py::array(v->size(), v->data(), capsule);
	return py::array_t<real>(shape, strides, v, capsule);
}

auto sys_get_velocities() {

	real *ptr = (real*)sys11->particles.vel.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles, NDIM }, // shape
	{ NDIM * sizeof(real), sizeof(real) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer
}

auto sys_get_accelerations() {
	real *ptr = (real*)sys11->particles.acc.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles, NDIM }, // shape
	{ NDIM * sizeof(real), sizeof(real) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer
}

auto sys_get_forces() {
	real *ptr = (real*)sys11->particles.force.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles, NDIM }, // shape
	{ NDIM * sizeof(real), sizeof(real) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer
}

auto sys_get_unfolded_positions() {

	real *ptr = (real*)sys11->elements.unfoldedPos.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles, NDIM }, // shape
	{ NDIM * sizeof(real), sizeof(real) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer
}

auto sys_get_tetras() {
	//py::buffer_info info;
	//info.ptr = sys11->elements.tetras.data();
	//info.itemsize = sizeof(int);
	//info.ndim = 2;
	//info.shape = { (long long)sys11->elements.tetras.size(), NDIM + 1 };
	//info.strides = { sizeof(int)*(NDIM+1), sizeof(int) };
	//info.format = py::format_descriptor<int>::value;
	//return py::array(info);

	int* ptr = (int*)sys11->elements.tetras.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<int>(
	{ sys11->elements.numTetras , NDIM + 1 }, // shape
	{ (NDIM+1) * sizeof(int), sizeof(int) }, // C-style contiguous strides for double
	ptr,
	dummy); // the data pointer

}

auto sys_get_bonds() {
	int* ptr = (int*)sys11->bonds.bondList.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<int>(
	{ sys11->bonds.bondListLen/2 , 2 }, // shape
	{ 2 * sizeof(int), sizeof(int) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer

}

auto sys_get_particle_sizes() {
	real* ptr = (real*)sys11->particles.sigma.data();
	py::capsule dummy(ptr, [](void *f) {});
	return py::array_t<real>(
	{ sys11->particles.n_particles }, // shape
	{ sizeof(real) }, // C-style contiguous strides for double
		ptr,
		dummy); // the data pointer
}

auto sys_get_elements_volumes() {
	real* ptr1 = (real*)sys11->elements.refVol.data();
	real* ptr2 = (real*)sys11->elements.currVol.data();
	py::capsule dummy1(ptr1, [](void *f) {});
	py::capsule dummy2(ptr2, [](void *f) {});

	py::array_t<real> refVol = py::array_t<real>(
	{ sys11->elements.numTetras }, // shape
	{ sizeof(real) }, // C-style contiguous strides for double
		ptr1,
		dummy1); // the data pointer

	py::array_t<real> currVol = py::array_t<real>(
	{ sys11->elements.numTetras }, // shape
	{ sizeof(real) }, // C-style contiguous strides for double
		ptr2,
		dummy1); // the data pointer

	return py::make_tuple(refVol, currVol);
}

auto sys_get_elements_invariants() {
	real* ptr1 = (real*)sys11->elements.I1.data();
	real* ptr2 = (real*)sys11->elements.I2.data();
	py::capsule dummy1(ptr1, [](void *f) {});
	py::capsule dummy2(ptr2, [](void *f) {});

	py::array_t<real> I1 = py::array_t<real>(
	{ sys11->elements.numTetras }, // shape
	{ sizeof(real) }, // C-style contiguous strides for double
		ptr1,
		dummy1); // the data pointer

	py::array_t<real> I2 = py::array_t<real>(
	{ sys11->elements.numTetras }, // shape
	{ sizeof(real) }, // C-style contiguous strides for double
		ptr2,
		dummy1); // the data pointer

	return py::make_tuple(I1, I2);
}

py::list sys_get_avgs()
{
	return py::make_tuple(
		sys11->steps,
		sys11->props.uKinetic.avgval,
		sys11->props.uPotential.avgval,
		sys11->props.uTotal.avgval,
		sys11->props.uPair.avgval,
		sys11->props.uFea.avgval
		//sys->props.feaVirial.avgval,
		//sys->props.pressure.avgval,
		//sys->props.clusterUKinetic.avgval,
		//sys->props.clusterVirial.avgval,
		//sys->props.molPressure.avgval
	);
}

void sys_print_avgs()
{
	printf("%d %f %f %f %f %f\n",
		sys11->steps,
		sys11->props.uKinetic.avgval,
		sys11->props.uPotential.avgval,
		sys11->props.uTotal.avgval,
		sys11->props.uPair.avgval,
		sys11->props.uFea.avgval
		//sys11->props.pressure.avgval,
		//sys11->props.clusterUKinetic.avgval,
		//sys11->props.clusterVirial.avgval,
		//sys11->props.molPressure.avgval
	);
}

#ifdef PYBIND11_INTERFACE

PYBIND11_MODULE(femmd, m) {
	m.doc() = "femmd module, molecular dynamics and finite elements code"; // optional module docstring
	
	m.def("system_init", &sys_init, "Set up MD initial system");
	m.def("system_reset",  []() { sys11->reset(); }, "reset MD system");
	m.def("system_set_dt", [](real dt) {sys11->dt = dt; });
	m.def("system_get_dt", []() {return sys11->dt; });
	m.def("system_add_particle", &sys_add_particle);
	m.def("system_add_cluster", &sys_add_cluster);
	m.def("system_set_potential", &sys_set_potential);
	m.def("system_set_velocities", [](real v0) {sys11->particles.setParticlesVel(v0); });
	m.def("system_set_zero_center_of_mass_vel", []() {sys11->particles.zeroCM(); });
	m.def("system_set_gravity", &sys_set_gravity);
	m.def("system_set_friction", &sys_set_friction);
	m.def("system_set_bond_spring_constant", &sys_set_bond_spring_constant);
	m.def("system_set_youngs_modulus", &sys_set_youngs_modulus);
	m.def("system_set_swelling_energy_constants", &sys_set_swelling_energy_constants);
	m.def("system_set_equi_steps", [](int s) { sys11->props.step_equi = s; });
	m.def("system_set_avg_steps", [](int s) { sys11->props.step_avg = s; });
	m.def("system_set_box", [](VectorR b) { sys11->setBox(b); });
	m.def("system_get_box", []() {return sys11->box; });
	m.def("system_get_walls_pos", []() { return py::make_tuple(sys11->walls.pos[0],sys11->walls.pos[1]); });
	//m.def("system_get_walls_forces", []() { return sys11->walls.forces; });
	m.def("system_get_walls_forces", []() { return py::make_tuple(sys11->walls.forces[0], sys11->walls.forces[1]); });
	m.def("system_set_walls_pos", [](VectorR p0, VectorR p1) {sys11->walls.pos[0] = p0; sys11->walls.pos[1] = p1; });
	m.def("system_set_walls_moving_rate", [](VectorR r0, VectorR r1) 
		{sys11->walls.motion_rate[0] = r0; sys11->walls.motion_rate[1] = r1; });
	m.def("system_set_moving_walls", [](int s) { sys11->walls.moveWalls = s; });
	m.def("system_set_boundary_conditions", [](int s) { sys11->pbcType = (PBCTYPE)s; });
	m.def("system_get_positions1", &sys_get_positions1, py::return_value_policy::reference);
	m.def("system_get_positions", &sys_get_positions);
	m.def("system_get_velocities", &sys_get_velocities);
	m.def("system_get_accelerations", &sys_get_accelerations);
	m.def("system_get_forces", &sys_get_forces);
	m.def("system_get_unfolded_positions", &sys_get_unfolded_positions);
	m.def("system_get_tetras", &sys_get_tetras);
	m.def("system_get_bonds", &sys_get_bonds);
	m.def("system_get_particle_sizes", &sys_get_particle_sizes);
	m.def("system_get_elements_volumes", &sys_get_elements_volumes);
	m.def("system_get_elements_invariants", &sys_get_elements_invariants);
	m.def("system_unfold_positions", []() { sys11->elements.unfoldPos(sys11->box); });
	m.def("system_init_neighbor_list", []() {sys11->useNN = true; sys11->neighborlist.init(); });
	m.def("system_run", [](int s) {sys11->integrator.run(s); });
	m.def("system_get_avgs", &sys_get_avgs);
	m.def("system_print_avgs", &sys_print_avgs);
}

#endif
