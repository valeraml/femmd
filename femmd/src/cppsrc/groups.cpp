#include "groups.h"


int Groups::addGroup() {
	int ind = groups.size();
	std::vector<int> g;
	groups.push_back(g);
	len.push_back(0);
	return ind;
}

void Groups::insertParticleInGroup(int pid, int g) {
	id_of_particle[pid] = g;
	groups[g].push_back(pid);
	len[g] += 1;
}