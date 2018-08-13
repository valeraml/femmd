#ifndef __GROUPS_H__
#define __GROUPS_H__


#include "def.h"

class Groups {
public:

	std::vector< std::vector<int> > groups;
	std::vector<int> len;
	std::vector<int> id_of_particle; //array with ids of group that particle belong to

	int addGroup();
	void insertParticleInGroup(int pid, int g);
	void clear() {
		groups.clear();
		len.clear();
		id_of_particle.clear();
	};

};

#endif