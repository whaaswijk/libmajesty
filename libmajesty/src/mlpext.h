#ifndef MLPEXT_H
#define MLPEXT_H

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "mlp.h"

using namespace std;

extern void readVerilog(FILE *file, MIG **mig);

static inline bool has_input(MAJ3* n, const MAJ3* const input) {
	return n->in3 == input || n->in2 == input || n->in1 == input;
}

inline unsigned int maxlevel(MAJ3* n1, MAJ3* n2, MAJ3* n3) {
	if (n1->level > n2->level) {
		return (n1->level > n3->level) ? n1->level : n3->level;
	} else {
		return (n2->level > n3->level) ? n2->level : n3->level;
	}
}

vector<MAJ3*> mig_topsort(MIG *mig);

vector<unsigned int> outidx(MIG* m, MAJ3* n);
//vector<unsigned int> outidx(const xmg&, MAJ3*);

#endif
