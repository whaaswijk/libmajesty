#include <iostream>
#include <time.h>
#include "mlpext.h"
#include "xmg.h"
#include <queue>
#include <cassert>
#include <ctime>

#define set_prop_flag(n) n->flag |= 1
#define clear_prop_flag(n) n->flag &= (~1)
#define is_prop_set(n) (n->flag & 1)
#define set_rel_flag(n) n->flag |= 2
#define clear_rel_flag(n) n->flag &= (~2)
#define is_rel_set(n) (n->flag & 2)

// Sort the nodes of the MIG in topological order.
vector<MAJ3*> mig_topsort(MIG *mig) {
	resetflag(mig);

	vector<MAJ3*> v;
	queue<MAJ3*> q;

	for (auto i = 0u; i < mig->Nnodes; i++) {
		mig->nodes[i]->flag = 3;
	}

	q.push(mig->one);
	for (auto i = 0u; i < mig->Nin; i++) {
		q.push(mig->in[i]);
	}

	while (!q.empty()) {
		auto node = q.front(); q.pop();
		for (auto i = 0u; i < node->fanout; i++) {
			auto outnode = node->outEdges[i];
			--outnode->flag;
			if (outnode->flag == 0) {
				q.push(outnode);
			}
		}
		v.push_back(node);
	}

	resetflag(mig);
	return v;
}

static pair<MAJ3*,bool>
new_node(MAJ3* in1, bool c1, MAJ3* in2, bool c2, MAJ3* in3, bool c3, vector<MAJ3*>& nn) {
	if (in1 == in2) {
		if (c1 == c2)
			return make_pair(in1, c1);
		else
			return make_pair(in3, c3);
	} else if (in1 == in3) {
		if (c1 == c3)
			return make_pair(in1, c1);
		else
			return make_pair(in2, c2);
	} else if (in2 == in3) {
		if (c2 == c3)
			return make_pair(in2, c2);
		else
			return make_pair(in1, c1);
	}
	MAJ3* n = (MAJ3*)malloc(sizeof(MAJ3));
	n->in1 = in1;
	n->compl1 = c1;
	n->in2 = in2;
	n->compl2 = c2;
	n->in3 = in3;
	n->compl3 = c3;
	n->level = maxlevel(in1, in2, in3) + 1;
	n->flag = n->label = n->value = n->fanout = n->PI = n->PO = 0u;
	n->aux = NULL;
	n->outEdges = NULL;
	addfanout(in1, n);
	addfanout(in2, n);
	addfanout(in3, n);
	nn.push_back(n);
	return make_pair(n, false);
}

inline vector<pair<MAJ3*,bool>> select_lowest(MAJ3* n) {
	vector<pair<MAJ3*,bool>> buf;
	if (n->in1->level <= n->in2->level && n->in1->level <= n->in3->level)
		buf.push_back(make_pair(n->in1, n->compl1));
	if (n->in2->level <= n->in1->level && n->in2->level <= n->in3->level)
		buf.push_back(make_pair(n->in2, n->compl2));
	if (n->in3->level <= n->in1->level && n->in3->level <= n->in2->level)
		buf.push_back(make_pair(n->in3, n->compl3));
	
	return buf;
}

inline pair<MAJ3*,bool> select_highest(MAJ3* n) {
	if (n->in1->level > n->in2->level && n->in1->level > n->in3->level) {
		return make_pair(n->in1, n->compl1) ;
	} else if (n->in2->level > n->in1->level && n->in2->level > n->in3->level) {
		return make_pair(n->in2, n->compl2);
	} else if (n->in3->level > n->in1->level && n->in3->level > n->in2->level) {
		return make_pair(n->in3, n->compl3);
	} else {
		return make_pair(nullptr, false);
	}
}

inline pair<MAJ3*,bool> drop_two(MAJ3* n, MAJ3* one, MAJ3* two) {
	if (n->in1 != one && n->in1 != two) {
		return make_pair(n->in1, n->compl1);
	} else if (n->in2 != one && n->in2 != two) {
		return make_pair(n->in2, n->compl2);
	} else if (n->in3 != one && n->in3 != two) {
		return make_pair(n->in3, n->compl3);
	} else {
		assert(false);
	}
}


vector<unsigned int> outidx(MIG* m, MAJ3* n) {
	vector<unsigned int> v;
	for (auto i = 0u; i < m->Nout; i++) {
		if (m->out[i] == n) {
			v.push_back(i);
		}
	}
	return v;
}

