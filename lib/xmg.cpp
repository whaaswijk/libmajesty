#include <xmg.h>
#include "strashmap.h"
#include "primes.h"
#include <iostream>
#include <random>
#include "mlpext.h"
#include <functional>
#include <boost/optional.hpp>
#include "mlpext.h"
#include <sstream>
#include <maj_io.h>
#include <minisat/Solver.h>
#include <minisat/SolverTypes.h>

using namespace std;
using boost::optional;

static default_random_engine generator;
static uniform_int_distribution<unsigned int>
distribution(1,numeric_limits<unsigned int>::max());
static auto rvector = bind(distribution, generator);

MAJ3* hashnode(unsigned int h, hashmap& hnmap) {
	auto it = hnmap.find(h);
	if (it != hnmap.end()) {
		return it->second;
	} else {
		return nullptr;
	}
}

optional<nodeid> hashnode(unsigned int h, xhashmap& hnmap) {
	optional<nodeid> res;
	auto it = hnmap.find(h);
	if (it != hnmap.end()) {
		res = it->second;
	}
	return res;
}

bv*
new_bv(const bv* v1, bool c1, const bv* v2, bool c2, const bv* v3, bool c3) {
	assert(v1->size() == v2->size());
	assert(v1->size() == v3->size());
	auto res = new bv(v1->size());
	for (auto i = 0u; i < v1->size(); ++i) {
		auto b1 = v1->at(i); if (c1) b1 = ~b1;
		auto b2 = v2->at(i); if (c2) b2 = ~b2;
		auto b3 = v3->at(i); if (c3) b3 = ~b3;
		res->at(i) = (b1 & b2) | (b1 & b3) | (b2 & b3);
	}

	return res;
}

bv*
new_bv(const bv* v1, bool c1, const bv* v2, bool c2) {
	assert(v1->size() == v2->size());
	auto res = new bv(v1->size());
	for (auto i = 0u; i < v1->size(); ++i) {
		auto b1 = v1->at(i); if (c1) b1 = ~b1;
		auto b2 = v2->at(i); if (c2) b2 = ~b2;
		res->at(i) = (b1 ^ b2);
	}

	return res;
}


unsigned int hashbv(const bv* v) {
	assert(v->size() <= PRIMES.size());
	auto val = 0u;

	for (auto i = 0u; i < v->size(); i++) {
		val ^= v->at(i) * PRIMES[i];
	}

	return val;
}

unsigned int chashbv(const bv* v) {
	assert(v->size() <= PRIMES.size());
	auto val = 0u;

	for (auto i = 0u; i < v->size(); i++) {
		val ^= (~v->at(i)) * PRIMES[i];
	}

	return val;
}

inline size_t bv_size(unsigned int nbits) {
	auto bits_in_int = sizeof(unsigned int)*8;
	auto nints = ceil((1.0*nbits)/bits_in_int);
	return static_cast<unsigned int>(nints);
}

inline void randomize(bv* v) {
	for (auto i = 0u; i < v->size(); ++i) {
		v->at(i) = rvector();
	}
}

void
inithashmap(bvmap& bv_map, hashmap& hnmap, vector<MAJ3*>& inputs, MAJ3* one) {
	assert(hnmap.size() == 0);
	auto ov = bv_map[one];
	auto oh = hashbv(ov.get());
	hnmap[oh] = one;
	for (auto input : inputs) {
		auto v = bv_map[input];
		auto h = hashbv(v.get());
		hnmap[h] = input;
	}
}

void inithashmap(xbvmap& bv_map, xhashmap& hnmap, vector<nodeid>& inputs) {
	assert(hnmap.size() == 0);
	auto ov = bv_map[0];
	auto oh = hashbv(ov.get());
	hnmap[oh] = 0;
	for (auto input : inputs) {
		auto v = bv_map[input];
		auto h = hashbv(v.get());
		hnmap[h] = input;
	}
}

void initsim(bvmap& bv_map, hashmap& hnmap, vector<MAJ3*>& inputs, MAJ3* one) {
	auto size = bv_size(BITVECTOR_SIZE);
	shared_ptr<bv> ov(new bv(size));
	for (auto i = 0u; i < size; i++) {
		ov->at(i) = -1;
	}
	bv_map[one] = ov;
	for (auto input : inputs) {
		shared_ptr<bv> v(new bv(size));
		randomize(v.get());
		bv_map[input] = v;
	}
	inithashmap(bv_map, hnmap, inputs, one);
}

void initsim(xbvmap& bv_map, xhashmap& hnmap, vector<nodeid>& inputs) {
	auto size = bv_size(BITVECTOR_SIZE);
	shared_ptr<bv> ov(new bv(size));
	for (auto i = 0u; i < size; i++) {
		ov->at(i) = -1;
	}
	bv_map[0] = ov;
	for (auto input : inputs) {
		shared_ptr<bv> v(new bv(size));
		randomize(v.get());
		bv_map[input] = v;
	}
	inithashmap(bv_map, hnmap, inputs);
}

bv*
find_or_create_bv(const bv* v1, bool c1, const bv* v2, bool c2, const bv* v3, 
		bool c3, MAJ3* node, bvmap& bvmap) {
	auto it = bvmap.find(node);
	if (it == bvmap.end()) {
		auto v = shared_ptr<bv>(new_bv(v1, c1, v2, c2, v3, c3));
		bvmap[node] = v;
		return v.get();
	}
	auto v = it->second.get();
	assert(v->size() == (v1->size()-1));
	auto b1 = v1->at(v->size()); if (c1) b1 = ~b1;
	auto b2 = v2->at(v->size()); if (c2) b2 = ~b2;
	auto b3 = v3->at(v->size()); if (c3) b3 = ~b3;
	v->push_back((b1 & b2) | (b1 & b3) | (b2 & b3));
	return v;
}

bv*
find_or_create_bv(const bv* v1, bool c1, const bv* v2, bool c2, const bv* v3, 
		bool c3, nodeid nodeid, xbvmap& bvmap) {
	auto v = bvmap.at(nodeid);
	if (v == nullptr) {
		auto v = shared_ptr<bv>(new_bv(v1, c1, v2, c2, v3, c3));
		bvmap[nodeid] = v;
		return v.get();
	}
	assert(v->size() == (v1->size()-1));
	auto b1 = v1->at(v->size()); if (c1) b1 = ~b1;
	auto b2 = v2->at(v->size()); if (c2) b2 = ~b2;
	auto b3 = v3->at(v->size()); if (c3) b3 = ~b3;
	v->push_back((b1 & b2) | (b1 & b3) | (b2 & b3));
	return v.get();
}

bv*
find_or_create_bv(const bv* v1, bool c1, const bv* v2, bool c2, 
		nodeid nodeid, xbvmap& bvmap) {
	auto v = bvmap.at(nodeid);
	if (v == nullptr) {
		v = shared_ptr<bv>(new_bv(v1, c1, v2, c2));
		bvmap[nodeid] = v;
		return v.get();
	}
	assert(v->size() == (v1->size()-1));
	auto b1 = v1->at(v->size()); if (c1) b1 = ~b1;
	auto b2 = v2->at(v->size()); if (c2) b2 = ~b2;
	v->push_back(b1 ^ b2); 
	return v.get();
}

void
simulate(bvmap& m, hashmap& hnmap, simmap& simrep, vector<MAJ3*>& nodes) {
	for (auto i = 0u; i < nodes.size(); i++) {
		auto node = nodes[i];
		if (node->PI) {
			continue;
		}
		auto bv = find_or_create_bv(
			m[node->in1].get(), static_cast<bool>(node->compl1),
			m[node->in2].get(), static_cast<bool>(node->compl2),
			m[node->in3].get(), static_cast<bool>(node->compl3), node, m);
		auto hash = hashbv(bv);
		auto hashc = chashbv(bv);

		auto hnode = hashnode(hash, hnmap);
		if (hnode != nullptr) {
			simrep[node] = hnode;
			continue;
		}
		hnode = hashnode(hashc, hnmap);
		if (hnode != nullptr) {
			simrep[node] = hnode;
		} else {
			hnmap[hash] = node;
			simrep[node] = node;
		}
	}
}

void 
simulate(xbvmap& m, xhashmap& hnmap, xsimmap& simrep, const majesty::xmg& xmg) {
	const auto& nodes = xmg.nodes();
	for (auto i = 0u; i < nodes.size(); i++) {
		auto node = nodes[i];
		if (is_pi(node)) {
			continue;
		}
		bv* bv;
		if (is_xor(node)) {
			bv = find_or_create_bv(m[node.in1].get(), is_c1(node),
					m[node.in2].get(), is_c1(node), i, m);
		} else {
			bv = find_or_create_bv(m[node.in1].get(), is_c1(node),
					m[node.in2].get(), is_c2(node),
					m[node.in3].get(), is_c3(node), i, m);
		}
		auto hash = hashbv(bv);
		auto hashc = chashbv(bv);

		auto hnode = hashnode(hash, hnmap);
		if (hnode) {
			simrep[i] = *hnode;
			continue;
		}
		hnode = hashnode(hashc, hnmap);
		if (hnode) {
			simrep[i] = *hnode;
		} else {
			hnmap[hash] = i;
			simrep[i] = i;
		}
	}
}

using namespace Minisat;

namespace majesty {

	lbool equivalent(nodeid i1, nodeid i2, bool c, unsigned nr_backtracks, 
			const varmap& v, Solver& solver) {
		auto lit1 = mkLit(v.at(i1), false);
		auto lit2 = mkLit(v.at(i2), c ? true : false);
		auto onelit = mkLit(v.at(0), false);

		vec<Lit> assumps;
		assumps.push(~lit1);
		assumps.push(lit2);
		assumps.push(onelit);

		// Check if both variables could ever be different
		auto res = solver.frsolve(assumps, nr_backtracks);
		if (res == l_True) {
			return l_False;
		} else if (res == l_Undef) {
			return l_Undef;
		}
		assumps[0] = lit1;
		assumps[1] = ~lit2;
		res = solver.frsolve(assumps, nr_backtracks);
		if (res == l_True) {
			return l_False;
		} else if (res == l_Undef) {
			return l_Undef;
		} else {
			return l_True;
		}
	}

	boost::dynamic_bitset<> counterexample(unsigned nin, varmap& v, Solver& s) {
		boost::dynamic_bitset<> ctr;
		for (auto i = 1u; i <= nin; i++) {
			ctr.push_back(s.modelValue(v[i]) == l_True ? true : false);
		}
		return ctr;
	}

	// Append the SAT solver's counterexample
	void appendcounter(xbvmap& m, const boost::dynamic_bitset<>& counter, 
			vector<nodeid>& inputs) {
		for (auto i = 0u; i < inputs.size(); i++) {
			auto input = inputs[i];
			auto bv = m[input];
			auto newi = rvector();
			newi &= ~1;
			if (counter[i]) {
				newi |= 1;
			}
			bv->push_back(newi);
		}
		auto ov = m[0];
		ov->push_back(-1);
	}

	// Append the SAT solver's counterexample
	void appendcounter(bvmap& m, const boost::dynamic_bitset<>& counter, 
			vector<MAJ3*>& inputs, MAJ3* one) {
		for (auto i = 0u; i < inputs.size(); i++) {
			auto input = inputs[i];
			auto bv = m[input];
			auto newi = rvector();
			newi &= ~1;
			if (counter[i]) {
				newi |= 1;
			}
			bv->push_back(newi);
		}
		auto ov = m[one];
		ov->push_back(-1);
	}

	unsigned xmg::nin() const {
		auto res = 0u;
		for (auto i = 1u; i < _nodes.size(); i++) { // Skip one node
			const auto& n = _nodes[i];
			if (is_pi(n)) {
				++res;
			} else {
				break;
			}
		}
		return res;
	}

	int xmg::depth() const {
		auto max_depth = -1;
		vector<int> depth(_nodes.size());
		for (auto i = 0u; i < _nodes.size(); i++) {
			const auto& n = _nodes[i];
			if (is_pi(n)) {
				depth[i] = 0;
			} else {
				const auto in1_depth = depth[n.in1];
				const auto in2_depth = depth[n.in2];
				const auto in3_depth = depth[n.in3];
				auto max_in = (in1_depth > in2_depth) ? in1_depth : in2_depth;
				max_in = (max_in > in3_depth) ? max_in : in3_depth;
				depth[i] = max_in + 1;
			}
		}
		for (auto i = 0u; i < _nodes.size(); i++) {
			if (depth[i] > max_depth)
				max_depth = depth[i];
		}
		return max_depth;
	}

	nodeid xmg::create_input() {
		node in;
		in.flag = in.in1 = in.in2 = in.in3 = 0;
		in.ecnext = EC_NULL;
		set_pi(in);
		auto idx = _nodes.size();
		in.ecrep = idx;
		_nodes.push_back(in);
		return idx;
	}

	nodeid xmg::create_input(bool c) {
		node in;
		in.flag = in.in1 = in.in2 = in.in3 = 0;
		in.ecnext = EC_NULL;
		set_pi(in);
		if (c) {
			set_c1(in);
		}
		auto idx = _nodes.size();
		in.ecrep = idx;
		_nodes.push_back(in);
		return idx;
	}

	nodeid xmg::create_input(const string& name) {
		node in;
		in.flag = in.in1 = in.in2 = in.in3 = 0;
		in.ecnext = EC_NULL;
		set_pi(in);
		auto idx = _nodes.size();
		in.ecrep = idx;
		_nodes.push_back(in);
		_innames.push_back(name);
		return idx;
	}

	nodeid xmg::create_input(varmap& var_map, Solver& solver) {
		node in;
		in.flag = in.in1 = in.in2 = in.in3 = 0;
		in.ecnext = EC_NULL;
		set_pi(in);
		auto idx = _nodes.size();
		in.ecrep = idx;
		_nodes.push_back(in);
		// Make a SAT variable for this input
		var_map[idx] = solver.newVar();
		return idx;
	}

	nodeid xmg::create_input(const string& name, varmap& var_map, Solver& solver) {
		node in;
		in.flag = in.in1 = in.in2 = in.in3 = 0;
		in.ecnext = EC_NULL;
		set_pi(in);
		auto idx = _nodes.size();
		in.ecrep = idx;
		_nodes.push_back(in);
		_innames.push_back(name);
		// Make a SAT variable for this input
		var_map[idx] = solver.newVar();
		return idx;
	}

	nodeid xmg::create_node(maj3inputs) {
		node n;
		n.flag = 0;
		n.in1 = in1; if (c1) { set_c1(n); }
		n.in2 = in2; if (c2) { set_c2(n); }
		n.in3 = in3; if (c3) { set_c3(n); }
		n.ecnext = EC_NULL;
		auto idx = _nodes.size();
		n.ecrep = idx;
		_nodes.push_back(n);

		return idx;
	}

	inline void add_clauses(maj3inputs, nodeid idx, varmap& var_map, Solver& solver) {
		// Add SAT variable and clauses for this node
		auto new_var = solver.newVar();
		auto nlit = mkLit(new_var, false);
		var_map[idx] = new_var;

		auto lit1 = mkLit(var_map[in1], c1);
		auto lit2 = mkLit(var_map[in2], c2);
		auto lit3 = mkLit(var_map[in3], c3);

		// Otherwise: generate MAJ3 clauses
		// D = MAJ(A,B,C) ==>
		// (¬a||¬b||d)&&(¬a||¬c||d)&&(a||b||¬d)&&(a||c||¬d)&&(¬b||¬c||d)&&(b||c||¬d)
		solver.addClause(~lit1, ~lit2, nlit);
		solver.addClause(~lit1, ~lit3, nlit);
		solver.addClause(lit1, lit2, ~nlit);
		solver.addClause(lit1, lit3, ~nlit);
		solver.addClause(~lit2, ~lit3, nlit);
		solver.addClause(lit2, lit3, ~nlit);
	}

	nodeid xmg::create_node(maj3inputs, 
			varmap& var_map, Solver& solver, fanoutmap& f) {
		node n;
		n.flag = 0;
		n.in1 = in1; if (c1) { set_c1(n); }
		n.in2 = in2; if (c2) { set_c2(n); }
		n.in3 = in3; if (c3) { set_c3(n); }
		n.ecnext = EC_NULL;
		auto idx = _nodes.size();
		n.ecrep = idx;
		_nodes.push_back(n);

		f[in1].push_back(idx);
		f[in2].push_back(idx);
		f[in3].push_back(idx);

		add_clauses(in1, c1, in2, c2, in3, c3, idx, var_map, solver);
		
		return idx;
	}

	nodeid xmg::create_node(xorinputs) {
		node n;
		n.flag = 0;
		set_xor(n);
		n.in1 = in1; if (c1) { set_c1(n); }
		n.in2 = in2; if (c2) { set_c2(n); }
		n.ecnext = EC_NULL;
		auto idx = _nodes.size();
		n.ecrep = idx;
		_nodes.push_back(n);

		return idx;
	}
	
	inline void add_clauses(xorinputs, nodeid idx, varmap& var_map, Solver& solver) {
		// Add SAT variable and clauses for this node
		auto new_var = solver.newVar();
		auto nlit = mkLit(new_var, false);
		var_map[idx] = new_var;

		auto lit1 = mkLit(var_map[in1], c1);
		auto lit2 = mkLit(var_map[in2], c2);

		// Generate XOR clauses
		solver.addClause(~lit1, ~lit2, ~nlit);
		solver.addClause(lit1, lit2, ~nlit);
		solver.addClause(lit1, ~lit2, nlit);
		solver.addClause(~lit1, lit2, nlit);
	}

	nodeid xmg::create_node(xorinputs, varmap& var_map, 
			Minisat::Solver& solver, fanoutmap& f) {
		node n;
		n.flag = 0;
		set_xor(n);
		n.in1 = in1; if (c1) { set_c1(n); }
		n.in2 = in2; if (c2) { set_c2(n); }
		n.ecnext = EC_NULL;
		auto idx = _nodes.size();
		n.ecrep = idx;
		_nodes.push_back(n);

		f[in1].push_back(idx);
		f[in2].push_back(idx);

		add_clauses(in1, c1, in2, c2, idx, var_map, solver);

		return idx;

	}


	pair<nodeid,bool> xmg::create(xorinputs) {
		sort_inputs(in1, c1, in2, c2);
		return make_pair(create_node(in1, c1, in2, c2), false);
	}

	pair<nodeid, bool> xmg::create(input in1, input in2) {
		return create(in1.first, in1.second, in2.first, in2.second);
	}

	pair<nodeid,bool> xmg::find_or_create(xorinputs, strashmap& shmap) {
		if (in1 == in2) {
			return make_pair(0u, c1 == c2);
		} else if (in1 == 0) {
			return make_pair(in2, c1 == c2);
		} else if (in2 == 0) {
			return make_pair(in1, c1 == c2);
		}
		// Check if an trivially equivalent node already exists
		// (strashing)
		// In the case of XOR nodes we can ensure that inputs are
		// never complemented
		if (c1 && !c2) {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this),
					true);
		} else if (!c1 && c2) {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this),
					true);
		} else {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this),
					false);
		}
	}

	// Returns the index of a node corresponding to the
	// specified inputs.
	// Performs XOR propagation and 1-level strashing.
	pair<nodeid,bool> 
		xmg::find_or_create(xorinputs, strashmap& shmap, 
				varmap& v, Minisat::Solver& s, fanoutmap& f) {
		if (in1 == in2) {
			return make_pair(0u, c1 == c2);
		} else if (in1 == 0) {
			return make_pair(in2, c1 == c2);
		} else if (in2 == 0) {
			return make_pair(in1, c1 == c2);
		}
		// Check if an trivially equivalent node already exists
		// (strashing)
		// In the case of XOR nodes we can ensure that inputs are
		// never complemented
		if (c1 && !c2) {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this, v, s, f), 
					true);
		} else if (!c1 && c2) {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this, v, s, f), 
					true);
		} else {
			return make_pair(
					shmap.find_or_add(in1, false, in2, false, *this, v, s, f), 
					false);
		}
	}

	pair<nodeid,bool> xmg::create(maj3inputs) {
		sort_inputs(in1, c1, in2, c2, in3, c3);
		return make_pair(create_node(in1, c1, in2, c2, in3, c3), false);
	}
	
	pair<nodeid, bool> xmg::create(input in1, input in2, input in3) {
		return create(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second);
	}

	pair<nodeid,bool> xmg::prop_create(maj3inputs) {
		sort_inputs(in1, c1, in2, c2, in3, c3);
		if (in1 == in2) {
			if (c1 == c2) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in3, c3);
			}
		} else if (in1 == in3) {
			if (c1 == c3) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in2, c2);
			}
		} else if (in2 == in3) {
			if (c2 == c3) {
				return make_pair(in2, c2);
			} else {
				return make_pair(in1, c1);
			}
		}
		// Ensure that at most one input is complemented
		if ((c1 && !c2 && !c3) || (!c1 && c2 && !c3) || (!c1 && !c2 && c3)
				|| (!c1 && !c2 && !c3)) {
			return make_pair(create_node(in1, c1, in2, c2, in3, c3), false);
		} else {
			return make_pair(create_node(in1, c1 != true, in2, c2 != true, in3, c3 != true), true);
		}
	}
	
	pair<nodeid, bool> xmg::find_or_create(pair<nodeid,bool> in1, pair<nodeid,bool> in2, pair<nodeid,bool> in3, strashmap& shmap) {
		return find_or_create(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second, shmap);
	}

	pair<nodeid,bool> xmg::find_or_create(maj3inputs, strashmap& shmap) {
		if (in1 == in2) {
			if (c1 == c2) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in3, c3);
			}
		} else if (in1 == in3) {
			if (c1 == c3) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in2, c2);
			}
		} else if (in2 == in3) {
			if (c2 == c3) {
				return make_pair(in2, c2);
			} else {
				return make_pair(in1, c1);
			}
		}
		// Ensure that at most one input is complemented
		if ((c1 && !c2 && !c3) || (!c1 && c2 && !c3) || (!c1 && !c2 && c3)
				|| (!c1 && !c2 && !c3)) {
			auto idx = 
				shmap.find_or_add(in1, c1, in2, c2, in3, c3, *this);
			return make_pair(idx, false);
		} else {
			auto idx = shmap.find_or_add(in1, c1 != true, 
					in2, c2 != true, in3, c3 != true, *this);
			return make_pair(idx, true);
		}
	}

	// Strashe without ensuring at most one input is complemented
	pair<nodeid, bool> xmg::find_or_create_no_compl(maj3inputs, strashmap& shmap) {
		if (in1 == in2) {
			if (c1 == c2) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in3, c3);
			}
		} else if (in1 == in3) {
			if (c1 == c3) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in2, c2);
			}
		} else if (in2 == in3) {
			if (c2 == c3) {
				return make_pair(in2, c2);
			} else {
				return make_pair(in1, c1);
			}
		}
		auto idx = shmap.find_or_add(in1, c1, in2, c2, in3, c3, *this);
		return make_pair(idx, false);
	}

	pair<nodeid, bool> xmg::find_or_create_no_compl(pair<nodeid, bool> in1, pair<nodeid, bool> in2, pair<nodeid, bool> in3, strashmap& shmap) {
		return find_or_create_no_compl(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second, shmap);
	}
	
	pair<nodeid, bool> xmg::find_or_create_no_prop(pair<nodeid, bool> in1, pair<nodeid, bool> in2, pair<nodeid, bool> in3, strashmap& shmap) {
		return find_or_create_no_prop(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second, shmap);
	}

	// Returns a strashed node without applying MAJ3 propagation rules
	pair<nodeid, bool> xmg::find_or_create_no_prop(maj3inputs, strashmap& shmap) {
		auto idx =
			shmap.find_or_add(in1, c1, in2, c2, in3, c3, *this);
		return make_pair(idx, false);
	}

	pair<nodeid, bool> xmg::prop_create(input in1, input in2, input in3) {
		return prop_create(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second);
	}

	// Returns the index of a node corresponding to the
	// specified inputs.
	// Performs M_3 propagation and 1-level strashing.
	pair<nodeid,bool> xmg::find_or_create(maj3inputs, strashmap& shmap, 
			varmap& v, Solver& s, fanoutmap& f) {
		if (in1 == in2) {
			if (c1 == c2) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in3, c3);
			}
		} else if (in1 == in3) {
			if (c1 == c3) {
				return make_pair(in1, c1);
			} else {
				return make_pair(in2, c2);
			}
		} else if (in2 == in3) {
			if (c2 == c3) {
				return make_pair(in2, c2);
			} else {
				return make_pair(in1, c1);
			}
		}
		// Ensure that at most one input is complemented
		if ((c1 && !c2 && !c3) || (!c1 && c2 && !c3) || (!c1 && !c2 && c3)
				|| (!c1 && !c2 && !c3)) {
			auto idx = 
				shmap.find_or_add(in1, c1, in2, c2, in3, c3, *this, v, s, f);
			return make_pair(idx, false);
		} else {
			auto idx = shmap.find_or_add(in1, c1 != true, 
					in2, c2 != true, in3, c3 != true, *this, v, s, f);
			return make_pair(idx, true);
		}
	}

	bool xmg::xorrecinfan(nodeid i1, nodeid i2) {
		auto& n2 = _nodes[i2]; 
		if (is_flag_set(n2)) {
			return false;
		}
		set_flag(n2);
		if (i1 == i2) {
			return true;
		} 

		if (is_pi(n2)) {
			return false;
		} else if (ecrecinfan(i1, n2.in1)) {
			return true;
		} else if (ecrecinfan(i1, n2.in2)) {
			return true;
		} else {
			return false;
		}
	}

	bool xmg::majrecinfan(nodeid i1, nodeid i2) {
		auto& n2 = _nodes[i2]; 
		if (is_flag_set(n2)) {
			return false;
		}
		set_flag(n2);
		if (i1 == i2) {
			return true;
		} 

		if (is_pi(n2)) {
			return false;
		} else if (ecrecinfan(i1, n2.in1)) {
			return true;
		} else if (ecrecinfan(i1, n2.in2)) {
			return true;
		} else if (ecrecinfan(i1, n2.in3)) {
			return true;
		} else {
			return false;
		}
	}

	bool xmg::ecrecinfan(nodeid i1, nodeid i2) {
		const auto& n2 = _nodes[i2];
		auto nextid = n2.ecrep;
		assert(nextid == i2);
		do {
			const auto nextnode = _nodes[nextid];
			bool res;
			if (is_xor(nextnode)) {
				res = xorrecinfan(i1, nextid);
			} else {
				res = majrecinfan(i1, nextid);
			}
			if (res) {
				return true;
			}
			const auto& nnode = _nodes[nextid];
			nextid = nnode.ecnext;
		} while (nextid != EC_NULL);
		return false;
	}

	void xmg::resetflag() {
		for (auto& n : _nodes) {
			reset_flag(n);
		}
	}

	// Checks if n1 is in the transitive fanin of n2
	bool xmg::infan(nodeid i1, nodeid i2) {
		bool res;
		const auto& n2 = _nodes[i2];
		if (n2.ecrep == i2) { 
			res = ecrecinfan(i1, i2);
		} else {
			if (is_xor(n2)) {
				res = xorrecinfan(i1, i2);
			} else {
				res = majrecinfan(i1, i2);
			}
		}
		resetflag();
		return res;
	}

	void xmg::set_choice(nodeid i1, nodeid i2, bool c, fanoutmap& f) {
		if (i1 == i2) {
			// r is already the choice node for n
			// this can happen due to structural hashing
			return;
		}
		bool addtoeqclass = !infan(i1, i2);
		if (addtoeqclass) {
			const auto& r = _nodes[i2];
			auto nodeid = r.ecrep;
			while (true) {
				auto& node = _nodes[nodeid];
				if (node.ecnext == EC_NULL) {
					node.ecnext = i1;
					break;
				}
				nodeid = node.ecnext;
			}
		}

		auto& i1fanout = f.at(i1);
		if (i1fanout.size() > 0) {
			for (auto nodeid : i1fanout) {
				auto& outnode = _nodes[nodeid];
				f[i2].push_back(nodeid);
				if (outnode.in1 == i1) {
					outnode.in1 = i2;
					if (c) {
						if (is_c1(outnode)) {
							reset_c1(outnode);
						} else {
							set_c1(outnode);
						}
					}
				} else if (outnode.in2 == i1) {
					outnode.in2 = i2;
					if (c) {
						if (is_c2(outnode)) {
							reset_c2(outnode);
						} else {
							set_c2(outnode);
						}
					}
				} else if (outnode.in3 == i1) {
					outnode.in3 = i2;
					if (c) {
						if (is_c3(outnode)) {
							reset_c3(outnode);
						} else {
							set_c3(outnode);
						}
					}
				} else {
					assert(false);
				}
			}
			i1fanout.clear();
		}

		auto& n = _nodes[i1];
		n.ecrep = i2;
		if (c) {
			set_c(n);
		}
	}

	xmg::xmg(const xmg& xmg) {
		_nodes = xmg._nodes;
		_outputs = xmg._outputs;
		_outcompl = xmg._outcompl;
		_innames = xmg._innames;
		_outnames = xmg._outnames;
	}
	
	xmg::xmg(xmg&& xmg) {
		_nodes = std::move(xmg._nodes);
		_outputs = std::move(xmg._outputs);
		_outcompl = std::move(xmg._outcompl);
		_innames = std::move(xmg._innames);
		_outnames = std::move(xmg._outnames);
	}

	xmg& xmg::operator=(xmg&& xmg) {
		_nodes = std::move(xmg._nodes);
		_outputs = std::move(xmg._outputs);
		_outcompl = std::move(xmg._outcompl);
		_innames = std::move(xmg._innames);
		_outnames = std::move(xmg._outnames);
		return *this;
	}

    xmg::xmg(MIG* mig) {
		unordered_map<MAJ3*, pair<nodeid, bool>> nodemap;

		auto torder = mig_topsort(mig);
		for (auto node : torder) {
			if (node->PI || node == mig->one) {
				nodemap[node] = make_pair(create_input(), false);
			} else {
				const auto& p1 = nodemap[node->in1];
				const auto& p2 = nodemap[node->in2];
				const auto& p3 = nodemap[node->in3];
				nodemap[node] = create(
					p1.first, p1.second != node->compl1,
					p2.first, p2.second != node->compl2,
					p3.first, p3.second != node->compl3);
			}
        }
        for (auto i = 0u; i < mig->Nin; i++) {
			_innames.push_back(string(mig->innames[i]));
		}

		for (auto i = 0u; i < mig->Nout; i++) {
			const auto& np = nodemap[mig->out[i]];
			const auto nodeid = np.first; const auto c = np.second;
			const auto migc = static_cast<bool>(mig->outcompl[i]);
			_outputs.push_back(nodeid);
			_outcompl.push_back(c != migc);
			auto& outnode = _nodes[nodeid];
			set_po(outnode);
			_outnames.push_back(string(mig->outnames[i]));
		}
    }

	xmg::xmg(MIG* mig, const xmg_params* p) {
		unordered_map<MAJ3*,pair<nodeid,bool>> nodemap;

		Minisat::Solver solver;
		// A map of nodes to their SAT variables
		varmap var_map;
		// A map from MAJ3 nodes to their corresponding bitvectors
		bvmap bv_map;
		// A map from nodes to their simulation representatives.
		simmap simrep;
		// A map of bitvector hashes to MAJ3 nodes.
		hashmap hashnode_map;

		fanoutmap fanout(mig->Nin+mig->Nnodes+1);

		xmg_stats stats {
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		strashmap shmap(mig->Nnodes/2, stats);

		// Perform random simulation on MIG, computing
		//simrep[n] for every n in the MIG
		auto torder = mig_topsort(mig);
		vector<MAJ3*> inputs(mig->Nin);
		for (auto i = 0u; i < mig->Nin; i++) {
			inputs[i] = mig->in[i];
		}
		initsim(bv_map, hashnode_map, inputs, mig->one);
		simulate(bv_map, hashnode_map, simrep, torder);

		// Create the "one" input
		nodemap[mig->one] = make_pair(create_input(var_map, solver), false);
		for (auto node : torder) {
			if (node == mig->one) {
				continue;
			} else if (node->PI) {
				nodemap[node] = make_pair(create_input(var_map, solver), false);
				continue;
			}
			const auto& p1 = nodemap[node->in1];
			const auto& p1node = _nodes[p1.first];
			const auto& p2 = nodemap[node->in2];
			const auto& p2node = _nodes[p2.first];
			const auto& p3 = nodemap[node->in3];
			const auto& p3node = _nodes[p3.first];
			auto np = find_or_create(
					p1node.ecrep, (is_c(p1node) != p1.second) != node->compl1,
					p2node.ecrep, (is_c(p2node) != p2.second) != node->compl2, 
					p3node.ecrep, (is_c(p3node) != p3.second) != node->compl3, 
					shmap, var_map, solver, fanout);
			if (is_pi(_nodes[np.first])) { 
				// A majority propagation took place
				nodemap[node] = np;
				continue;
			}
			auto rep = simrep[node];
			if (node != rep)  {
				bool c = false;
				const auto rp = nodemap[rep];
				const auto& rpnode = _nodes[rp.first];
				if (np.first != rp.first) {
					++stats.nr_potentials;
					const auto& bv_n = bv_map[node];
					const auto& bv_r = bv_map[rep];
					if (hashbv(bv_n.get()) != hashbv(bv_r.get())) {
						c = true;
					}
					// Check for functional equivalence
					auto equiv = equivalent(np.first, rp.first, c, 
							p->nr_backtracks, var_map, solver);
					if (equiv == l_True) {
						++stats.nr_matches;
						set_choice(np.first, rpnode.ecrep, 
								(c != is_c(rpnode)), fanout);
					} else if (equiv == l_False) {
						// Resimulate MIG with counter example
						++stats.nr_misses;
						hashnode_map.clear();
						inithashmap(bv_map, hashnode_map, inputs, mig->one);
						auto counter = 
							counterexample(mig->Nin, var_map, solver);
						appendcounter(bv_map, counter, inputs, mig->one);
						simulate(bv_map, hashnode_map, simrep, torder);
					} else {
						//hashnode_map.clear();
						++stats.nr_undefined;
					}
				}
			}
			nodemap[node] = np;
		}

		for (auto i = 0u; i < mig->Nin; i++) {
			_innames.push_back(string(mig->innames[i]));
		}

		for (auto i = 0u; i < mig->Nout; i++) {
			const auto& np = nodemap[mig->out[i]];
			const auto nodeid = np.first; const auto c = np.second;
			const auto& node = _nodes[nodeid];
			const auto noderep = node.ecrep; const auto nodec = is_c(node);
			const auto migc = static_cast<bool>(mig->outcompl[i]);
			_outputs.push_back(noderep);
			_outcompl.push_back((c != migc) != nodec);
			auto& outnode = _nodes[noderep];
			set_po(outnode);
			_outnames.push_back(string(mig->outnames[i]));
		}

		// Deallocate memory if possible
		_nodes.resize(_nodes.size());

		cout << "Nr. inputs: " << nin() << endl;
		cout << "Nr. outputs: " << nout() << endl;
		cout << "Nr. nodes: " << nnodes() << endl;
		cout << "Nr. strash hits: " << stats.strash_hits << endl;
		cout << "Nr. potentials: " << stats.nr_potentials << endl;
		cout << "Nr. matches: " << stats.nr_matches << endl;
		cout << "Nr. misses: " << stats.nr_misses << endl;
		cout << "Nr. undefined: " << stats.nr_undefined << endl;
	}

	xmg::xmg(const xmg& xmg, const xmg_params* p) {
		unordered_map<nodeid,pair<nodeid,bool>> nodemap;

		Minisat::Solver solver;
		// A map of nodes to their SAT variables
		varmap var_map;
		// A map from MAJ3 nodes to their corresponding bitvectors
		xbvmap bv_map(xmg.nnodes(), nullptr);
		// A map from nodes to their simulation representatives.
		xsimmap simrep(xmg.nnodes());
		// A map of bitvector hashes to MAJ3 nodes.
		xhashmap hashnode_map;

		fanoutmap fanout(xmg.nnodes());

		xmg_stats stats {
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		strashmap shmap(xmg.nnodes()/2, stats);

		// Perform random simulation on XMG, computing
		//simrep[n] for every n in the XMG
		const auto& nodes = xmg.nodes();
		vector<nodeid> inputs(xmg.nin());
		for (auto i = 0u; i < xmg.nin(); i++) {
			inputs[i] = i+1;
		}
		initsim(bv_map, hashnode_map, inputs);
		simulate(bv_map, hashnode_map, simrep, xmg);

		const auto& innames = xmg.innames();

		for (auto i = 0u; i < xmg.nnodes(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				if (i == 0u) {
					nodemap[node.ecrep] = make_pair(create_input(var_map, solver), false);
				} else {
					nodemap[node.ecrep] = make_pair(create_input(innames[i - 1], var_map, solver), false);
				}
				continue;
			}
			const auto& p1 = nodemap[node.in1];
			const auto& p1node = _nodes[p1.first];
			const auto& p2 = nodemap[node.in2];
			const auto& p2node = _nodes[p2.first];
			pair<nodeid,bool> np;
			if (!is_xor(node)) {
				const auto& p3 = nodemap[node.in3];
				const auto& p3node = _nodes[p3.first];
				np = find_or_create(
						p1node.ecrep, 
						(is_c(p1node) != p1.second) != is_c1(node),
						p2node.ecrep, 
						(is_c(p2node) != p2.second) != is_c2(node),
						p3node.ecrep, 
						(is_c(p3node) != p3.second) != is_c3(node), 
						shmap, var_map, solver, fanout);
			} else {
				np = find_or_create(
						p1node.ecrep, 
						(is_c(p1node) != p1.second) != is_c1(node),
						p2node.ecrep, 
						(is_c(p2node) != p2.second) != is_c2(node),
						shmap, var_map, solver, fanout);
			}
			if (is_pi(_nodes[np.first])) { 
				// A majority propagation took place
				nodemap[i] = np;
				continue;
			}
			auto rep = simrep[i];
			if (static_cast<nodeid>(i) != rep)  {
				bool c = false;
				const auto rp = nodemap[rep];
				const auto& rpnode = _nodes[rp.first];
				if (np.first != rp.first) {
					++stats.nr_potentials;
					const auto& bv_n = bv_map[i];
					const auto& bv_r = bv_map[rep];
					if (hashbv(bv_n.get()) != hashbv(bv_r.get())) {
						c = true;
					}
					// Check for functional equivalence
					auto equiv = equivalent(np.first, rp.first, c, 
							p->nr_backtracks, var_map, solver);
					if (equiv == l_True) {
						++stats.nr_matches;
						set_choice(np.first, rpnode.ecrep, 
								(c != is_c(rpnode)), fanout);
					} else if (equiv == l_False) {
						// Resimulate xmg with counter example
						++stats.nr_misses;
						hashnode_map.clear();
						inithashmap(bv_map, hashnode_map, inputs);
						auto counter = 
							counterexample(xmg.nin(), var_map, solver);
						appendcounter(bv_map, counter, inputs);
						simulate(bv_map, hashnode_map, simrep, xmg);
					} else {
						//hashnode_map.clear();
						++stats.nr_undefined;
					}
				}
			}
			nodemap[i] = np;
		}

		for (auto i = 0u; i < xmg.nin(); i++) {
			_innames.push_back(xmg._innames[i]);
		}

		for (auto i = 0u; i < xmg.nout(); i++) {
			const auto& np = nodemap[xmg._outputs[i]];
			const auto nodeid = np.first; const auto c = np.second;
			const auto& node = _nodes[nodeid];
			const auto noderep = node.ecrep; const auto nodec = is_c(node);
			const auto xmgc = xmg._outcompl[i];
			_outputs.push_back(noderep);
			_outcompl.push_back((c != xmgc) != nodec);
			auto& outnode = _nodes[noderep];
			set_po(outnode);
			_outnames.push_back(xmg._outnames[i]);
		}

		// Deallocate memory if possible
		_nodes.resize(_nodes.size());

		cout << "Nr. inputs: " << nin() << endl;
		cout << "Nr. outputs: " << nout() << endl;
		cout << "Nr. nodes: " << nnodes() << endl;
		cout << "Nr. strash hits: " << stats.strash_hits << endl;
		cout << "Nr. potentials: " << stats.nr_potentials << endl;
		cout << "Nr. matches: " << stats.nr_matches << endl;
		cout << "Nr. misses: " << stats.nr_misses << endl;
		cout << "Nr. undefined: " << stats.nr_undefined << endl;
	}

	xmg strash(const xmg& sxmg) {
		xmg res;

		nodemap nodemap;
		const auto& nodes = sxmg.nodes();
		const auto nnodes = sxmg.nnodes();

		xmg_stats stats{
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		strashmap shmap(nnodes / 2, stats);
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				if (i == 0u) {
					nodemap[i] = make_pair(res.create_input(), is_c);
				} else {
					nodemap[i] = make_pair(res.create_input(), is_c);
				}
			} else {
				auto in1 = nodemap[node.in1];
				in1.second = (in1.second != is_c1(node));
				auto in2 = nodemap[node.in2];
				in2.second = (in2.second != is_c2(node));
				auto in3 = nodemap[node.in3];
				in3.second = (in3.second != is_c3(node));
				nodemap[i] = res.find_or_create(in1, in2, in3, shmap);
			}
		}

		const auto& outputs = sxmg.outputs();
		const auto& outcompl = sxmg.outcompl();
		const auto nouts = outputs.size();
		for (auto i = 0u; i < nouts; i++) {
			auto outid = outputs[i];
			auto outc = outcompl[i];
			auto outnode = nodemap[outid];
			res.create_output(outnode.first, outnode.second != outc);
		}

        const auto& innames = sxmg.innames();
        for (const auto& name : innames) {
            res.add_inname(name);
        }
        
		const auto& outnames = sxmg.outnames();
        for (const auto& name : outnames) {
            res.add_outname(name);
        }

		return res;
	}

	xmg strash_no_compl(const xmg& sxmg) {
		xmg res;

		nodemap nodemap;
		const auto& nodes = sxmg.nodes();
		const auto nnodes = sxmg.nnodes();

		xmg_stats stats{
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		strashmap shmap(nnodes / 2, stats);
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				if (i == 0u) {
					nodemap[i] = make_pair(res.create_input(), is_c);
				} else {
					nodemap[i] = make_pair(res.create_input(), is_c);
				}
			} else {
				auto in1 = nodemap[node.in1];
				in1.second = (in1.second != is_c1(node));
				auto in2 = nodemap[node.in2];
				in2.second = (in2.second != is_c2(node));
				auto in3 = nodemap[node.in3];
				in3.second = (in3.second != is_c3(node));
				nodemap[i] = res.find_or_create_no_compl(in1, in2, in3, shmap);
			}
		}

		const auto& outputs = sxmg.outputs();
		const auto& outcompl = sxmg.outcompl();
		const auto nouts = outputs.size();
		for (auto i = 0u; i < nouts; i++) {
			auto outid = outputs[i];
			auto outc = outcompl[i];
			auto outnode = nodemap[outid];
			res.create_output(outnode.first, outnode.second != outc);
		}

        const auto& innames = sxmg.innames();
        for (const auto& name : innames) {
            res.add_inname(name);
        }
        
		const auto& outnames = sxmg.outnames();
        for (const auto& name : outnames) {
            res.add_outname(name);
        }

		return res;
	}

	xmg rdup(const xmg& sxmg) {
		xmg res;

		nodemap nodemap;
		const auto& nodes = sxmg.nodes();
		const auto nnodes = sxmg.nnodes();

		xmg_stats stats{
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		strashmap shmap(nnodes / 2, stats);
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res.create_input(), is_c);
			} else {
				auto in1 = nodemap[node.in1];
				in1.second = (in1.second != is_c1(node));
				auto in2 = nodemap[node.in2];
				in2.second = (in2.second != is_c2(node));
				auto in3 = nodemap[node.in3];
				in3.second = (in3.second != is_c3(node));
				nodemap[i] = res.find_or_create_no_prop(in1, in2, in3, shmap);
			}
		}

		const auto& outputs = sxmg.outputs();
		const auto& outcompl = sxmg.outcompl();
		const auto nouts = outputs.size();
		for (auto i = 0u; i < nouts; i++) {
			auto outid = outputs[i];
			auto outc = outcompl[i];
			auto outnode = nodemap[outid];
			res.create_output(outnode.first, outnode.second != outc);
		}
        
		const auto& innames = sxmg.innames();
        for (const auto& name : innames) {
            res.add_inname(name);
        }
        
		const auto& outnames = sxmg.outnames();
        for (const auto& name : outnames) {
            res.add_outname(name);
        }

		return res;
	}

	void xmg::create_output(nodeid nodeid, bool c) {
		auto& outnode = _nodes[nodeid];
		set_po(outnode);
		_outputs.push_back(nodeid);
		_outcompl.push_back(c);
    }

    void xmg::create_output(nodeid nodeid, bool c, const string& name) {
        create_output(nodeid, c);
		_outnames.push_back(name);
	}


	unique_ptr<xmg_params> default_xmg_params() {
		unique_ptr<xmg_params> res(new xmg_params);
		res->nr_backtracks = DEFAULT_BACKTRACKS;
		return res;
	}

	bool xmg::is_mig() const {
		for (const auto& node : _nodes) {
			if (!(is_pi(node) || is_maj(node))) {
				return false;
			}
		}
		return true;
	}

	MIG* xmg::extractmig() const {
		MIG* mig = (MIG*)malloc(sizeof(MIG));
		mig->Nin = nin();
		mig->Nout = nout();
		mig->in = (MAJ3**)malloc(sizeof(MAJ3*)*mig->Nin);
		mig->innames = (char**)malloc(sizeof(char*)*mig->Nin);
		mig->out = (MAJ3**)malloc(sizeof(MAJ3*)*mig->Nout);
		mig->outnames = (char**)malloc(sizeof(char*)*mig->Nout);
		mig->outcompl = (unsigned*)malloc(sizeof(unsigned)*mig->Nout);
		mig->Nnodes = 0;
		mig->nodes = NULL;

		vector<MAJ3*> nodes;
		unordered_map<nodeid,MAJ3*> nodemap;
		for (auto i = 0u; i < _nodes.size(); i++) {
			const auto& node = _nodes[i];
			if (is_pi(node)) {
				auto n = (MAJ3*)malloc(sizeof(MAJ3));
				assert(n != NULL);
				n->in1=n->in2=n->in3=NULL;
				n->outEdges=NULL;
				n->value=2;
				n->aux=NULL;
				n->fanout=0;
				n->label=i;
				n->flag=0;
				n->compl1=n->compl2=n->compl3=0;
				n->PI=1;
				n->PO=0;
				n->level = 0;
				if (i == 0) { // One node
					mig->one = n;
				} else {
					--n->label;
					mig->innames[i-1] = (char*)malloc(sizeof(char)*MAXCAR);
					assert(mig->innames[i - 1] != NULL);
					strcpy(mig->innames[i-1], _innames[i-1].c_str());
				}
				nodemap[node.ecrep] = n;
			} else if (node.ecrep != static_cast<nodeid>(i)) {
				continue;
			} else {
				auto n = new_node_nop(mig, nodemap[node.in1], is_c1(node),
						nodemap[node.in2], is_c2(node), 
						nodemap[node.in3], is_c3(node));
				nodemap[node.ecrep] = n;
			}
		}

		for (auto i = 0u; i < mig->Nout; i++) {
			mig->out[i] = nodemap[_outputs[i]];
			mig->outcompl[i] = _outcompl[i];
			mig->outnames[i] = (char*)malloc(sizeof(char)*MAXCAR);
			assert(mig->outnames[i] != NULL);
			strcpy(mig->outnames[i], _outnames[i].c_str());
		}

		return mig;
	}

	bool xmg::equals(const xmg& other) const {
		if (other.nin() != nin() || other.nout() != nout()) {
			return false;
		}
		Solver solver;
		varmap vm1, vm2;
		const auto nin = other.nin();
		Lit onelit;
		for (auto i = 0u; i <= nin; i++) {
			const auto& node1 = nodes()[i];
			const auto& node2 = other.nodes()[i];
			auto newvar = solver.newVar();
			auto newlit = mkLit(newvar, false);
			if (i == 0) {
				onelit = newlit;
			}
			auto complvar = solver.newVar();
			auto compllit = mkLit(complvar, false);
			solver.addClause(~newlit, ~compllit);
			solver.addClause(newlit, compllit);
			if (is_pi_c(node1)) {
				vm1[i] = complvar;
			} else {
				vm1[i] = newvar;
			}
			if (is_pi_c(node2)) {
				vm2[i] = complvar;
			} else {
				vm2[i] = newvar;
			}
		}

		for (auto i = 0u; i < nnodes(); i++) {
			const auto& node = nodes()[i];
			if (is_pi(node)) {
				continue;
			}
			if (is_xor(node)) {
				add_clauses(node.in1, is_c1(node), node.in2, is_c2(node), i, vm1, solver);
			} else {
				add_clauses(node.in1, is_c1(node), node.in2, is_c2(node), node.in3, is_c3(node), i, vm1, solver);
			}
		}

		for (auto i = 0u; i < other.nnodes(); i++) {
			const auto& node = other.nodes()[i];
			if (is_pi(node)) {
				continue;
			}
			if (is_xor(node)) {
				add_clauses(node.in1, is_c1(node), node.in2, is_c2(node), i, vm2, solver);
			} else {
				add_clauses(node.in1, is_c1(node), node.in2, is_c2(node), node.in3, is_c3(node), i, vm2, solver);
			}
		}

		// XOR the ouputs 
		varmap outvar_map;
		for (auto i = 0u; i < nout(); i++) {
			auto var1 = vm1[_outputs[i]];
			auto c1 = _outcompl[i];
			auto var2 = vm2[other._outputs[i]];
			auto c2 = other._outcompl[i];

			auto new_var = solver.newVar();
			outvar_map[i] = new_var;
			auto nlit = mkLit(new_var, false);
			auto lit1 = mkLit(var1, c1);
			auto lit2 = mkLit(var2, c2);

			// Generate XOR clauses
			solver.addClause(~lit1, ~lit2, ~nlit);
			solver.addClause(lit1, lit2, ~nlit);
			solver.addClause(lit1, ~lit2, nlit);
			solver.addClause(~lit1, lit2, nlit);
		}

		// OR all the outputs together
		Lit miterlit;
		varmap orvar_map;
		if (nout() > 1) {
			for (auto i = 1u; i < nout(); i++) {
				auto or_var = solver.newVar();
				auto nlit = mkLit(or_var, false);
				orvar_map[i] = or_var;

				auto lit1 = mkLit(outvar_map[i - 1], false);
				auto lit2 = mkLit(outvar_map[i], false);

				solver.addClause(~lit1, ~lit2, nlit);
				solver.addClause(~lit1, ~onelit, nlit);
				solver.addClause(lit1, lit2, ~nlit);
				solver.addClause(lit1, onelit, ~nlit);
				solver.addClause(~lit2, ~onelit, nlit);
				solver.addClause(lit2, onelit, ~nlit);

				if (i == (nout() - 1)) { // The final output variable will be the miter output
					miterlit = mkLit(or_var, false);
				}
			}
		} else {
			assert(nout() == 1);
			miterlit = mkLit(outvar_map[0], false);
		}

		// Solve the miter
		vec<Lit> assumps;
		assumps.push(miterlit);
		assumps.push(onelit);
		if (solver.solve(assumps)) {
			return false;
		} else {
			return true;
		}
	}

	bool simulate_node(const node& n, unordered_map<nodeid, bool>& simval) {
		if (is_xor(n)) {
			auto in1val = simval[n.in1];
			auto truein1 = (is_c1(n) ? !in1val : in1val);
			auto in2val = simval[n.in2];
			auto truein2 = (is_c2(n) ? !in2val : in2val);
			return truein1 != truein2;
		} else {
			auto in1val = simval[n.in1];
			auto truein1 = (is_c1(n) ? !in1val : in1val);
			auto in2val = simval[n.in2];
			auto truein2 = (is_c2(n) ? !in2val : in2val);
			auto in3val = simval[n.in3];
			auto truein3 = (is_c3(n) ? !in3val : in3val);
			return (truein1 && truein2) || (truein1 && truein3) || (truein2 && truein3);
		}
	}

	boost::dynamic_bitset<> simulate_xmg(const xmg& xmg) {
		assert(xmg.nin() <= 6);
		assert(xmg.nout() == 1);
		boost::dynamic_bitset<> func;

		const auto nin = xmg.nin();
		const auto& nodes = xmg.nodes();
		const auto& nnodes = xmg.nnodes();
		const auto& outputs = xmg.outputs();
		const auto& outcompl = xmg.outcompl();
		const auto nsimvectors = (1u << nin);
		unordered_map<nodeid, bool> simval;
		for (auto j = 0u; j < nsimvectors; j++) {
			boost::dynamic_bitset<> invec(nin, j);
			for (auto i = 0u; i < nnodes; i++) {
				const auto& node = nodes[i];
				if (is_pi(node)) {
					auto is_complemented = is_pi_c(node);
					if (i == 0u) { // Const 1
						simval[i] = (is_complemented != true);
					} else {
						simval[i] = (is_complemented != invec.test(i - 1));
					}
				} else {
					simval[i] = simulate_node(node, simval);
				}
			}
			auto outval = outcompl[0] ? !simval[outputs[0]] : simval[outputs[0]];
			func.push_back(outval);
		}
		assert(func.size() == nsimvectors);

		return func;
	}

    void xmg::create_dummy_innames() {
        _innames.clear();
		auto count = 0u;
		for (const auto& node : _nodes) {
			if (!is_pi(node)) {
				break;
			}
			_innames.push_back("x" + std::to_string(count));
			++count;
		}
	}

	void xmg::create_dummy_outnames() {
        _outnames.clear();
		for (auto i = 0u; i < _outputs.size(); i++) {
			_outnames.push_back("f[" + std::to_string(i) + "]");
		}
	}

	void xmg::create_dummy_names() {
		create_dummy_innames();
		create_dummy_outnames();
	}

	string xmg::to_verilog() const {
		stringstream s;
		write_verilog(*this, s);
		return s.str();
	}
}

