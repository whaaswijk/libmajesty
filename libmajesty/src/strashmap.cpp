#include "strashmap.h"
#include <cassert>
#include <iostream>
#include "mlpext.h"

/* function pair
Synopsis [NxN->N pairing function]

Description [Cantor pairing function to uniquely map 2 integer into another
integer]

Side effects [-]
*/ 
static inline unsigned cantor_pair(unsigned i, unsigned j) { 
	unsigned int p;

	p=((i+j)*(i+j+1)/2)+i;

	return p;
}


#define SWAP(x, xc, y, yc) if (y < x) {\
	tn = x; tc = xc; x = y; xc = yc; y = tn; yc = tc; }

static inline void sort_inputs(
		int32_t& in1, bool& c1, 
		int32_t& in2, bool& c2,
		int32_t& in3, bool& c3) {
		int32_t tn; bool tc;
		SWAP(in2, c2, in3, c3);
		SWAP(in1, c1, in3, c3);
		SWAP(in1, c1, in2, c2);
}

static inline void 
sort_inputs(int32_t& in1, bool& c1, int32_t& in2, bool& c2) {
		int32_t tn; bool tc;
		SWAP(in1, c1, in2, c2);
}

namespace majesty {

	unsigned strashmap::hash(maj3inputs) const {
			auto bigprime = 15485863u;// 1000000th prime

			auto xy = cantor_pair(in1+c1, in2+c2) % bigprime;
			auto h= cantor_pair(xy, in3+c3) % bigprime;

			return h % _size;
	}

	unsigned strashmap::hash(xorinputs) const {
			auto bigprime = 15485863u;// 1000000th prime

			auto xy = cantor_pair(in1+c1, in2+c2) % bigprime;

			return xy % _size;
	}

	strashmap::strashmap(size_t size, xmg_stats& stats) : _stats(stats) {
		_size = size;
		_table.resize(size);
	}

	int32_t strashmap::find(maj3inputs, xmg& m) const {
		sort_inputs(in1, c1, in2, c2, in3, c3);
		auto key = hash(in1, c1, in2, c2, in3, c3);
		const auto& entries = _table[key];
		for (auto idx : entries) {
			const auto& n = m._nodes[idx];
			if (is_xor(n)) {
				continue;
			}
			auto nc1 = is_c1(n);
			auto nc2 = is_c2(n);
			auto nc3 = is_c3(n);
			if (n.in1 == in1 && nc1 == c1 && n.in2 == in2 &&
					nc2 == c2 && n.in3 == in3 && nc3 == c3) {
				++_stats.strash_hits;
				return idx;
			}
		}
		return -1;
	}

	int32_t strashmap::find(xorinputs, xmg& m) const {
		sort_inputs(in1, c1, in2, c2);
		auto key = hash(in1, c1, in2, c2);
		const auto& entries = _table[key];
		for (auto idx : entries) {
			const auto& n = m._nodes[idx];
			if (!is_xor(n)) {
				continue;
			}
			auto nc1 = is_c1(n);
			auto nc2 = is_c2(n);
			if (n.in1 == in1 && nc1 == c1 && n.in2 == in2 && nc2 == c2) {
				++_stats.strash_hits;
				return idx;
			}
		}
		return -1;
	}
	
	int32_t strashmap::find_or_add(maj3inputs, xmg& m) {
		auto found = find(in1, c1, in2, c2, in3, c3, m);
		if (found != -1) {
			return found;
		}
		sort_inputs(in1, c1, in2, c2, in3, c3);
		auto idx = m.create_node(in1, c1, in2, c2, in3, c3);
		auto key = hash(in1, c1, in2, c2, in3, c3);
		auto& entries = _table[key];
		entries.push_back(idx);
		return idx;
	}

	int32_t strashmap::find_or_add(maj3inputs, xmg& m, 
			varmap& v, Minisat::Solver& s, fanoutmap& f) {
		auto found = find(in1, c1, in2, c2, in3, c3, m);
		if (found != -1) {
			return found;
		}
		sort_inputs(in1, c1, in2, c2, in3, c3);
		auto idx = m.create_node(in1, c1, in2, c2, in3, c3, v, s, f);
		auto key = hash(in1, c1, in2, c2, in3, c3);
		auto& entries = _table[key];
		entries.push_back(idx);
		return idx;
	}

	int32_t strashmap::find_or_add(xorinputs, xmg& m) {
		auto found = find(in1, c1, in2, c2, m);
		if (found != -1) {
			return found;
		}
		sort_inputs(in1, c1, in2, c2);
		auto idx = m.create_node(in1, c1, in2, c2);
		auto key = hash(in1, c1, in2, c2);
		auto& entries = _table[key];
		entries.push_back(idx);
		return idx;
	}

	int32_t strashmap::find_or_add(xorinputs, xmg& m, 
			varmap& v, Minisat::Solver& s, fanoutmap& f) {
		auto found = find(in1, c1, in2, c2, m);
		if (found != -1) {
			return found;
		}
		sort_inputs(in1, c1, in2, c2);
		auto idx = m.create_node(in1, c1, in2, c2, v, s, f);
		auto key = hash(in1, c1, in2, c2);
		auto& entries = _table[key];
		entries.push_back(idx);
		return idx;
	}
}

