#include <memory>
#include <algorithm>
#include <bitset>
#include <math.h>
#include "cut.h"

using namespace std;
using namespace cirkit;

inline unsigned hamming_weight(unsigned v) {
#ifdef WIN32
	return __popcnt(v);
#else
	return __builtin_popcount(v);
#endif
}

namespace majesty {

	cutvec eqclass_cuts(const vector<node>& nodes, const node& N, 
			const cutmap& cut_map, const cut_params* p) {
		// Compute the union of the k-feasible cuts of this equivalence class,
		// while filtering out duplicates and dominated cuts.
		cutvec res;
		auto nextidx = N.ecrep;
		do {
			const auto& node = nodes[nextidx];
			auto cuts = node_cuts(node, cut_map, p);
			for (auto& cut : cuts) {
				// Check equivalence and redundancy
				safeinsert(res, move(cut));
			}
			nextidx = node.ecnext;
		} while (nextidx != EC_NULL);

		// The trivial cut of the functional representative
		unique_ptr<cut> c(new cut(N.ecrep));
		res.push_back(move(c));

		return res;
	}

	cutvec filtered_eqclass_cuts(const vector<node>& nodes, const node& N,
			const cutmap& cut_map, const cut_params* p, funcmap& fm, const unordered_set<unsigned long>& timeoutfuncs) {
		// Compute the union of the k-feasible cuts of this equivalence class,
		// while filtering out duplicates and dominated cuts.
		cutvec res;
		auto nextidx = N.ecrep;
		do {
			const auto& node = nodes[nextidx];
			auto cuts = node_cuts(node, cut_map, p);
			for (auto& cut : cuts) {
				cut->computefunction(node, fm);
				auto& cutfunc = fm[cut.get()];
				if (timeoutfuncs.find(cutfunc->to_ulong()) == timeoutfuncs.end()) { // Make sure it's not a function we've timed out on
					// Check equivalence and redundancy
					safeinsert(res, move(cut));
				}
			}
			nextidx = node.ecnext;
		} while (nextidx != EC_NULL);

		// The trivial cut of the functional representative
		unique_ptr<cut> c(new cut(N.ecrep));
		res.push_back(move(c));
		c->computefunction(N, fm);

		return res;
	}

	cutvec node_cuts(const node& n, 
			const cutmap& cut_map, const cut_params *p) {
		const auto& cuts1 = cut_map.at(n.in1);
		const auto& cuts2 = cut_map.at(n.in2);
		if (is_xor(n)) {
			return crossproduct(cuts1, cuts2, p);
		} else {
			const auto& cuts3 = cut_map.at(n.in3);
			return crossproduct(cuts1, cuts2, cuts3, p);
		}
	}

	cutvec crossproduct(
			const cutvec& cuts1, 
			const cutvec& cuts2, 
			const cut_params *p) {
		cutvec res;

		for (const auto& cut1 : cuts1) {
			for (const auto& cut2 : cuts2) {
				auto sig1 = cut1->sig();
				auto sig2 = cut2->sig();
				if (hamming_weight(sig1 | sig2) > p->klut_size) {
					// Merged cut size will be > k
					continue;
				}
				if (cut1->size() == p->klut_size && 
						cut2->size() == p->klut_size) {
					if (sig1 != sig2) {
						// Merged cut cannot be k-feasible
						continue;
					}

					// Cut can be added iff cuts are equal
					if (!cut1->equal(cut2.get())) {
						continue;
					}
				}
				// Create the merged cut and check it's size
				unique_ptr<cut> merged(new cut(cut1.get(), cut2.get()));
				if (merged->size() <= p->klut_size) {
					safeinsert(res, move(merged));
				}
			}
		}

		return res;
	}

	cutvec crossproduct(
			const cutvec& cuts1,
			const cutvec& cuts2,
			const cutvec& cuts3,
			const cut_params *p) {
		cutvec res;

		for (const auto& cut1 : cuts1) {
			for (const auto& cut2 : cuts2) {
				for (const auto& cut3 : cuts3) {
					auto sig1 = cut1->sig();
					auto sig2 = cut2->sig();
					auto sig3 = cut3->sig();
					if (hamming_weight(sig1 | sig2 | sig3) > p->klut_size) {
						// Merged cut size will be > k
						continue;
					}
					if (cut1->size() == p->klut_size && 
							cut2->size() == p->klut_size) {
						if (sig1 != sig2) {
							// Merged cut cannot be k-feasible
							continue;
						} else if (!cut1->equal(cut2.get())) {
							continue;
						}
					} else if (cut1->size() == p->klut_size &&
							cut3->size() == p->klut_size) {
						if (sig1 != sig3) {
							continue;
						} else if (!cut1->equal(cut3.get())) {
							continue;
						}
					} else if (cut2->size() == p->klut_size &&
							cut3->size() == p->klut_size) {
						if (sig2 != sig3) {
							continue;
						} else if (!cut2->equal(cut3.get())) {
							continue;
						}
					}
					// Create the merged cut and check it's size
					unique_ptr<cut> merged(
							new cut(cut1.get(), cut2.get(), cut3.get()));
					if (merged->size() <= p->klut_size) {
						safeinsert(res, move(merged));
					}
				}
			}
		}

		return res;
	}

	// Safely inserts a cut into a vector of cuts, avoiding duplicated
	// and redundant cuts.
	void safeinsert(cutvec& v, unique_ptr<cut> c) {
		auto sig1 = c->sig();
		for (auto i = 0u; i < v.size(); i++) {
			auto cp = v[i].get();
			auto sig2 = cp->sig();
			// Check for equality first
			if (sig1 == sig2 && cp->size() == c->size()) { 
				if (c->equal(cp)) {
					return;
				}
			}

			// Check for redundancy
			if ((sig1 | sig2) == sig2) {
				// Requires detailed domination check
				if (c->dominates(cp)) {
					v[i].reset(nullptr);
					continue;
				}
			}
			if ((sig1 | sig2) == sig1) {
				// Requires detailed domination check
				if (cp->dominates(c.get())) {
					return;
				}
			}
		}
		cutvec newv;
		for (auto i = 0u; i < v.size(); i++) {
			if (v[i] != nullptr) {
				newv.push_back(move(v[i]));
			}
		}
		newv.push_back(move(c));
		v.clear();
		for (auto& c : newv) {
			v.push_back(move(c));
		}
	}

	cut::cut() {
		_sig = 0u;
	}

	cut::cut(nodeid n) {
		_nodes.push_back(n);
		_sig = (1u) << (n % CHASH_BUF);
	}

	inline vector<nodeid> 
		mergevectors(const vector<nodeid>& v1, const vector<nodeid>& v2) {
			vector<nodeid> v;
			auto i1 = 0u; auto i2 = 0u;
			auto s1 = v1.size(); auto s2 = v2.size();
			while (i1 < s1 && i2 < s2) {
				auto node1 = v1[i1]; auto node2 = v2[i2];
				if (node1 < node2) {
					v.push_back(node1);
					++i1;
				} else if (node1 > node2) {
					v.push_back(node2);
					++i2;
				} else {
					v.push_back(node1);
					++i1; ++i2;
				}
			}
			while (i1 < s1) {
				v.push_back(v1[i1++]);
			} 
			while (i2 < s2) {
				v.push_back(v2[i2++]);
			}
			return v;
		}

	// Creates a cut (possibly not k-feasible) by merging two existing ones
	cut::cut(cut* c1, cut* c2) {
		_nodes = mergevectors(c1->_nodes, c2->_nodes);
		_c1 = c1; _c2 = c2;
		_sig = c1->_sig | c2->_sig;
	}

	// Creates a cut (possibly not k-feasible) by merging three existing ones
	cut::cut(cut* c1, cut* c2, cut* c3) {
		auto m1 = mergevectors(c1->_nodes, c2->_nodes);
		_nodes = mergevectors(m1, c3->_nodes);
		_c1 = c1; _c2 = c2; _c3 = c3;
		_sig = c1->_sig | c2->_sig | c3->_sig;
	}

	// NOTE: assumes both cuts have the same size!
	bool cut::equal(const cut* c) const {
		for (auto i = 0u; i < _nodes.size(); i++) {
			if (_nodes[i] != c->_nodes[i]) {
				return false;
			}
		}
		return true;
	}

	bool cut::dominates(const cut* c) const {
		if (c->size() > this->size()) {
			for (auto node : _nodes) {
				if (find(c->_nodes.begin(), c->_nodes.end(), node) == 
						c->_nodes.end())
					return false;
			}
			return true;
		} else {
			return false;
		}
	}

	cutmap enumerate_cuts(const xmg& m, const cut_params *p) {
		cout << "Enumerating cuts..." << endl;
		cutmap cut_map(m.nnodes());

		// We assume that the nodes are stored in topological order
		auto nodes = m.nodes();
		auto total_nodes = nodes.size();
		auto processed_nodes = 0u;
		for (auto i = 0u; i < total_nodes; i++) {
			cutvec res;
			const auto& n = nodes[i];
			if (i == 0) { // One node
				unique_ptr<cut> oc(new cut());
				res.push_back(move(oc));
			} else if (is_pi(n)) {
				unique_ptr<cut> c(new cut(n.ecrep));
				res.push_back(move(c));
			} else if (n.ecrep != static_cast<nodeid>(i)) { 
				// Consider only equivalence class representatives
				++processed_nodes;
				continue;
			} else {
				res = eqclass_cuts(nodes, n, cut_map, p);
			}
			cut_map[n.ecrep] = move(res);
			cout << "Progress: (" << ++processed_nodes << "/" << total_nodes;
			cout << ")\r";
		}
		cout << endl;

		return cut_map;
	}

	cutmap filtered_enumerate_cuts(const xmg& m, const cut_params *p, funcmap& fm, const unordered_set<unsigned long>& timeoutfuncs) {
		cout << "Enumerating cuts..." << endl;
		cutmap cut_map(m.nnodes());
		fm.clear();

		// We assume that the nodes are stored in topological order
		auto nodes = m.nodes();
		auto total_nodes = nodes.size();
		auto processed_nodes = 0u;
		for (auto i = 0u; i < total_nodes; i++) {
			cutvec res;
			const auto& n = nodes[i];
			if (i == 0) { // One node
				unique_ptr<cut> oc(new cut());
				res.push_back(move(oc));
				unique_ptr<tt> f(new tt(tt_const1()));
				fm[oc.get()] = std::move(f);
			} else if (is_pi(n)) {
				unique_ptr<cut> c(new cut(n.ecrep));
				res.push_back(move(c));
				unique_ptr<tt> f(new tt(tt_nth_var(0)));
				fm[c.get()] = std::move(f);
			} else if (n.ecrep != static_cast<nodeid>(i)) { 
				// Consider only equivalence class representatives
				++processed_nodes;
				continue;
			} else {
				res = filtered_eqclass_cuts(nodes, n, cut_map, p, fm, timeoutfuncs);
			}
			cut_map[n.ecrep] = move(res);
			cout << "Progress: (" << ++processed_nodes << "/" << total_nodes;
			cout << ")\r";
		}
		cout << endl;

		return cut_map;
	}

	void cut::computefunction(const node& node, funcmap& fm) { 
		if (_nodes.size() == 1) {
			// This is a trivial cut, we cannot compute its function
			// based on its child cuts as it has none
			unique_ptr<tt> f(new tt(tt_nth_var(0)));
			fm[this] = std::move(f);
		} else {
			if (is_xor(node)) {
				compute_xor_function(node, fm);
			} else {
				compute_maj3_function(node, fm);
			}
		}
	}

	void cut::compute_maj3_function(const node& node, funcmap& m) { 
		map<nodeid,unsigned> sigma;
		for (auto i = 0u; i < _nodes.size(); i++) {
			sigma[_nodes[i]] = i;
		}

		// Express parent cut functions in terms of new variable assignment
		auto f1 = *m.at(_c1);
		for (int i = _c1->size()-1; i >= 0; i--) {
			auto s_m = tt_nth_var(sigma[_c1->_nodes[i]]);
			auto plus = tt_cof1(f1, i);
			auto min = tt_cof0(f1, i);
			tt_align(s_m, plus); tt_align(s_m, min);
			f1 = (s_m & plus) | (~s_m & min);
		}
		if (is_c1(node)) {
			f1 = ~f1;
		}
		auto f2 = *m.at(_c2);
		for (int i = _c2->size()-1; i >= 0; i--) {
			auto s_m = tt_nth_var(sigma[_c2->_nodes[i]]);
			auto plus = tt_cof1(f2, i);
			auto min = tt_cof0(f2, i);
			tt_align(s_m, plus); tt_align(s_m, min);
			f2 = (s_m & plus) | (~s_m & min);
		}
		if (is_c2(node)) {
			f2 = ~f2;
		}
		auto f3 = *m.at(_c3);
		for (int i = _c3->size()-1; i >= 0; i--) {
			auto s_m = tt_nth_var(sigma[_c3->_nodes[i]]);
			auto plus = tt_cof1(f3, i);
			auto min = tt_cof0(f3, i);
			tt_align(s_m, plus); tt_align(s_m, min);
			f3 = (s_m & plus) | (~s_m & min);
		}
		if (is_c3(node)) {
			f3 = ~f3;
		}
		tt_align(f1, f2); tt_align(f1, f3); tt_align(f1, f2);

		unique_ptr<tt> function(new tt((f1 & f2) | (f1 & f3) | (f2 & f3)));
		if (is_c(node)) {
			*function = ~(*function);
		}
		tt support;
		tt_to_minbase(*function, &support);
		vector<nodeid> nonvacs;
		for (auto i = 0u; i < _nodes.size(); i++) {
			if (support.test(i))
				nonvacs.push_back(_nodes[i]);
		}
		_nodes = nonvacs;
		computesignature();

		m[this] = std::move(function);
	}

	void cut::compute_xor_function(const node& node, funcmap& m) { 
		map<nodeid,unsigned> sigma;
		for (auto i = 0u; i < _nodes.size(); i++) {
			sigma[_nodes[i]] = i;
		}

		// Express parent cut functions in terms of new variable assignment
		auto f1 = *m.at(_c1);
		for (int i = _c1->size()-1; i >= 0; i--) {
			auto s_m = tt_nth_var(sigma[_c1->_nodes[i]]);
			auto plus = tt_cof1(f1, i);
			auto min = tt_cof0(f1, i);
			tt_align(s_m, plus); tt_align(s_m, min);
			f1 = (s_m & plus) | (~s_m & min);
		}
		if (is_c1(node)) {
			f1 = ~f1;
		}
		auto f2 = *m.at(_c2);
		for (int i = _c2->size()-1; i >= 0; i--) {
			auto s_m = tt_nth_var(sigma[_c2->_nodes[i]]);
			auto plus = tt_cof1(f2, i);
			auto min = tt_cof0(f2, i);
			tt_align(s_m, plus); tt_align(s_m, min);
			f2 = (s_m & plus) | (~s_m & min);
		}
		if (is_c2(node)) {
			f2 = ~f2;
		}
		tt_align(f1, f2); 
		unique_ptr<tt> function(new tt(f1 ^ f2));
		if (is_c(node)) {
			*function = ~(*function);
		}
		tt support;
		tt_to_minbase(*function, &support);
		vector<nodeid> nonvacs;
		for (auto i = 0u; i < _nodes.size(); i++) {
			if (support.test(i))
				nonvacs.push_back(_nodes[i]);
		}
		_nodes = nonvacs;
		computesignature();

		m[this] = std::move(function);
	}

	void cut::computesignature() {
		auto sig = 0u;
		auto one = 1u;
		for (auto nodeid : _nodes) {
			sig |= one << (nodeid % CHASH_BUF);
		}
		_sig = sig;
	}

	unique_ptr<cut_params> default_cut_params() {
		unique_ptr<cut_params> res(new cut_params);
		res->klut_size = DEFAULT_KLUT_SIZE;
		res->num_cuts = DEFAULT_CUT_LIMIT;
		return res;
	}

	void cut::set_required() {
		if (_is_required) {
			return;
		}
		_is_required = true;
		if (_c1 == nullptr) {
			return;
		}
		_c1->set_required();
		_c2->set_required();
		if (_c3 == nullptr) {
			return;
		}
		_c3->set_required();
	}
}
