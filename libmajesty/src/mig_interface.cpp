#include <mig_interface.h>
#include <xmg.h>
#include "strashmap.h"
#include <algorithm>
#include <iostream>
#include <truth_table_utils.hpp>
#include <convert.h>
#include <boost/pending/integer_log2.hpp>

using namespace std;

namespace majesty {
	
	xmg* random_mig(mt19937_64& rng, unsigned ninputs, unsigned nnodes, unsigned noutputs) {
		static uniform_int_distribution<mt19937_64::result_type> rule_dist(0, 9);
		static uniform_int_distribution<mt19937_64::result_type> rule_sel(0, 1);
		static vector<nodeid> n1(3);
		static vector<bool> c1(3);

		auto mig = new xmg();
		// Create the "one" input
		mig->create_input();
		for (auto i = 0u; i < ninputs; i++) {
			mig->create_input();
		}
		while (mig->nnodes() < nnodes) {
			uniform_int_distribution<mt19937_64::result_type> node_dist(0, mig->nnodes() - 1);
			if (rule_dist(rng) == 0) { // Ensure that one of the rules applies
				// For associativity we need a random non-PI node, so check that
				// such a node exists. If not, make random Maj3 node instead.
				if (rule_sel(rng) == 0 || mig->nin() + 1 == mig->nnodes()) { // Maj3
					n1[0] = static_cast<nodeid>(node_dist(rng));
					c1[0] = rule_sel(rng) == 0;
					n1[1] = n1[0];
					c1[1] = rule_sel(rng) == 0;
					n1[2] = static_cast<nodeid>(node_dist(rng));
					c1[2] = rule_sel(rng) == 0;
					random_shuffle(n1.begin(), n1.end());
					random_shuffle(c1.begin(), c1.end());
					mig->create(n1[0], c1[0], n1[1], c1[1], n1[2], c1[2]);
				} else {
					uniform_int_distribution<mt19937_64::result_type>
						nonpi_dist(mig->nin() + 1, mig->nnodes() - 1);
					n1[0] = static_cast<nodeid>(nonpi_dist(rng));
					assert(!is_pi(mig->nodes()[n1[0]]));
					n1[1] = mig->nodes()[n1[0]].in1;
					n1[2] = static_cast<nodeid>(node_dist(rng));
					c1[0] = rule_sel(rng) == 0;
					c1[1] = rule_sel(rng) == 0;
					c1[2] = rule_sel(rng) == 0;
					random_shuffle(n1.begin(), n1.end());
					random_shuffle(c1.begin(), c1.end());
					mig->create(n1[0], c1[0], n1[1], c1[1], n1[2], c1[2]);
				}
			} else { // Make a random node
				mig->create(
					static_cast<nodeid>(node_dist(rng)), rule_sel(rng) == 0,
					static_cast<nodeid>(node_dist(rng)), rule_sel(rng) == 0,
					static_cast<nodeid>(node_dist(rng)), rule_sel(rng) == 0);
			}
		}

		uniform_int_distribution<mt19937_64::result_type> node_dist(0, mig->nnodes() - 1);
		for (auto i = 0u; i < noutputs; i++) {
			auto outidx = static_cast<nodeid>(node_dist(rng));
			auto outcompl = rule_sel(rng) == 0;
			mig->create_output(outidx, outcompl, "out" + i);
		}

		return mig;
	}

	xmg* mig_manager::create_random_graph(unsigned ninputs, unsigned nnodes) {
		return random_mig(rng, ninputs, nnodes, ninputs);
	}
	
	xmg* mig_manager::create_random_graph(unsigned ninputs, unsigned nnodes, unsigned noutputs) {
		return random_mig(rng, ninputs, nnodes, noutputs);
	}

	xmg* mig_manager::random_mig_decomposition(unsigned ninputs) {
		auto num_funcs = (1u << (1u << ninputs));
		static uniform_int_distribution<mt19937_64::result_type> rule_dist(0, num_funcs-1);
		auto rand_func = rule_dist(rng);
		return new xmg(mig_decompose(ninputs, rand_func));
	}

	xmg* mig_string_decompose(const string& func) {
		auto ninputs = boost::integer_log2(func.size());
		return new xmg(mig_decompose(ninputs, func));
	}
	
	xmg* mig_expression_decompose(unsigned ninputs, const string& expr) {
		return new xmg(xmg_from_string(ninputs, expr));
	}

	xmg* mig_int_decompose(unsigned ninputs, unsigned func) {
		return new xmg(mig_decompose(ninputs, func));
	}

	float compute_reward(const majesty::xmg& mig_orig, const majesty::xmg& mig_final) {
		float nonodes = static_cast<float>(mig_orig.nnodes());
		float nfnodes = static_cast<float>(mig_final.nnodes());
		return nfnodes / nonodes;
	}

	inline bool maj3_applies(const node& n) {
		return (n.in1 == n.in2) || (n.in1 == n.in3) || (n.in2 == n.in3);
	}

	inline pair<nodeid, bool> maj3_prop(const node& n) {
		if (n.in1 == n.in2) {
			if (is_c1(n) == is_c2(n)) {
				return make_pair(n.in1, is_c1(n));
			} else {
				return make_pair(n.in3, is_c3(n));
			}
		} else if (n.in1 == n.in3) {
			if (is_c1(n) == is_c3(n)) {
				return make_pair(n.in1, is_c1(n));
			} else {
				return make_pair(n.in2, is_c2(n));
			}
		} else if (n.in2 == n.in3) {
			if (is_c2(n) == is_c3(n)) {
				return make_pair(n.in2, is_c2(n));
			} else {
				return make_pair(n.in1, is_c1(n));
			}
		} else {
			return make_pair(n.ecrep, false);
		}
	}

	xmg* apply_maj3(const xmg& mig, nodeid id) {
		auto res = new xmg();


		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		unordered_map<nodeid, pair<nodeid, bool>> nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				res->create_input(is_c);
				nodemap[i] = make_pair(i, false);
			} else if (i == id) {
				nodemap[i] = maj3_prop(node);
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}

		const auto& outputs = mig.outputs();
		const auto& outcompl = mig.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res->create_output(np.first, np.second != outcompl[i], "out" + i);
		}

		return res;
	}

	xmg* apply_inv_prop(const xmg& mig, nodeid id) {
		auto res = new xmg();

		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		unordered_map<nodeid, pair<nodeid, bool>> nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (i == id) {
				if (is_pi(node)) {
					auto is_c = (is_pi_c(node) != true);
					nodemap[i] = make_pair(res->create_input(is_c), true);
				} else if (i == id) {
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					auto invnode = res->create(
						in1.first, (in1.second != is_c1(node)) != true,
						in2.first, (in2.second != is_c2(node)) != true,
						in3.first, (in3.second != is_c3(node)) != true);
					invnode.second = true;
					nodemap[i] = invnode;
				}
			} else if (is_pi(node)) {
				nodemap[i] = make_pair(res->create_input(is_pi_c(node)), false);
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}

		const auto& outputs = mig.outputs();
		const auto& outcompl = mig.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res->create_output(np.first, np.second != outcompl[i], "out" + i);
		}

		return res;
	}

	

	// Checks if two nodes share the same input with the same polarity.
	// Used to implement the associativity axiom.
	inline bool share_input_polarity(const node& n1, const node& n2) {
		return ( 
			(n1.in1 == n2.in1 && is_c1(n1) == is_c1(n2)) ||
			(n1.in1 == n2.in2 && is_c1(n1) == is_c2(n2)) ||
			(n1.in1 == n2.in3 && is_c1(n1) == is_c3(n2)) ||
			(n1.in2 == n2.in1 && is_c2(n1) == is_c1(n2)) ||
			(n1.in2 == n2.in2 && is_c2(n1) == is_c2(n2)) ||
			(n1.in2 == n2.in3 && is_c2(n1) == is_c3(n2)) ||
			(n1.in3 == n2.in1 && is_c3(n1) == is_c1(n2)) ||
			(n1.in3 == n2.in2 && is_c3(n1) == is_c2(n2)) ||
			(n1.in3 == n2.in3 && is_c3(n1) == is_c3(n2)) 
			);
	}

	// Returns the id of the common child shared by two nodes. NOTE: fails if no such child exists!
	inline pair<nodeid,bool> shared_input_polarity(const node& n1, const node& n2) {
		if (n1.in1 == n2.in1 && is_c1(n1) == is_c1(n2)) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in1 == n2.in2 && is_c1(n1) == is_c2(n2)) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in1 == n2.in3 && is_c1(n1) == is_c3(n2)) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in2 == n2.in1 && is_c2(n1) == is_c1(n2)) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in2 == n2.in2 && is_c2(n1) == is_c2(n2)) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in2 == n2.in3 && is_c2(n1) == is_c3(n2)) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in3 == n2.in1 && is_c3(n1) == is_c1(n2)) {
			return make_pair(n1.in3, is_c3(n1));
		} else if (n1.in3 == n2.in2 && is_c3(n1) == is_c2(n2)) {
			return make_pair(n1.in3, is_c3(n1));
		} else if (n1.in3 == n2.in3 && is_c3(n1) == is_c3(n2)) {
			return make_pair(n1.in3, is_c3(n1));
		} else {
			cerr << "Error: shared input not found";
			exit(1);
		}
	}

	xmg* swap(const xmg& mig, nodeid id1, nodeid id2) {
		auto res = new xmg();
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		
		// First, check if the call is not ambiguous.
		const auto& gp = nodes[id1];
		if (is_pi(gp)) { // Grandparent obviously may not be a PI
			return NULL;
		}
		auto oldgpchildren = get_children(gp);
		auto filt_parents = filter_nodes(oldgpchildren, [&nodes, &gp, id2](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (!is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == id2 || parent.in2 == id2 || parent.in3 == id2) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return NULL;
		} else if (filt_parents.size() > 1) {
			// The call is ambiguous: there are multiple parents for which the axiom applies. This means that there
			// are multiple options to swap. We do not allow for ambiguity.
			return NULL;
		}
		// If we reach this point, there is one unambigious child of the grandparent to swap with:
		// the child that is not the shared child and not the parent.
		const auto parentnp = filt_parents[0];
		// If the parent is complemented, we need to apply inverter propagation first.
		if (parentnp.second) {
			return NULL;
		}
		const auto& parent = nodes[parentnp.first];
		auto common_childnp = shared_input_polarity(gp, parent);
		auto swap_childnp = drop_child(drop_child(oldgpchildren, parentnp), common_childnp)[0];
		// Thew new inner (parent) node should retain the same nodes except for the grandchild. We should also add the swap child to it.
		// Similarly, the new outer (grandparent) should retain the same nodes except for the swap node. We add to grandchild to it.
		auto oldgrandchildren = get_children(parent);
		auto filtgrandchildren = filter_nodes(drop_child(oldgrandchildren, common_childnp), [id2](pair<nodeid, bool> child) {
			if (child.first == id2) {
				return true;
			}
			return false;
		});
		if (filtgrandchildren.size() == 0) { 
			// After removing the common child, there is no grandchild z (id2) left. This means that the network is trying to swap
			// the common child itself which is not allowed.
			return NULL;
		} else if (filtgrandchildren.size() == 2) {
			// NOTE: This may be ambiguous: removing the common child from the parent leaves us with M(y, - , z) where
			// z is id2. Now, if y = z, we're not sure which node we're referring to if they have opposite polarities.
			auto filtgp1 = filtgrandchildren[0];
			auto filtgp2 = filtgrandchildren[1];
			if (filtgp1.second != filtgp2.second) { 
				// To avoid ambiguity we use the non-complemented child. The system can use inverter propagation to select the
				// other child instead.
				if (filtgp1.second) {
					filtgrandchildren = { filtgp2 };
				} else {
					filtgrandchildren = { filtgp1 };
				}
			}
		}
		auto grandchildnp = filtgrandchildren[0];
		auto ynodenp = drop_child(drop_child(oldgrandchildren, common_childnp), grandchildnp)[0];
		vector<pair<nodeid,bool>> newgrandchildren = { common_childnp, swap_childnp, ynodenp };

		// Count the fanout of the parent node. If it's > 1, we need to duplicate it.
		auto fanout = 0u;
		for (const auto& node : nodes) {
			if (is_pi(node)) {
				continue;
			}
			if (node.in1 == parentnp.first || node.in2 == parentnp.first || node.in3 == parentnp.first) {
				++fanout;
			}
		}
		assert(fanout >= 1);

		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == parentnp.first) {
				if (fanout > 1 || is_po(node)) { // Duplicate
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					nodemap[i] = res->create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node)
					);
				} 
			} else if (i == id1) {
				// Create the new parent node first.
				auto pin1 = newgrandchildren[0];
				auto newpin1 = nodemap[pin1.first];
				auto pin2 = newgrandchildren[1];
				auto newpin2 = nodemap[pin2.first];
				auto pin3 = newgrandchildren[2];
				auto newpin3 = nodemap[pin3.first];
				auto newparent = res->create(
					newpin1.first, newpin1.second != pin1.second,
					newpin2.first, newpin2.second != pin2.second,
					newpin3.first, newpin3.second != pin3.second
				);
				auto newgpin1 = nodemap[common_childnp.first];
				auto newgpin2 = nodemap[grandchildnp.first];
				nodemap[i] = res->create(
					newgpin1.first, newgpin1.second != common_childnp.second,
					newgpin2.first, newgpin2.second != grandchildnp.second,
					newparent.first, newparent.second != parentnp.second
				);
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}

		const auto& outputs = mig.outputs();
		const auto& outcompl = mig.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res->create_output(np.first, np.second != outcompl[i], "out" + i);
		}

		return res;
	}

	xmg* swap_ternary(const xmg& mig, nodeid gpid, nodeid pid, nodeid gcid) {
		auto res = new xmg();
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		
		// First, check if the call is not ambiguous.
		const auto& gp = nodes[gpid];
		if (is_pi(gp)) { // Grandparent obviously may not be a PI
			return NULL;
		}
		auto oldgpchildren = get_children(gp);
		auto filt_parents = filter_nodes(oldgpchildren, [&nodes, &gp, pid, gcid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == pid && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return NULL;
		} else if (filt_parents.size() > 1) {
			// The call is ambiguous: there are multiple parents for which the axiom applies. This means that there
			// are multiple options to swap. We do not allow for ambiguity.
			return NULL;
		}

		// If we reach this point, there is one unambigious child of the grandparent to swap with:
		// the child that is not the shared child and not the parent.
		const auto parentnp = filt_parents[0];
		// If the parent is complemented, we need to apply inverter propagation first.
		if (parentnp.second) {
			return NULL;
		}
		const auto& parent = nodes[parentnp.first];
		auto common_childnp = shared_input_polarity(gp, parent);
		auto swap_childnp = drop_child(drop_child(oldgpchildren, parentnp), common_childnp)[0];
		// Thew new inner (parent) node should retain the same nodes except for the grandchild. We should also add the swap child to it.
		// Similarly, the new outer (grandparent) should retain the same nodes except for the swap node. We add to grandchild to it.
		auto oldgrandchildren = get_children(parent);
		auto filtgrandchildren = filter_nodes(drop_child(oldgrandchildren, common_childnp), [gcid](pair<nodeid, bool> child) {
			if (child.first == gcid) {
				return true;
			}
			return false;
		});
		if (filtgrandchildren.size() == 0) { 
			// After removing the common child, there is no grandchild z (id2) left. Either the common child was the grandchild
			// to be swapped, or the given grandchild was never actually a grandchild. 
			return NULL;
		} else if (filtgrandchildren.size() == 2) {
			// NOTE: This may be ambiguous: removing the common child from the parent leaves us with M(y, - , z) where
			// z is id2. Now, if y = z, we're not sure which node we're referring to if they have opposite polarities.
			auto filtgp1 = filtgrandchildren[0];
			auto filtgp2 = filtgrandchildren[1];
			if (filtgp1.second != filtgp2.second) {
				// To avoid ambiguity we use the non-complemented child. The system can use inverter propagation to select the
				// other child instead.
				if (filtgp1.second) {
					filtgrandchildren = { filtgp2 };
				} else {
					filtgrandchildren = { filtgp1 };
				}
			}
		}
		auto grandchildnp = filtgrandchildren[0];
		auto ynodenp = drop_child(drop_child(oldgrandchildren, common_childnp), grandchildnp)[0];
		vector<pair<nodeid,bool>> newgrandchildren = { common_childnp, swap_childnp, ynodenp };

		// Count the fanout of the parent node. If it's > 1, we need to duplicate it.
		auto fanout = 0u;
		for (const auto& node : nodes) {
			if (is_pi(node)) {
				continue;
			}
			if (node.in1 == parentnp.first || node.in2 == parentnp.first || node.in3 == parentnp.first) {
				++fanout;
			}
		}
		assert(fanout >= 1);

		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == parentnp.first) {
				if (fanout > 1 || is_po(node)) { // Duplicate
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					nodemap[i] = res->create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node)
					);
				} 
			} else if (i == gpid) {
				// Create the new parent node first.
				auto pin1 = newgrandchildren[0];
				auto newpin1 = nodemap[pin1.first];
				auto pin2 = newgrandchildren[1];
				auto newpin2 = nodemap[pin2.first];
				auto pin3 = newgrandchildren[2];
				auto newpin3 = nodemap[pin3.first];
				auto newparent = res->create(
					newpin1.first, newpin1.second != pin1.second,
					newpin2.first, newpin2.second != pin2.second,
					newpin3.first, newpin3.second != pin3.second
				);
				auto newgpin1 = nodemap[common_childnp.first];
				auto newgpin2 = nodemap[grandchildnp.first];
				nodemap[i] = res->create(
					newgpin1.first, newgpin1.second != common_childnp.second,
					newgpin2.first, newgpin2.second != grandchildnp.second,
					newparent.first, newparent.second != parentnp.second
				);
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}

		const auto& outputs = mig.outputs();
		const auto& outcompl = mig.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res->create_output(np.first, np.second != outcompl[i], "out" + i);
		}

		return res;
	}

	xmg* apply_move(const xmg& mig, move& move) {
		switch (move.type) {
		case MAJ3_PROP:
			if (maj3_applies(mig.nodes()[move.nodeid1])) {
				return apply_maj3(mig, move.nodeid1);
			} else {
				return NULL;
			}
		case INVERTER_PROP:
			return apply_inv_prop(mig, move.nodeid1);
			break;
		case SWAP:
			return swap(mig, move.nodeid1, move.nodeid2);
			break;
		case SWAP_TERNARY:
			return swap_ternary(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case MAJ3_XXY:
		case MAJ3_XYY:
		default:
			return NULL;
			break;
		}
	}

	bool swap_applies(const vector<node>& nodes, nodeid gpid, nodeid z) {
		const auto& gp = nodes[gpid];
		if (is_pi(gp)) { // Grandparent obviously may not be a PI
			return false;
		}
		auto gpchildren = get_children(gp);
		auto filt_parents = filter_nodes(gpchildren, [&nodes, &gp, z](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (!is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == z || parent.in2 == z || parent.in3 == z) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return false;
		} else if (filt_parents.size() > 1) {
			// The call is ambiguous: there are multiple parents for which the axiom applies. This means that there
			// are multiple options to swap. We do not allow for ambiguity.
			return false;
		}
		const auto parentnp = filt_parents[0];
		// If the parent is complemented, we need to apply inverter propagation first.
		if (parentnp.second) {
			return false;
		}

		const auto& parent = nodes[parentnp.first];
		auto common_childnp = shared_input_polarity(gp, parent);
		auto swap_childnp = drop_child(drop_child(gpchildren, parentnp), common_childnp)[0];
		// Thew new inner (parent) node should retain the same nodes except for the grandchild. We should also add the swap child to it.
		// Similarly, the new outer (grandparent) should retain the same nodes except for the swap node. We add to grandchild to it.
		auto oldgrandchildren = get_children(parent);
		auto filtgrandchildren = filter_nodes(drop_child(oldgrandchildren, common_childnp), [z](pair<nodeid, bool> child) {
			if (child.first == z) {
				return true;
			}
			return false;
		});
		if (filtgrandchildren.size() == 0) { 
			// After removing the common child, there is no grandchild z (id2) left. This means that the network is trying to swap
			// the common child itself which is not allowed.
			return false;
		} else if (filtgrandchildren.size() > 1) {
			// NOTE: This may be ambiguous: removing the common child from the parent leaves us with M(y, - , z) where
			// z is id2. Now, if y = z, we're not sure which node we're referring to if they have opposite polarities.
			// For now we just return the non-complemented one.
			auto filtgp1 = filtgrandchildren[0];
			auto filtgp2 = filtgrandchildren[1];
			if (filtgp1.second != filtgp2.second) {
				return true;
			}
		}

		return true;
	}



	bool swap_ternary_applies(const vector<node>& nodes, nodeid gpid, nodeid pid, nodeid gcid) {
		const auto& gp = nodes[gpid];
		if (is_pi(gp)) { // Grandparent obviously may not be a PI
			return false;
		}
		auto gpchildren = get_children(gp);
		auto filt_parents = filter_nodes(gpchildren, [&nodes, &gp, pid, gcid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == pid && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return false;
		} else if (filt_parents.size() > 1) {
			// The call is ambiguous: there are multiple parents for which the axiom applies. This means that there
			// are multiple options to swap. We do not allow for ambiguity.
			return false;
		}
		const auto parentnp = filt_parents[0];
		// If the parent is complemented, we need to apply inverter propagation first.
		if (parentnp.second) {
			return false;
		}

		const auto& parent = nodes[parentnp.first];
		auto common_childnp = shared_input_polarity(gp, parent);
		auto swap_childnp = drop_child(drop_child(gpchildren, parentnp), common_childnp)[0];
		// Thew new inner (parent) node should retain the same nodes except for the grandchild. We should also add the swap child to it.
		// Similarly, the new outer (grandparent) should retain the same nodes except for the swap node. We add to grandchild to it.
		auto oldgrandchildren = get_children(parent);
		auto filtgrandchildren = filter_nodes(drop_child(oldgrandchildren, common_childnp), [gcid](pair<nodeid, bool> child) {
			if (child.first == gcid) {
				return true;
			}
			return false;
		});
		if (filtgrandchildren.size() == 0) { 
			// After removing the common child, there is no grandchild z (id2) left. Either the common child was the grandchild
			// to be swapped, or the given grandchild was never actually a grandchild.
			return false;
		} else if (filtgrandchildren.size() > 1) {
			// NOTE: This may be ambiguous: removing the common child from the parent leaves us with M(y, - , z) where
			// z is id2. Now, if y = z, we're not sure which node we're referring to if they have opposite polarities.
			// For now we just return the non-complemented one.
			auto filtgp1 = filtgrandchildren[0];
			auto filtgp2 = filtgrandchildren[1];
			if (filtgp1.second != filtgp2.second) {
				return true;
			}
		}

		return true;
	}

	bool dist_left_right_applies(const vector<node>& nodes, nodeid outnodeid, nodeid innodeid, nodeid distnodeid) {
		const auto& outnode = nodes[outnodeid];
		if (is_pi(outnode)) { // Grandparent obviously may not be a PI
			return false;
		}
		auto outnodechildren = get_children(outnode);
		auto filtered_innernodes = filter_nodes(outnodechildren, [&nodes, &outnode, innodeid, distnodeid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == innodeid && !is_pi(parent)) {
				if (parent.in1 == distnodeid || parent.in2 == distnodeid || parent.in3 == distnodeid) {
					return true;
				}
			}
			return false;
		});
		if (filtered_innernodes.size() == 0) {
			// The specified inner node is not both a child of the outer node and a parent of the grandchild.
			return false;
		} 
		// The call may be ambiguous: there may be multiple inner nodes for which the axiom applies. We do not 
		// 	allow for ambiguity. The axiom only applies if there is  a non-complemented inner nodes.
		bool have_non_complemented = false;
		for (const auto& np : filtered_innernodes) {
			if (!np.second) {
				have_non_complemented = true;
			}
		}
		if (!have_non_complemented) {
			return false;
		}
		const auto innernodep = make_pair(innodeid, false);

		const auto& innernode = nodes[innernodep.first];
		auto oldinnerchildren = get_children(innernode);
		auto filtinnerchildren = filter_nodes(oldinnerchildren, [distnodeid](pair<nodeid, bool> child) {
			if (child.first == distnodeid) {
				return true;
			}
			return false;
		});
		if (filtinnerchildren.size() > 1) {
			// NOTE: This may be ambiguous: the distchild occurs multiple times in the inner node. In this case
			// there is ambiguity if the child occurs in different polarities. We select the non-complemented version
			// of the child. It should be possible to select the other one with inverter propgation.
			return true;
		}
		return true;
	}

	xmg* apply_dist_left_right(const xmg& mig, nodeid outnodeid, nodeid innodeid, nodeid distnodeid) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		const auto& outnode = nodes[outnodeid];
		if (is_pi(outnode)) { // Grandparent obviously may not be a PI
			return NULL;
		}
		auto outnodechildren = get_children(outnode);
		auto filtered_innernodes = filter_nodes(outnodechildren, [&nodes, &outnode, innodeid, distnodeid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == innodeid && !is_pi(parent)) {
				if (parent.in1 == distnodeid || parent.in2 == distnodeid || parent.in3 == distnodeid) {
					return true;
				}
			}
			return false;
		});
		if (filtered_innernodes.size() == 0) {
			// The specified inner node is not both a child of the outer node and a parent of the grandchild.
			return NULL;
		}
		// The call may be ambiguous: there may be multiple inner nodes for which the axiom applies. We do not 
		// 	allow for ambiguity. The axiom only applies if there is  a non-complemented inner nodes.
		bool have_non_complemented = false;
		for (const auto& np : filtered_innernodes) {
			if (!np.second) {
				have_non_complemented = true;
			}
		}
		if (!have_non_complemented) {
			return NULL;
		}
		const auto innernodep = make_pair(innodeid, false);
		const auto& innode = nodes[innodeid];
		const auto outer_remainder_nodes = drop_child(outnodechildren, innernodep);
		assert(outer_remainder_nodes.size() == 2);

		auto oldinnerchildren = get_children(innode);
		auto filtinnerchildren = filter_nodes(oldinnerchildren, [distnodeid](pair<nodeid, bool> child) {
			if (child.first == distnodeid) {
				return true;
			}
			return false;
		});
		pair<nodeid, bool> distnodep;
		if (filtinnerchildren.size() > 1) {
			// NOTE: This may be ambiguous: the distchild occurs multiple times in the inner node. In this case
			// we need to decide which one to select. For now we just pick the first one.
			distnodep = filtinnerchildren[0];
		} else {
			distnodep = filtinnerchildren[0];
		}
		const auto inner_remainder_nodes = drop_child(oldinnerchildren, distnodep);
		assert(inner_remainder_nodes.size() == 2);

		// Count the fanout of the parent node. If it's > 1, we need to duplicate it.
		auto fanout = 0u;
		for (const auto& node : nodes) {
			if (is_pi(node)) {
				continue;
			}
			if (node.in1 == innodeid|| node.in2 == innodeid || node.in3 == innodeid) {
				++fanout;
			}
		}
		assert(fanout >= 1);
		auto duplicate = (fanout > 1 || is_po(innode));

		auto res = new xmg();
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == innodeid) {
				if (duplicate) {
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					nodemap[i] = res->create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node)
					);
				} 
			} else if (i == outnodeid) {
				// Create the two new inner children first.
				auto new_outer_remainderp1 = nodemap[outer_remainder_nodes[0].first];
				auto new_outer_remainderp2 = nodemap[outer_remainder_nodes[1].first];

				auto new_inner_remainderp1 = nodemap[inner_remainder_nodes[0].first];
				auto new_inner_remainderp2 = nodemap[inner_remainder_nodes[1].first];

				auto new_inner1 = res->create(
					new_outer_remainderp1.first, outer_remainder_nodes[0].second != new_outer_remainderp1.second,
					new_outer_remainderp2.first, outer_remainder_nodes[1].second != new_outer_remainderp2.second,
					new_inner_remainderp1.first, new_inner_remainderp1.second != inner_remainder_nodes[0].second
				);

				auto new_inner2 = res->create(
					new_outer_remainderp1.first, outer_remainder_nodes[0].second != new_outer_remainderp1.second,
					new_outer_remainderp2.first, outer_remainder_nodes[1].second != new_outer_remainderp2.second,
					new_inner_remainderp2.first, new_inner_remainderp2.second != inner_remainder_nodes[1].second
				);

				auto new_distnodep = nodemap[distnodep.first];
				nodemap[i] = res->create( new_inner1, new_inner2, new_distnodep );
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}

		const auto& outputs = mig.outputs();
		const auto& outcompl = mig.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res->create_output(np.first, np.second != outcompl[i], "out" + i);
		}

		return res;
	}

	vector<move> compute_moves(const xmg& mig) {
		vector<move> moves;

		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			{
				// Note that inverter propagation always applies.
				move move = {};
				move.type = INVERTER_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}
			if (is_pi(node)) {
				continue;
			}
			if (maj3_applies(node)) {
				move move = {};
				move.type = MAJ3_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}
			for (auto j = 0u; j < nnodes; j++) {
				if (swap_applies(nodes, i, j)) {
					move move = {};
					move.type = SWAP;
					move.nodeid1 = i;
					move.nodeid2 = j;
					moves.push_back(move);
				}
				for (auto k = 0u; k < nnodes; k++) {
					if (swap_ternary_applies(nodes, i, j, k)) {
						move move = {};
						move.type = SWAP_TERNARY;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					if (dist_left_right_applies(nodes, i, j, k)) {
						move move = {};
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
				}
			}
		}

		return moves;
	}
	
	xmg* get_optimum_mig(const xmg& mig) {
		auto func = simulate_xmg(mig);
		return new xmg(exact_mig(func));
	}

	xmg* strash_xmg(const xmg& mig) {
		return new xmg(strash(mig));
	}
}
