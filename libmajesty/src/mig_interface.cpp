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
			mig->create_output(outidx, outcompl);
		}

		return mig;
	}

	xmg* mig_manager::create_random_graph(unsigned ninputs, unsigned nnodes) {
		auto res = random_mig(rng, ninputs, nnodes, ninputs);
        res->create_dummy_names();
        return res;
	}
	
	xmg* mig_manager::create_random_graph(unsigned ninputs, unsigned nnodes, unsigned noutputs) {
		auto res = random_mig(rng, ninputs, nnodes, noutputs);
        res->create_dummy_names();
        return res;
	}

	xmg* mig_manager::random_mig_decomposition(unsigned ninputs) {
		auto num_funcs = (1ull << (1ull << ninputs));
		static uniform_int_distribution<mt19937_64::result_type> rule_dist(0, num_funcs-1);
		auto rand_func = rule_dist(rng);
        auto res = new xmg(mig_decompose(ninputs, rand_func));
        res->create_dummy_names();
        return res;
	}

	xmg* mig_string_decompose(const string& func) {
		auto ninputs = boost::integer_log2(func.size());
		return new xmg(mig_decompose(ninputs, func));
	}
	
	xmg* mig_expression_decompose(unsigned ninputs, const string& expr) {
		return new xmg(xmg_from_string(ninputs, expr));
	}

	xmg* mig_int_decompose(unsigned ninputs, unsigned func) {
		auto res = new xmg(mig_decompose(ninputs, func));
        res->create_dummy_names();
        return res;
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
			throw runtime_error("Error: shared input not found");
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
			// Note that we can have multiple fanouts to the same node!
			if (node.in1 == parentnp.first) {
				++fanout;
			}
			if (node.in2 == parentnp.first) {
				++fanout;
			}
			if (node.in3 == parentnp.first) {
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

	bool constructive_maj_applies(const vector<node>& nodes, nodeid parentnodeid, nodeid childnodeid, nodeid cnodeid) {
		// We need to make sure of the following:
		// (1) the parent node is not in the transitive fanin of the new child node
		// (2) the parent node is not a PI
		// (3) the specified child is indeed a child of the parent
		// Note that for now we don't allow construction of nodes when cnodeid > parentnodeid because it's tricky
		// to mainain topological order.
		const auto& parentnode = nodes[parentnodeid];
		if (is_pi(parentnode)) {
			return false;
		}
		if (parentnode.in1 == childnodeid || parentnode.in2 == childnodeid || parentnode.in3 == childnodeid) {
			return cnodeid < parentnodeid;
		} else {
			return false;
		}
	}

	// NOTE: We assume that the move has already been validated and do not check it again.
	// (Specifically we do not check for cycles.
	xmg* maj_xxy(const xmg& mig, nodeid parentnodeid, nodeid childnodeid, nodeid cnodeid) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();

		// A small point of ambiguity: the parent node may have multiple incoming edges from the child node.
		// At this moment we simply select the first one although better strategies may exist.
		auto parent_node = nodes[parentnodeid];
		auto child_nodes = get_children(parent_node);
		auto filt_child_nodes = filter_nodes(child_nodes, [&child_nodes, childnodeid](pair<nodeid, bool> ch) {
			if (ch.first == childnodeid) {
				return true;
			} 
			return false;
		});
		auto child_nodep = child_nodes[0];
		auto remaining_child_nodes = drop_child(child_nodes, child_nodep);
		assert(remaining_child_nodes.size() == 2);
		auto remaininc_child1p = remaining_child_nodes[0];
		auto remaininc_child2p = remaining_child_nodes[1];

		auto res = new xmg();
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == parentnodeid) {
				// Note that if cnodeid > parentnodeid we need to be careful to maintain topological order.
				// Create the new inner node first.
				auto new_child_nodep = nodemap[child_nodep.first];
				auto new_cnodep = nodemap[cnodeid];
				auto newinner = res->create(
					new_child_nodep.first, new_child_nodep.second != child_nodep.second,
					new_child_nodep.first, new_child_nodep.second != child_nodep.second,
					new_cnodep.first, new_cnodep.second
				);
				auto new_remaining_child1p = nodemap[remaininc_child1p.first];
				auto new_remaining_child2p = nodemap[remaininc_child2p.first];
				nodemap[i] = res->create(
					newinner.first, newinner.second,
					new_remaining_child1p.first, new_remaining_child1p.second != remaininc_child1p.second,
					new_remaining_child2p.first, new_remaining_child2p.second != remaininc_child2p.second
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

	// NOTE: We assume that the move has already been validated and do not check it again.
	// (Specifically we do not check for cycles.
	xmg* maj_xyy(const xmg& mig, nodeid parentnodeid, nodeid childnodeid, nodeid cnodeid) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();

		// A small point of ambiguity: the parent node may have multiple incoming edges from the child node.
		// At this moment we simply select the first one although better strategies may exist.
		auto parent_node = nodes[parentnodeid];
		auto child_nodes = get_children(parent_node);
		auto filt_child_nodes = filter_nodes(child_nodes, [&child_nodes, childnodeid](pair<nodeid, bool> ch) {
			if (ch.first == childnodeid) {
				return true;
			} 
			return false;
		});
		auto child_nodep = child_nodes[0];
		auto remaining_child_nodes = drop_child(child_nodes, child_nodep);
		assert(remaining_child_nodes.size() == 2);
		auto remaininc_child1p = remaining_child_nodes[0];
		auto remaininc_child2p = remaining_child_nodes[1];

		auto res = new xmg();
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == parentnodeid) {
				// Note that if cnodeid > parentnodeid we need to be careful to maintain topological order.
				// Create the new inner node first.
				auto new_child_nodep = nodemap[child_nodep.first];
				auto new_cnodep = nodemap[cnodeid];
				auto newinner = res->create(
					new_child_nodep.first, new_child_nodep.second != child_nodep.second,
					new_cnodep.first, new_cnodep.second,
					new_cnodep.first, new_cnodep.second != true
				);
				auto new_remaining_child1p = nodemap[remaininc_child1p.first];
				auto new_remaining_child2p = nodemap[remaininc_child2p.first];
				nodemap[i] = res->create(
					newinner.first, newinner.second,
					new_remaining_child1p.first, new_remaining_child1p.second != remaininc_child1p.second,
					new_remaining_child2p.first, new_remaining_child2p.second != remaininc_child2p.second
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
			if (np.first == pid && np.second == false && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return NULL;
		}

		// If we reach this point, there is one unambigious child of the grandparent to swap with:
		// the child that is not the shared child and not the parent.
		const auto parentnp = filt_parents[0];
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
			// Note that we make that the parent is not complemented. If it is, associativity doesn't apply and
			// we should try applying inverter propagation first.
			if (np.first == pid && !np.second && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() == 0) {
			// There is no child of the grandparent that has both a child in common and is a parent of the grandchild.
			return false;
		}
		const auto parentnp = filt_parents[0];

		const auto& parent = nodes[parentnp.first];
		auto common_childnp = shared_input_polarity(gp, parent);
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
			// For now we just return the first one.
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
			if (np.first == innodeid && !np.second && !is_pi(parent)) {
				if (parent.in1 == distnodeid || parent.in2 == distnodeid || parent.in3 == distnodeid) {
					return true;
				}
			}
			return false;
		});
		if (filtered_innernodes.size() == 0) {
			// The specified inner node is not a non-complemented child of the outer node and a parent of the grandchild.
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
			// there is ambiguity if the child occurs in different polarities. We select the first one
			// appear in the list of child nodes. 
			return true;
		}
		return true;
	}

	bool dist_right_left_applies(const vector<node>& nodes, nodeid outnodeid, nodeid distnodeid) {
		const auto& outnode = nodes[outnodeid];
		if (is_pi(outnode)) { // Grandparent obviously may not be a PI
			return false;
		}
		auto outnodechildren = get_children(outnode);

		// Note that this call is potentially ambiguous: distR->L applies iff the two inner children of the outer node
		// share two nodes. However there are (3 choose 2) = 3 ways for the inner nodes to share 2 children. We always
		// pick the first pair in this case.
		// Futhermore: distnodeid does not unambiguously identify the node over which to distribute: the outer node
		// may have both a complemented and a non-complemented node with the same id. We again pick the first one.

		auto filtered_outernodes = filter_nodes(outnodechildren, [&nodes, &outnode, distnodeid](pair<nodeid,bool> np) {
			if (np.first == distnodeid) {
				return true;
			}
			return false;
		});
		if (filtered_outernodes.size() == 0) {
			// The specified nodeid is not a child of the outer node.
			return false;
		} 
		auto distnodep = filtered_outernodes[0];
		auto remainingchildren = drop_child(outnodechildren, distnodep);
		assert(remainingchildren.size() == 2);
		for (const auto& ch : remainingchildren) {
			// We cannot allow the child to be complemented or to be a PI.
			auto chnode = nodes[ch.first];
			if (ch.second || is_pi(chnode)) {
				return false;
			}
		}
		vector<input> shared_children;
		
		auto remainingchild1p = remainingchildren[0];
		auto remainingchild1 = nodes[remainingchild1p.first];
		auto remainingchild2p = remainingchildren[1];
		auto remainingchild2 = nodes[remainingchild2p.first];
		auto remainingchild1children = get_children(remainingchild1);
		auto remainingchild2children = get_children(remainingchild2);
		for (const auto& ch1 : remainingchild1children) {
			for (const auto ch2 : remainingchild2children) {
				if (ch1.first == ch2.first && ch1.second == ch2.second) {
					shared_children.push_back(ch1);
					// Avoid duplicates.
					remainingchild2children = drop_child(remainingchild2children, ch1);
					break;
				}
			}
		}
		if (shared_children.size() < 2) {
			// The two inner nodes don't actually share 2 children to distribute.
			return false;
		}
		return true;
	}

	xmg* dist_right_left(const xmg& mig, nodeid outnodeid, nodeid distnodeid) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		const auto& outnode = nodes[outnodeid];
		if (is_pi(outnode)) { // Grandparent obviously may not be a PI
			return NULL;
		}
		auto outnodechildren = get_children(outnode);

		// Note that this call is potentially ambiguous: distR->L applies iff the two inner children of the outer node
		// share two nodes. However there are (3 choose 2) = 3 ways for the inner nodes to share 2 children. We always
		// pick the first pair in this case.
		// Futhermore: distnodeid does not unambiguously identify the node over which to distribute: the outer node
		// may have both a complemented and a non-complemented node with the same id. We again pick the first one.

		auto filtered_outernodes = filter_nodes(outnodechildren, [&nodes, &outnode, distnodeid](pair<nodeid,bool> np) {
			if (np.first == distnodeid) {
				return true;
			}
			return false;
		});
		if (filtered_outernodes.size() == 0) {
			// The specified nodeid is not a child of the outer node.
			return NULL;
		} 
		auto distnodep = filtered_outernodes[0];
		auto remainingchildren = drop_child(outnodechildren, distnodep);
		assert(remainingchildren.size() == 2);
		for (const auto& ch : remainingchildren) {
			// We cannot allow the child to be complemented or to be a PI.
			auto chnode = nodes[ch.first];
			if (ch.second || is_pi(chnode)) {
				return NULL;
			}
		}
		vector<input> shared_children;
		auto remainingchild1p = remainingchildren[0];
		auto remainingchild1 = nodes[remainingchild1p.first];
		auto remainingchild2p = remainingchildren[1];
		auto remainingchild2 = nodes[remainingchild2p.first];
		auto remainingchild1children = get_children(remainingchild1);
		auto remainingchild2children = get_children(remainingchild2);
		for (const auto& ch1 : remainingchild1children) {
			for (const auto ch2 : remainingchild2children) {
				if (ch1.first == ch2.first && ch1.second == ch2.second) {
					shared_children.push_back(ch1);
					// Avoid duplicates.
					remainingchild2children = drop_child(remainingchild2children, ch1);
					break;
				}
			}
		}
		// Reset after removing potential duplicates.
		remainingchild2children = get_children(remainingchild2);
		if (shared_children.size() < 2) {
			// The two inner nodes don't actually share 2 children to distribute.
			return NULL;
		}
		auto non_distchild1p = drop_child(drop_child(remainingchild1children, shared_children[0]), shared_children[1])[0];
		auto non_distchild2p = drop_child(drop_child(remainingchild2children, shared_children[0]), shared_children[1])[0];
		
		// Count the fanouts of the two remaining children. If they're > 1 or PO, we need to duplicate them.
		auto fanout1 = 0u, fanout2 = 0u;
		for (const auto& node : nodes) {
			if (is_pi(node)) {
				continue;
			}
			// Note that we can have multiple fanouts to the same node!
			if (node.in1 == remainingchild1p.first) {
				++fanout1;
			}
			if (node.in2 == remainingchild1p.first) {
				++fanout1;
			}
			if (node.in3 == remainingchild1p.first) {
				++fanout1;
			}
			if (node.in1 == remainingchild2p.first) {
				++fanout2;
			}
			if (node.in2 == remainingchild2p.first) {
				++fanout2;
			}
			if (node.in3 == remainingchild2p.first) {
				++fanout2;
			}
		}
		assert(fanout1 >= 1 && fanout2 >= 1);
		auto duplicate1 = (fanout1 > 1 || is_po(remainingchild1));
		auto duplicate2 = (fanout2 > 1 || is_po(remainingchild2));

		auto res = new xmg();
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == remainingchild1p.first) {
				if (duplicate1) {
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					nodemap[i] = res->create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node)
					);
				}
			} else if (i == remainingchild2p.first) {
				if (duplicate2) {
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
				auto newnondistchild1p = nodemap[non_distchild1p.first];
				auto newnondistchild2p = nodemap[non_distchild2p.first];
				auto newdistnodep = nodemap[distnodep.first];
				auto newdistchildp = res->create(
					newnondistchild1p.first, newnondistchild1p.second != non_distchild1p.second,
					newnondistchild2p.first, newnondistchild2p.second != non_distchild2p.second,
					newdistnodep.first, newdistnodep.second != distnodep.second
				);

				auto newremainingchild1p = nodemap[shared_children[0].first];
				auto newremainingchild2p = nodemap[shared_children[1].first];
				nodemap[i] = res->create(
					newremainingchild1p.first, newremainingchild1p.second != shared_children[0].second, 
					newremainingchild2p.first, newremainingchild2p.second != shared_children[1].second, 
					newdistchildp.first, newdistchildp.second 
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

	xmg* dist_left_right(const xmg& mig, nodeid outnodeid, nodeid innodeid, nodeid distnodeid) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		const auto& outnode = nodes[outnodeid];
		if (is_pi(outnode)) { // Grandparent obviously may not be a PI
			return NULL;
		}
		auto outnodechildren = get_children(outnode);
		auto filtered_innernodes = filter_nodes(outnodechildren, [&nodes, &outnode, innodeid, distnodeid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == innodeid && !np.second && !is_pi(parent)) {
				if (parent.in1 == distnodeid || parent.in2 == distnodeid || parent.in3 == distnodeid) {
					return true;
				}
			}
			return false;
		});
		if (filtered_innernodes.size() == 0) {
			// The specified inner node is not a non-complemented child of the outer node and a parent of the grandchild.
			return NULL;
		}
		const auto innernodep = filtered_innernodes[0];
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
			// Note that we can have multiple fanouts to the same node!
			if (node.in1 == innodeid) {
				++fanout;
			}
			if (node.in2 == innodeid) {
				++fanout;
			}
			if (node.in3 == innodeid) {
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
				new_distnodep.second = (new_distnodep.second != distnodep.second);
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
			if (is_pi(node)) {
				continue;
			}
			{
				// Note that inverter propagation always applies.
				move move;
				move.type = INVERTER_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}
			if (maj3_applies(node)) {
				move move;
				move.type = MAJ3_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}
			for (auto j = 0u; j < nnodes; j++) {
				/** 
				 * NOTE: Disabling binary swap for now.
				if (swap_applies(nodes, i, j)) {
					move move;
					move.type = SWAP;
					move.nodeid1 = i;
					move.nodeid2 = j;
					moves.push_back(move);
				}
				 */
				if (dist_right_left_applies(nodes, i, j)) {
					move move;
					move.type = DIST_RIGHT_LEFT;
					move.nodeid1 = i;
					move.nodeid2 = j;
					moves.push_back(move);
				}
				for (auto k = 0u; k < nnodes; k++) {
					if (swap_ternary_applies(nodes, i, j, k)) {
						move move;
						move.type = SWAP_TERNARY;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					if (dist_left_right_applies(nodes, i, j, k)) {
						move move;
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					if (constructive_maj_applies(nodes, i, j, k)) {
						move move;
						move.type = MAJ3_XXY;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
						move.type = MAJ3_XYY;
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

	xmg* get_optimum_xmg(const xmg& mig) {
		auto func = simulate_xmg(mig);
		return new xmg(exact_xmg(func));
	}

	xmg* strash_xmg(const xmg& mig) {
		return new xmg(strash(mig));
	}


	xmg* remove_duplicates(const xmg& mig) {
		return new xmg(rdup(mig));
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
		case SWAP_TERNARY:
			return swap_ternary(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case DIST_LEFT_RIGHT:
			return dist_left_right(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case DIST_RIGHT_LEFT:
			return dist_right_left(mig, move.nodeid1, move.nodeid2);
			break;
		case MAJ3_XXY:
			return maj_xxy(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case MAJ3_XYY:
			return maj_xyy(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		default:
			return NULL;
			break;
		}
	}

    xmg* verilog_to_xmg_ptr(const string& verilog_string) {
        return new xmg(verilog_to_xmg(verilog_string));
    }
}
