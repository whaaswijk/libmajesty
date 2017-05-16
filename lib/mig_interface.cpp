#include <mig_interface.h>
#include <xmg.h>
#include "strashmap.h"
#include <algorithm>
#include <iostream>
#include <truth_table_utils.hpp>
#include <convert.h>
#include <boost/pending/integer_log2.hpp>
#include <boost/optional.hpp>
#include <npn_canonization.hpp>
#include <string>
#include "maj_io.h"

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
		auto res = new xmg(mig_decompose(ninputs, func));
		res->create_dummy_names();
		return res;
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
		}

		return res;
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

	inline boost::optional<pair<nodeid,bool>> shared_input_polarity(const node& n1, const node& n2, const nodeid gcid) {
		if (n1.in1 == n2.in1 && is_c1(n1) == is_c1(n2) && n1.in1 != gcid) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in1 == n2.in2 && is_c1(n1) == is_c2(n2) && n1.in1 != gcid) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in1 == n2.in3 && is_c1(n1) == is_c3(n2) && n1.in1 != gcid) {
			return make_pair(n1.in1, is_c1(n1));
		} else if (n1.in2 == n2.in1 && is_c2(n1) == is_c1(n2) && n1.in2 != gcid) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in2 == n2.in2 && is_c2(n1) == is_c2(n2) && n1.in2 != gcid) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in2 == n2.in3 && is_c2(n1) == is_c3(n2) && n1.in2 != gcid) {
			return make_pair(n1.in2, is_c2(n1));
		} else if (n1.in3 == n2.in1 && is_c3(n1) == is_c1(n2) && n1.in3 != gcid) {
			return make_pair(n1.in3, is_c3(n1));
		} else if (n1.in3 == n2.in2 && is_c3(n1) == is_c2(n2) && n1.in3 != gcid) {
			return make_pair(n1.in3, is_c3(n1));
		} else if (n1.in3 == n2.in3 && is_c3(n1) == is_c3(n2) && n1.in3 != gcid) {
			return make_pair(n1.in3, is_c3(n1));
		} else {
			return boost::none;
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
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
		// It's possible that the grandparent has the parent as input twice! E.g. when the majority
		// rule applies at the grandparent. In that case, we cannot swap, because swapping would try
		// to make the parent a child of itself...
		auto nparentfanin = 0u;
		for (auto& gpchild : oldgpchildren) {
			if (gpchild.first == pid) {
				++nparentfanin;
			}
		}
		if (nparentfanin > 1u) {
			return NULL;
		}
		auto filt_parents = filter_nodes(oldgpchildren, [&nodes, &gp, pid, gcid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			if (np.first == pid && np.second == false && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});
		if (filt_parents.size() != 1) {
			// There is not exactly one uncomplemented child of the grandparent 
			// that has both a child in common and is a parent of the grandchild.
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
		}

		return res;
	}

	pair<nodeid, bool> node_relevance(xmg& xmg, input z, input x, input y, bool yc) {	
		if (z.first == x.first) { // Relevance applies!
			return make_pair(y.first, x.second != y.second != z.second != yc);
		}
		const auto node = xmg.nodes()[z.first];
		if (is_pi(node)) {
			return z;
		} else {
			auto in1 = node_relevance(xmg, make_pair(node.in1, is_c1(node)), x, y, yc);
			auto in2 = node_relevance(xmg, make_pair(node.in2, is_c2(node)), x, y, yc);
			auto in3 = node_relevance(xmg, make_pair(node.in3, is_c3(node)), x, y, yc);
			auto res = xmg.create(in1.first, in1.second, in2.first, in2.second, in3.first, in3.second);
			res.second = (res.second != z.second);
			return res;
		}
	}

	void select_node(const xmg& xmg, nodeid nid, vector<unsigned>& nref) {
		nref[nid]++;
		const auto& node = xmg.nodes()[nid];
		if (!is_pi(node)) {
			select_node(xmg, node.in1, nref);
			select_node(xmg, node.in2, nref);
			select_node(xmg, node.in3, nref);
		}
	}

	xmg* apply_substitution(const xmg& mig, nodeid nid, nodeid u, nodeid v) {
		const auto& node = mig.nodes()[nid];
		if (is_pi(node)) {
			return NULL;
		}
		if ((u > nid) || (v > nid) || (u == v)) {
			return NULL;
		}

		xmg tmp_res;
		vector<unsigned> nref;
		{
			const auto& nodes = mig.nodes();
			const auto nnodes = mig.nnodes();
			nodemap nodemap;
			for (auto i = 0u; i < nnodes; i++) {
				const auto& node = nodes[i];
				if (is_pi(node)) {
					nodemap[i] = make_pair(tmp_res.create_input(), false);
				} else if (i == nid) {
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					auto orig_node = tmp_res.create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node));

					auto unode = nodemap[u];
					auto vnode = nodemap[v];
					auto vu_node = node_relevance(tmp_res, orig_node, vnode, unode, false);
					auto vup_node = node_relevance(tmp_res, orig_node, vnode, unode, true);

					auto first_inner = tmp_res.create(vnode.first, vnode.second != true,
						vu_node.first, vu_node.second, unode.first, unode.second);
					auto second_inner = tmp_res.create(vnode.first, vnode.second != true,
						vup_node.first, vup_node.second, unode.first, unode.second != true);

					nodemap[i] = tmp_res.create(vnode, first_inner, second_inner);
				} else {
					const auto& in1 = nodemap[node.in1];
					const auto& in2 = nodemap[node.in2];
					const auto& in3 = nodemap[node.in3];
					nodemap[i] = tmp_res.create(
						in1.first, in1.second != is_c1(node),
						in2.first, in2.second != is_c2(node),
						in3.first, in3.second != is_c3(node));
				}
			}
			const auto& outputs = mig.outputs();
			const auto& outcompl = mig.outcompl();
			nref.resize(tmp_res.nnodes());
			for (auto i = 0u; i < outputs.size(); i++) {
				const auto nodeid = outputs[i];
				const auto np = nodemap[nodeid];
				select_node(tmp_res, np.first, nref);
				tmp_res.create_output(np.first, np.second != outcompl[i]);
			}
		}

		xmg res;
		const auto& nodes = tmp_res.nodes();
		const auto nnodes = tmp_res.nnodes();
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				nodemap[i] = make_pair(res.create_input(), false);
			} else if (nref[i] > 0) {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res.create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node));
			}
		}
		const auto& outputs = tmp_res.outputs();
		const auto& outcompl = tmp_res.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto nodeid = outputs[i];
			const auto np = nodemap[nodeid];
			res.create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res.add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res.add_outname(outname);
		}
	
		return remove_duplicates(res);
	}

	

	bool relevance_applies(const vector<node>& nodes, nodeid nid, nodeid x, nodeid y) {
		// We make sure of the following:
		// (1) the node is not a PI
		// (2) the specified x and y is indeed a child of the parent
		// (3) x != y, otherwise the operation is not interesting
		// Note that for now we don't allow construction of nodes when u > nid because it's tricky
		// to mainain topological order.
		const auto& node = nodes[nid];
		if (is_pi(node)) {
			return false;
		}

		if (x > nid || y > nid || x == y) {
			return false;
		}

		auto children = get_children(node);
		auto x_child = false, y_child = false;
		for (auto child : children) {
			if (child.first == x) {
				x_child = true;
			} else if (child.first == y) {
				y_child = true;
			}
		}
		if (!(x_child && y_child)) {
			return false;
		}
		return true;
	}

	xmg* apply_relevance(const xmg& mig, nodeid nid, nodeid x, nodeid y) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();

		const auto& node = nodes[nid];
		if (is_pi(node)) {
			return NULL;
		}
		if ((x > nid) || (y > nid) || (x == y)) {
			return NULL;
		}
		const auto children = get_children(node);
		auto x_childid = -1, y_childid = -1;
		for (auto i = 0u; i < 3; i++) {
			if (children[i].first == x) {
				x_childid = i;
				break;
			}
		}
		for (auto i = 0u; i < 3; i++) {
			if (children[i].first == y) {
				y_childid = i;
				break;
			}
		}
		if ((x_childid == -1) || (y_childid == -1)) {
			return NULL;
		}

		auto x_child = children[x_childid];
		auto y_child = children[y_childid];
		auto z_child = drop_child(drop_child(children, x_child), y_child)[0];

		auto res = new xmg;
		vector<unsigned> nref;
		nodemap nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				auto is_c = is_pi_c(node);
				nodemap[i] = make_pair(res->create_input(is_c), false);
			} else if (i == nid) {
				
				const auto& x_in = nodemap[x_child.first];
				const auto& y_in = nodemap[y_child.first];
				const auto& z_in = nodemap[z_child.first];
				const auto z_rel_in = node_relevance(*res, z_in, 
					make_pair(x_in.first, x_in.second != x_child.second), 
					make_pair(y_in.first, y_in.second != y_child.second), true);
				nodemap[i] = res->create(
					x_in.first, x_in.second != x_child.second, 
					y_in.first, y_in.second != y_child.second,
					z_rel_in.first, z_rel_in.second != z_child.second
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
		}
	

		return res;
	}

	inline bool substitution_applies(const vector<node>& nodes, nodeid nid, nodeid u, nodeid v) {
		// We make sure of the following:
		// (1) the node is not a PI
		// (2) the specified u and v is indeed a child of the parent
		// (3) u != v, otherwise the operation is not interesting
		// Note that for now we don't allow construction of nodes when u > nid because it's tricky
		// to mainain topological order.
		return (u < nid) && (v < nid) && (u != v);
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

	/*
	bool swap_ternary_applies_fast(const vector<node>& nodes, nodeid gpid, nodeid pid, nodeid gcid) {
		const auto& gp = nodes[gpid];
		const auto& pnode = nodes[pid];
		if (is_pi(pnode)) {
			return false;
		}

		// It's possible that the grandparent has the parent as input twice! E.g. when the majority
		// rule applies at the grandparent. In that case, we cannot swap, because swapping would try
		// to make the parent a child of itself...
		auto nparentfanin = 0;

		auto parentcmpl = true;

		if (gp.in1 == pid) {
			if (pnode.in1 == gcid || pnode.in2 == gcid || pnode.in3 == gcid) {
				++nparentfanin;
				parentcmpl = is_c1(gp);
			}
		}
		if (gp.in2 == pid) {
			if (pnode.in1 == gcid || pnode.in2 == gcid || pnode.in3 == gcid) {
				++nparentfanin;
				parentcmpl = is_c2(gp);
			}
		}
		if (gp.in3 == pid) {
			if (pnode.in1 == gcid || pnode.in2 == gcid || pnode.in3 == gcid) {
				++nparentfanin;
				parentcmpl = is_c3(gp);
			}
		}

		if (nparentfanin != 1 || parentcmpl) {
			return false;
		}

		auto shared = shared_input_polarity(gp, pnode, gcid);
		if (shared) {
			return true;
		} else {
			return false;
		}
		}*/

	bool swap_ternary_applies(const vector<node>& nodes, nodeid gpid, nodeid pid, nodeid gcid) {
		const auto& gp = nodes[gpid];
		if (is_pi(gp)) { // Grandparent obviously may not be a PI
			return false;
		}

		// It's possible that the grandparent has the parent as input twice! E.g. when the majority
		// rule applies at the grandparent. In that case, we cannot swap, because swapping would try
		// to make the parent a child of itself...
		auto gpchildren = get_children(gp);
		auto nparentfanin = 0;
		for (auto& gpchild : gpchildren) {
			if (gpchild.first == pid) {
				++nparentfanin;
			}
		}
		if (nparentfanin > 1) {
			return false;
		}

		auto filt_parents = filter_nodes(gpchildren, [&nodes, &gp, pid, gcid](pair<nodeid,bool> np) {
			auto parent = nodes[np.first];
			// Note that we make sure that the parent is not complemented. If it is, associativity doesn't apply and
			// we should try applying inverter propagation first.
			if (np.first == pid && !np.second && !is_pi(parent) && share_input_polarity(gp, parent)) {
				if (parent.in1 == gcid || parent.in2 == gcid || parent.in3 == gcid) {
					return true;
				}
			}
			return false;
		});

		if (filt_parents.size() != 1) {
			// There is not exactly one uncomplemented child of the grandparent that has both a child in 
			// common and is a parent of the grandchild.
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

	bool dist_left_right_applies_fast(const vector<node>& nodes, nodeid outnodeid, nodeid innodeid, nodeid distnodeid) {
		const auto& outnode = nodes[outnodeid];
		const auto& innode = nodes[innodeid];
		if (is_pi(innode))
			return false;
		
		auto ninnodes = 0;
		if (outnode.in1 == innodeid && !is_c1(outnode))
			++ninnodes;
		if (outnode.in2 == innodeid && !is_c2(outnode))
			++ninnodes;
		if (outnode.in3 == innodeid && !is_c3(outnode))
			++ninnodes;
		if (ninnodes != 1)
			return false;
		
		return innode.in1 == distnodeid || innode.in2 == distnodeid || innode.in3 == distnodeid;
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
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
			res->create_output(np.first, np.second != outcompl[i]);
		}

		for (const auto& inname : mig.innames()) {
			res->add_inname(inname);
		}
		for (const auto& outname : mig.outnames()) {
			res->add_outname(outname);
		}
		
		return res;
	}

	vector<move> compute_moves(const xmg& mig) {
		vector<move> moves;

        {
		    // Always add the identity
		    move move;
            move.type = IDENTITY;
            move.nodeid1 = 0;
			move.nodeid2 = 0;
			move.nodeid3 = 0;
            moves.push_back(move);
        }
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			move move;
			move.nodeid1 = 0;
			move.nodeid2 = 0;
			move.nodeid3 = 0;
			if (is_pi(node)) {
				continue;
			}
			{
				// Note that inverter propagation always applies.
				move.type = INVERTER_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}
			if (maj3_applies(nodes, node)) {
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
					move.type = DIST_RIGHT_LEFT;
					move.nodeid1 = i;
					move.nodeid2 = j;
					moves.push_back(move);
				}
				for (auto k = 0u; k < nnodes; k++) {
					if (swap_ternary_applies(nodes, i, j, k)) {
						move.type = SWAP_TERNARY;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					if (dist_left_right_applies(nodes, i, j, k)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}

					if (relevance_applies(nodes, i, j, k)) {
						move.type = RELEVANCE;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}

					if (substitution_applies(nodes, i, j, k)) {
						move.type = SUBSTITUTION;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					/*
					if (constructive_maj_applies(nodes, i, j, k)) {
						move.type = MAJ3_XXY;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
						move.type = MAJ3_XYY;
						moves.push_back(move);
					}
					*/
				}
			}
		}

		return moves;
	}

	vector<move> compute_partial_moves_exhaustive(const xmg& mig) {
		vector<move> moves;

		{
			move m;
			m.type = IDENTITY;
			m.nodeid1 = 0;
			moves.push_back(m);
		}

		const auto& nodes = mig.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& n = nodes[i];
			if (maj3_applies(nodes, n)) {
				move m;
				m.type = MAJ3_PROP;
				m.nodeid1 = i;
				moves.push_back(m);
			}
			if (pm_start_inv_prop(nodes, n)) {
				move m;
				m.type = INVERTER_PROP;
				m.nodeid1 = i;
				moves.push_back(m);
			}
			for (auto j = 0u; j < nodes.size(); j++) {
				if (pm_start_dist_right_left(nodes, n)) {
					partial_move pm;
					pm.c_move.type = DIST_RIGHT_LEFT;
					pm.c_move.nodeid1 = i;
					pm.filled = 1;
					if (partial_move_applies(nodes, j, pm)) {
						move m;
						m.type = DIST_RIGHT_LEFT;
						m.nodeid1 = pm.c_move.nodeid1;
						m.nodeid2 = j;
						moves.push_back(m);
					}
				}
				if (pm_start_ternary_swap(nodes, n)) {
					partial_move pm;
					pm.c_move.type = SWAP_TERNARY;
					pm.c_move.nodeid1 = i;
					pm.filled = 1;
					if (partial_move_applies(nodes, j, pm)) {
						partial_move pmm;
						pmm.c_move = pm.c_move;
						pmm.c_move.nodeid2 = j;
						pmm.filled = 2;
						for (auto k = 0u; k < nodes.size(); k++) {
							if (partial_move_applies(nodes, k, pmm)) {
								move m = pmm.c_move;
								m.nodeid3 = k;
								moves.push_back(m);
							}
						}
					}
				}
				if (pm_start_dist_left_right(nodes, n)) {
					partial_move pm;
					pm.c_move.type = DIST_LEFT_RIGHT;
					pm.c_move.nodeid1 = i;
					pm.filled = 1;
					if (partial_move_applies(nodes, j, pm)) {
						partial_move pmm;
						pmm.c_move = pm.c_move;
						pmm.c_move.nodeid2 = j;
						pmm.filled = 2;
						for (auto k = 0u; k < nodes.size(); k++) {
							if (partial_move_applies(nodes, k, pmm)) {
								move m = pmm.c_move;
								m.nodeid3 = k;
								moves.push_back(m);
							}
						}
					}
				}
				if (pm_start_substitution(nodes, n)) {
					partial_move pm;
					pm.c_move.type = SUBSTITUTION;
					pm.c_move.nodeid1 = i;
					pm.filled = 1;
					if (partial_move_applies(nodes, j, pm)) {
						partial_move pmm;
						pmm.c_move = pm.c_move;
						pmm.c_move.nodeid2 = j;
						pmm.filled = 2;
						for (auto k = 0u; k < nodes.size(); k++) {
							if (partial_move_applies(nodes, k, pmm)) {
								move m = pmm.c_move;
								m.nodeid3 = k;
								moves.push_back(m);
							}
						}
					}
				}
				if (pm_start_relevance(nodes, n)) {
					partial_move pm;
					pm.c_move.type = RELEVANCE;
					pm.c_move.nodeid1 = i;
					pm.filled = 1;
					if (partial_move_applies(nodes, j, pm)) {
						partial_move pmm;
						pmm.c_move = pm.c_move;
						pmm.c_move.nodeid2 = j;
						pmm.filled = 2;
						for (auto k = 0u; k < nodes.size(); k++) {
							if (partial_move_applies(nodes, k, pmm)) {
								move m = pmm.c_move;
								m.nodeid3 = k;
								moves.push_back(m);
							}
						}
					}
				}
			}
		}


		return moves;
	}

	vector<move> compute_moves_fast(const xmg& mig, const unsigned max_nr_moves) {
		vector<move> moves;
		auto max_nr_moves_exceeded = false;
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();

		{
		    // Always add the identity
		    move move;
            move.type = IDENTITY;
            move.nodeid1 = 0;
			move.nodeid2 = 0;
			move.nodeid3 = 0;
            moves.push_back(move);
        }
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			}
			move move;
			move.nodeid1 = 0;
			move.nodeid2 = 0;
			move.nodeid3 = 0;
			{
				// Note that inverter propagation always applies.
				move.type = INVERTER_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}

			if (maj3_applies(nodes, node)) {
				move.type = MAJ3_PROP;
				move.nodeid1 = i;
				moves.push_back(move);
			}

			if (dist_right_left_applies(nodes, i, node.in1)) {
				move.type = DIST_RIGHT_LEFT;
				move.nodeid1 = i;
				move.nodeid2 = node.in1;
				moves.push_back(move);
			}
			if (dist_right_left_applies(nodes, i, node.in2)) {
				move.type = DIST_RIGHT_LEFT;
				move.nodeid1 = i;
				move.nodeid2 = node.in2;
				moves.push_back(move);
			}
			if (dist_right_left_applies(nodes, i, node.in3)) {
				move.type = DIST_RIGHT_LEFT;
				move.nodeid1 = i;
				move.nodeid2 = node.in3;
				moves.push_back(move);
			}
			{
				auto innode1 = nodes[node.in1];
				if (swap_ternary_applies(nodes, i, node.in1, innode1.in1)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in1;
					move.nodeid3 = innode1.in1;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in1, innode1.in2)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in1;
					move.nodeid3 = innode1.in2;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in1, innode1.in3)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in1;
					move.nodeid3 = innode1.in3;
					moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in1, innode1.in1)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in1;
						move.nodeid3 = innode1.in1;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in1, innode1.in2)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in1;
						move.nodeid3 = innode1.in2;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in1, innode1.in3)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in1;
						move.nodeid3 = innode1.in3;
						moves.push_back(move);
				}
			}
			{
				auto innode2 = nodes[node.in2];
				if (swap_ternary_applies(nodes, i, node.in2, innode2.in1)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in2;
					move.nodeid3 = innode2.in1;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in2, innode2.in2)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in2;
					move.nodeid3 = innode2.in2;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in2, innode2.in3)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in2;
					move.nodeid3 = innode2.in3;
					moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in2, innode2.in1)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in2;
						move.nodeid3 = innode2.in1;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in2, innode2.in2)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in2;
						move.nodeid3 = innode2.in2;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in2, innode2.in3)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in2;
						move.nodeid3 = innode2.in3;
						moves.push_back(move);
				}
			}
			{
				auto innode3 = nodes[node.in3];
				if (swap_ternary_applies(nodes, i, node.in3, innode3.in1)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in3;
					move.nodeid3 = innode3.in1;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in3, innode3.in2)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in3;
					move.nodeid3 = innode3.in2;
					moves.push_back(move);
				}
				if (swap_ternary_applies(nodes, i, node.in3, innode3.in3)) {
					move.type = SWAP_TERNARY;
					move.nodeid1 = i;
					move.nodeid2 = node.in3;
					move.nodeid3 = innode3.in3;
					moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in3, innode3.in1)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in3;
						move.nodeid3 = innode3.in1;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in3, innode3.in2)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in3;
						move.nodeid3 = innode3.in2;
						moves.push_back(move);
				}
				if (dist_left_right_applies_fast(nodes, i, node.in3, innode3.in3)) {
						move.type = DIST_LEFT_RIGHT;
						move.nodeid1 = i;
						move.nodeid2 = node.in3;
						move.nodeid3 = innode3.in3;
						moves.push_back(move);
				}
			}
		}

		for (auto i = 0u; i < nnodes && !max_nr_moves_exceeded; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			}
			move move;
			move.nodeid1 = 0;
			move.nodeid2 = 0;
			move.nodeid3 = 0;

			for (auto j = 0u; j < nnodes && !max_nr_moves_exceeded; j++) {
				for (auto k = 0u; k < nnodes && !max_nr_moves_exceeded; k++) {
					if (max_nr_moves > 0 && moves.size() >= max_nr_moves) {
						max_nr_moves_exceeded = true;
						break;
					}
					if (substitution_applies(nodes, i, j, k)) {
						move.type = SUBSTITUTION;
						move.nodeid1 = i;
						move.nodeid2 = j;
						move.nodeid3 = k;
						moves.push_back(move);
					}
					if (relevance_applies(nodes, i, j, k)) {
						move.type = RELEVANCE;
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

	bool partial_move_applies(const vector<node>& nodes, const nodeid nid, const partial_move& pm) {
		if (pm.c_move.type == DIST_RIGHT_LEFT) {
				// Only applies if the specified nid is a child of the node that has already been selected
				// AND the other two (uncomplemented) children share 2 inputs
				if (pm.filled != 1)
					return false;

				const auto& snode = nodes[pm.c_move.nodeid1];
				auto children = get_children(snode);
				// To allow for more cases, try to drop the complemented version of the child first.
				// If another, uncomplemented child remains, we might be able use it, but we can never
				// use complemented children.
				auto dchildren = drop_child(children, make_pair(nid, true));
				if (dchildren.size() == 3)
					dchildren = drop_child(children, nid);					
				if (dchildren.size() != 2) // Specified nid is not a child of selected node
					return false;
				
				const auto& ochildp1 = dchildren[0];
				const auto& ochildp2 = dchildren[1];
				// The children may not be complemented
				if (ochildp1.second || ochildp2.second)
					return false;
				const auto& ochild1 = nodes[ochildp1.first];
				const auto& ochild2 = nodes[ochildp2.first];
				return !is_pi(ochild1) && !is_pi(ochild2) && share_two_input_polarity(ochild1, ochild2);
		} else if (pm.c_move.type == SWAP_TERNARY) {
			if (pm.filled == 1) {
				// A node has been selected to swap on, we want to check if the
				// specified nid corresponds to an uncomplemented(!) child of that node with which
				// it has a grandchild in common
				const auto& snode = nodes[pm.c_move.nodeid1];
				const auto& pnode = nodes[nid];
				if (is_pi(pnode)) {
					return false;
				}
				auto children = get_children(snode);
				auto ochildren = drop_child(children, make_pair(nid, false));
				if (ochildren.size() != 2) // Not an uncomplemented child
					return false;
				return share_input_polarity(snode, pnode);
			} else if (pm.filled == 2) {
				return swap_ternary_applies(nodes, pm.c_move.nodeid1, pm.c_move.nodeid2, nid);
			} else {
				return false;
			}
		} else if (pm.c_move.type == DIST_LEFT_RIGHT) {
			if (pm.filled == 1) {
				// The specified nid should correspond to a non-complemented non-PI child
				// of the selected node.
				const auto& snode = nodes[pm.c_move.nodeid1];
				auto children = get_children(snode);
				auto ochildren = drop_child(children, make_pair(nid, false));
				if (ochildren.size() != 2) // Not an uncomplemented child
					return false;
				const auto& innode = nodes[nid];
				return !is_pi(innode);
			} else if (pm.filled == 2) {
				return dist_left_right_applies(nodes, pm.c_move.nodeid1,
											   pm.c_move.nodeid2, nid);
			} else {
				return false;
			}
		} else if (pm.c_move.type == SUBSTITUTION) {
			if (pm.filled == 1) {
				// > 0 because otherwise we cannot select a third argument
				return nid < pm.c_move.nodeid1;
			} else if (pm.filled == 2) {
				return nid < pm.c_move.nodeid1 && nid != pm.c_move.nodeid2;
			} else {
				return false;
			}
		} else if (pm.c_move.type == RELEVANCE) {
			if (pm.filled == 1) {
				// We need the specified nid to actually be a child of the selected node
				const auto& snode = nodes[pm.c_move.nodeid1];
				return snode.in1 == nid || snode.in2 == nid || snode.in3 == nid;
			} else if (pm.filled == 2) {
				const auto& snode = nodes[pm.c_move.nodeid1];
				return (snode.in1 == nid || snode.in2 == nid || snode.in3 == nid) && nid != pm.c_move.nodeid2;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	
	xmg*
	get_optimum_mig(const xmg& mig)
	{
		auto func = simulate_xmg(mig);
		return new xmg(exact_mig(func));
	}

	xmg*
	get_depth_optimum_mig(const xmg& mig)
	{
		auto func = simulate_xmg(mig);
		return new xmg(exact_depth_mig(func));
	}

	xmg* get_optimum_xmg(const xmg& mig) {
		auto func = simulate_xmg(mig);
		return new xmg(exact_xmg(func));
	}

	xmg* strash_xmg(const xmg& mig, bool no_compl) {
		if (no_compl) {
			return new xmg(strash_no_compl(mig));
		} else {
			return new xmg(strash(mig));
		}
	}

	xmg* remove_duplicates(const xmg& mig) {
		return new xmg(rdup(mig));
	}

	xmg* apply_move(const xmg& mig, const move& move) {
		const auto& nodes = mig.nodes();
		switch (move.type) {
		case MAJ3_PROP:
			if (maj3_applies(nodes, nodes[move.nodeid1])) {
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
		case SUBSTITUTION:
			return apply_substitution(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case RELEVANCE:
			return apply_relevance(mig, move.nodeid1, move.nodeid2, move.nodeid3);
			break;
		case IDENTITY:
			return new xmg(mig);
		default:
			return NULL;
			break;
		}
	}

    xmg* verilog_to_xmg_ptr(const string& verilog_string) {
        return new xmg(verilog_to_xmg(verilog_string));
    }

	xmg* get_npn_representative(const xmg& xmg) {
		const auto tt = simulate_xmg(xmg);
		vector<unsigned> perm; cirkit::tt phase;
		auto npn_tt = cirkit::exact_npn_canonization(tt, phase, perm);
		return mig_int_decompose(xmg.nin(), npn_tt.to_ulong());
	}

	unsigned long get_truth_table(const xmg& xmg) {
		return simulate_xmg(xmg).to_ulong();
	}

	xmg*
	resyn2(const xmg& m)
	{
		write_verilog(m, "tmp.v");
		auto cmdstr = "abc -c \"r tmp.v; st; resyn2; w tmp.v; quit \" > /dev/null";
		auto success = system(cmdstr);
		if (success != 0) {
			throw runtime_error("Synthesis through ABC failed");
		}
		auto tmp_xmg = read_verilog("tmp.v");
		return new xmg(tmp_xmg);
	}

	xmg*
	resyn2(const xmg& m, const string& tmpfilename)
	{
		write_verilog(m, tmpfilename);
		auto cmdstr = "abc -c \"r " + tmpfilename + "; st; resyn2; w " + tmpfilename + "; quit \" > /dev/null";
		auto success = system(cmdstr.c_str());
		if (success != 0) {
			throw runtime_error("Synthesis through ABC failed");
		}
		auto tmp_xmg = read_verilog(tmpfilename);
		return new xmg(tmp_xmg);
	}
	
	
}
