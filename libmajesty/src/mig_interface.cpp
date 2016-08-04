#include <mig_interface.h>
#include <xmg.h>
#include "strashmap.h"
#include <algorithm>

using namespace std;

namespace majesty {

	xmg random_mig(mt19937_64& rng, unsigned ninputs, unsigned nnodes) {
		static uniform_int_distribution<mt19937_64::result_type> rule_dist(0, 9);
		static uniform_int_distribution<mt19937_64::result_type> rule_sel(0, 1);
		static vector<nodeid> n1(3), c1(3);

		xmg mig;
		xmg_stats stats{
				0u, // Nr. strash hits
				0u, // nr_potentials
				0u, // nr_matches
				0u, // nr_misses
				0u, // nr_undefined
		};
		strashmap shmap(nnodes / 2, stats);
		// Create the "one" input
		mig.create_input();
		for (auto i = 0u; i < ninputs; i++) {
			mig.create_input();
		}
		while (mig.nnodes() < nnodes) {
			uniform_int_distribution<mt19937_64::result_type> node_dist(0, mig.nnodes() - 1);
			if (rule_dist(rng) == 0) { // Ensure that one of the rules applies
				// For associativity we need a random non-PI node, so check that
				// such a node exists. If not, make random Maj3 node instead.
				if (rule_sel(rng) == 0 || mig.nin() + 1 == mig.nnodes()) { // Maj3
					n1[0] = node_dist(rng);
					c1[0] = rule_sel(rng) == 0;
					n1[1] = n1[0];
					c1[1] = rule_sel(rng) == 0;
					n1[2] = node_dist(rng);
					c1[2] = rule_sel(rng) == 0;
					random_shuffle(n1.begin(), n1.end());
					random_shuffle(c1.begin(), c1.end());
					mig.rfind_or_create(n1[0], c1[0], n1[1], c1[1], n1[2], c1[2],
						shmap);
				} else {
					uniform_int_distribution<mt19937_64::result_type>
						nonpi_dist(mig.nin() + 1, mig.nnodes() - 1);
					n1[0] = nonpi_dist(rng);
					assert(!is_pi(mig.nodes()[n1[0]]));
					n1[1] = mig.nodes()[n1[0]].in1;
					n1[2] = node_dist(rng);
					c1[0] = rule_sel(rng) == 0;
					c1[1] = rule_sel(rng) == 0;
					c1[2] = rule_sel(rng) == 0;
					random_shuffle(n1.begin(), n1.end());
					random_shuffle(c1.begin(), c1.end());
					mig.rfind_or_create(n1[0], c1[0], n1[1], c1[1], n1[2], c1[2],
						shmap);
				}
			} else { // Make a random node
				mig.rfind_or_create(
					node_dist(rng), rule_sel(rng) == 0,
					node_dist(rng), rule_sel(rng) == 0,
					node_dist(rng), rule_sel(rng) == 0, shmap);
			}
		}

		return mig;
	}

	xmg mig_manager::create_random_graph(unsigned ninputs, unsigned nnodes) {
		return random_mig(rng, ninputs, nnodes);
	}

	float compute_reward(const majesty::xmg& mig_orig, const majesty::xmg& mig_final) {
		float nonodes = mig_orig.nnodes();
		float nfnodes = mig_final.nnodes();
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

		xmg_stats stats{
				0u, // Nr. strash hits
				0u, // nr_potentials
				0u, // nr_matches
				0u, // nr_misses
				0u, // nr_undefined
		};
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		strashmap shmap(nnodes / 2, stats);
		unordered_map<nodeid, pair<nodeid, bool>> nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				res->create_input();
				nodemap[i] = make_pair(i, false);
			} else if (i == id) {
				nodemap[i] = maj3_prop(node);
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->rfind_or_create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node), shmap);
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

		xmg_stats stats{
				0u, // Nr. strash hits
				0u, // nr_potentials
				0u, // nr_matches
				0u, // nr_misses
				0u, // nr_undefined
		};
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		strashmap shmap(nnodes / 2, stats);
		unordered_map<nodeid, pair<nodeid, bool>> nodemap;
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				res->create_input();
				nodemap[i] = make_pair(i, false);
			} else if (i == id) {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				auto invnode = res->rfind_or_create(
					in1.first, in1.second != true,
					in2.first, in2.second != true,
					in3.first, in3.second != true, shmap);
				invnode.second = true;
				nodemap[i] = invnode;
			} else {
				const auto& in1 = nodemap[node.in1];
				const auto& in2 = nodemap[node.in2];
				const auto& in3 = nodemap[node.in3];
				nodemap[i] = res->rfind_or_create(
					in1.first, in1.second != is_c1(node),
					in2.first, in2.second != is_c2(node),
					in3.first, in3.second != is_c3(node), shmap);
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

	xmg* apply_unary_move(const xmg& mig, UNARY_MOVE move_type, nodeid id) {
		const auto& node = mig.nodes()[id];
		if (is_pi(node)) {
			return NULL;
		}
		switch (move_type) {
		case MAJ3_LEFT_RIGHT:
			if (maj3_applies(node)) {
				return apply_maj3(mig, id);
			} else {
				return NULL;
			}
			break;
		case INVERTER_PROP:
			return apply_inv_prop(mig, id);
			break;
		default:
			return NULL;
			break;
		}
	}

	void mig_to_array(const xmg& mig, nodeid* narray) {
		const auto& nodes = mig.nodes();
		const auto nnodes = mig.nnodes();
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				assert(node.in1 == 0 && node.in2 == 0 && node.in3 == 0);
			}
			narray[(i * 6) + 0] = node.in1;
			narray[(i * 6) + 1] = is_c1(node);
			narray[(i * 6) + 2] = node.in2;
			narray[(i * 6) + 3] = is_c2(node);
			narray[(i * 6) + 4] = node.in3;
			narray[(i * 6) + 5] = is_c3(node);
		}
	}
}
