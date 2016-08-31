#include <iostream>
#include <stack>
#include "lut_optimize.h"
#include "strashmap.h"
#include "npn_canonization.hpp"
#include <convert.h>

extern "C" {
#include <bool/lucky/luckyInt.h>
#include <misc/util/utilTruth.h>
}

using namespace std;
using namespace cirkit;
namespace majesty {

	xmg* ptr_lut_area_strategy(const xmg& m, unsigned lut_size, unsigned nr_backtracks) {
		auto frparams = default_xmg_params();
		frparams->nr_backtracks = nr_backtracks;
		return new xmg(gen_lut_area_strategy(m, frparams.get(), lut_size, ALL));
	}

	xmg lut_area_strategy(const xmg& m, 
			const xmg_params* frparams, unsigned lut_size ) {
		return gen_lut_area_strategy(m, frparams, lut_size, ALL);
	}
	
	xmg mig_lut_area_strategy(const xmg& m, 
			const xmg_params* frparams, unsigned lut_size) {
		return gen_lut_area_strategy(m, frparams, lut_size, MAJ);
	}

	xmg gen_lut_area_strategy(const xmg& m, 
			const xmg_params* frparams, unsigned lut_size, NODE_TYPE type) {
		xmg cmig(m, frparams);
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;

		while (true) {
			auto oldsize = cmig.nnodes();
			const auto cut_map = 
				enumerate_cuts(cmig, cut_params.get());
			auto best_area = eval_matches_area(cmig, cut_map);
			auto area_cover = build_cover(cmig, best_area);
			it_exact_cover(cmig, area_cover, 
					cut_map, best_area);
			auto fm = compute_functions(cmig, 
					area_cover, best_area, cut_map);
			auto lutxmg = xmg_from_luts(cmig, area_cover, 
					best_area, fm, type);
			//auto frlutxmg = xmg(lutxmg, frparams);
			auto newsize = lutxmg.nnodes();
			if (newsize < oldsize) {
				cmig = std::move(lutxmg);
			} else {
				break;
			}
		}

		return cmig;
	}
	
	/*
	void mine_functions(const xmg& m, const xmg_params* p, 
			unsigned lut_size) {
		xmg xmg(m, p);
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;
		while (true) {
			auto oldsize = xmg.nnodes();
			const auto cut_map = 
				enumerate_cuts(xmg, cut_params.get());
		}
	}
	*/

	inline void inv(const char* perm, char* invperm) {
		for ( auto i = 0u; i < 16; ++i ) { invperm[(int)perm[i]] = i; }
	}

	pair<nodeid,bool> frmaj3_from_string(const string& expr, 
			unsigned offset, const bracket_map_t& majbrackets, 
			const bracket_map_t& xorbrackets,
			const input_map_t& imap, xmg& xmg, strashmap& shmap) {
		assert(expr[offset] != '!');
		if (expr[offset] == '<') {
			pair<nodeid,bool> children[3];

			auto child_pos = offset + 1u;
			auto to = 0u;
			auto inv = false;
			for (auto i = 0u; i < 3; i++) {
				child_pos += to;
				inv = expr[child_pos] == '!';
				if (inv)
					++child_pos;
				if (expr[child_pos] == '<') 
					to = majbrackets.at(child_pos) - child_pos + 1u;
				else if (expr[child_pos] == '[')
					to = xorbrackets.at(child_pos) - child_pos + 1u;
				else
					to = 1u;
				children[i] = frmaj3_from_string(expr, child_pos,
						majbrackets, xorbrackets, imap, xmg, shmap);
				children[i].second ^= inv;
			}
			return xmg.find_or_create(children[0].first, children[0].second,
					children[1].first, children[1].second, children[2].first,
					children[2].second, shmap);
		} else if (expr[offset] == '[') {
			pair<nodeid,bool> children[2];

			auto child_pos = offset + 1u;
			auto to = 0u;
			auto inv = false;
			for (auto i = 0u; i < 2; i++) {
				child_pos += to;
				inv = expr[child_pos] == '!';
				if (inv)
					++child_pos;
				if (expr[child_pos] == '<') 
					to = majbrackets.at(child_pos) - child_pos + 1u;
				else if (expr[child_pos] == '[') 
					to = xorbrackets.at(child_pos) - child_pos + 1u;
				else
					to = 1u;
				children[i] = frmaj3_from_string(expr, child_pos,
						majbrackets, xorbrackets, imap, xmg, shmap);
				children[i].second ^= inv;
			}
			return xmg.find_or_create(children[0].first, children[0].second,
					children[1].first, children[1].second, shmap);

		} else if (expr[offset] == '0') {
			return make_pair(0, true);
		} else if (expr[offset] == '1') {
			return make_pair(0, false);
		} else {
			return imap.at(expr[offset]);
		}
	}

	pair<nodeid,bool> decompose_cut(
			xmg& xmg, const cut* cut, tt& cutfunction, strashmap& shmap, 
			nodemap& nodemap, function_store& fstore, NODE_TYPE type) {
		const auto& cutnodes = cut->nodes();
		//const auto npn = exact_npn_canonization(cutfunction, phase, perm);
		//auto npn = fstore.npn_canon(cutfunction, phase, perm);
		auto num_vars = tt_num_vars(cutfunction);
		unsigned uCanonPhase; char pCanonPerm[16]; char invperm[16];
		auto npn = jake_canon(cutfunction, &uCanonPhase, pCanonPerm);
		auto num_npn_vars = tt_num_vars(npn);
		cout  << num_npn_vars << endl;
		const auto min_xmg = fstore.min_size_depth_xmg(npn, type);
		input_map_t imap;
		inv(pCanonPerm, invperm);
		for (auto i = 0u; i < cutnodes.size(); i++) {
			auto inode = nodemap[cutnodes[i]];
			imap['a' + invperm[i]] = make_pair(
					inode.first, (uCanonPhase & (1u << i)) ^ inode.second);
		}
		const auto majbrackets = find_bracket_pairs(min_xmg, '<', '>');
		const auto xorbrackets = find_bracket_pairs(min_xmg, '[', ']');
		auto inv = false;
		auto offset = 0u;
		if (min_xmg[0] == '!') {
			inv = true;
			offset = 1u;
		}
		auto res = frmaj3_from_string(min_xmg, offset, 
				majbrackets, xorbrackets, imap, xmg, shmap);
		res.second = (res.second != inv);
		res.second = (res.second != (uCanonPhase & (1 << num_vars)));
		return res;
	}

	xmg xmg_from_luts(const xmg& m, const cover& cover, 
			const bestmap& best, const funcmap& funcmap, NODE_TYPE type) {
		xmg n;
		function_store fstore;

		xmg_stats stats {
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};
		
		nodemap nodemap;
		strashmap shmap(m.nnodes()/2, stats);

		const auto& nodes = m.nodes();
		const auto& innames = m.innames();
		nodemap[0] = make_pair(n.create_input(), false);
		auto total_nodes = nodes.size();
		auto progress = 0u;
		for (auto i = 1u; i < total_nodes; i++) {
			if (cover.at(i) == 0) {
				// Node is not in cover, we don't consider it
				++progress;
				continue;
			}
			const auto& node = nodes[i];
			if (is_pi(node)) {
				nodemap[i] = make_pair(
						n.create_input(innames[i-1]), false);
				++progress;
				continue;
			}
			const auto& cut = best.at(i);
			auto& f = *funcmap.at(cut);
			nodemap[i] = decompose_cut(n, cut, f, shmap, nodemap, fstore, type);
			cout << "Progress: (" << ++progress << "/" << total_nodes << ")\r";
		}
		cout << endl;
		const auto& outputs = m.outputs();
		const auto& outcompls = m.outcompl();
		const auto& outnames = m.outnames();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto np = nodemap[outputs[i]];
			n.create_output(np.first, outcompls[i] != np.second, outnames[i]);
		}

		return n;
	}

	// NPN canonization functions from ABC
	tt jake_canon(const tt& ttf, unsigned* uCanonPhase, char* pCanonPerm) {
		auto num_vars = tt_num_vars(ttf);

		vector<word> pTruth( ttf.num_blocks() );
		boost::to_block_range( ttf, &pTruth[0] );

		resetPCanonPermArray( pCanonPerm, num_vars );
		*uCanonPhase = luckyCanonicizer_final_fast(&pTruth[0], num_vars, pCanonPerm);

		tt tt_npn( pTruth.begin(), pTruth.end() );

		for ( auto i = 0u; i < num_vars; ++i ) {
			pCanonPerm[i] -= 'a';
		}
		Abc_TtImplementNpnConfig(&pTruth[0], num_vars, pCanonPerm, *uCanonPhase);
		tt ttf_new( pTruth.begin(), pTruth.end() );
		assert( ttf == ttf_new );

		return tt_npn;
	}
}
