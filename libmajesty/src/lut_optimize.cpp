#include <iostream>
#include <stack>
#include "lut_optimize.h"
#include "strashmap.h"
#include "npn_canonization.hpp"
#include <convert.h>
#include <exact.h>

using namespace std;
using namespace cirkit;
using boost::optional;

namespace majesty {

	xmg* ptr_lut_area_strategy(const xmg& m, unsigned lut_size, unsigned nr_backtracks) {
		auto frparams = default_xmg_params();
		frparams->nr_backtracks = nr_backtracks;
		return new xmg(lut_area_strategy(m, frparams.get(), lut_size));
	}

	xmg lut_area_strategy(const xmg& m, const xmg_params* frparams, unsigned lut_size) {
        return lut_area_timeout_strategy(m, frparams, lut_size, 0, rebuild_cover).value();
	}
   	
	logic_ntk lut_area_strategy(const logic_ntk& ntk, unsigned lut_size) {
		return lut_area_timeout_strategy(ntk, lut_size, 0).value();
	}

	xmg* ptr_lut_area_timeout_strategy(const xmg& m, unsigned lut_size, unsigned timeout, unsigned nr_backtracks) {
		auto frparams = default_xmg_params();
		frparams->nr_backtracks = nr_backtracks;
		auto oxmg = lut_area_timeout_strategy(m, frparams.get(), lut_size, timeout, rebuild_cover);
		if (oxmg) {
			return new xmg(std::move(oxmg.get()));
		} else {
			return NULL;
		}
	}
	
	optional<xmg> lut_area_timeout_strategy(const xmg& m, unsigned lut_size, unsigned timeout) {
		auto frparams = default_xmg_params();
		return lut_area_timeout_strategy(m, frparams.get(), lut_size, timeout, rebuild_cover);
	}
	
	optional<xmg> lut_area_timeout_strategy(const xmg& m, unsigned lut_size, unsigned timeout, timeout_behavior behavior) {
		auto frparams = default_xmg_params();
		return lut_area_timeout_strategy(m, frparams.get(), lut_size, timeout, behavior);
	}
	
	boost::optional<logic_ntk> 
		lut_area_timeout_strategy(const logic_ntk& ntk, unsigned lut_size, unsigned backtrack_limit) {
		return lut_area_timeout_strategy(ntk, lut_size, backtrack_limit, rebuild_cover);
	}

	boost::optional<logic_ntk> 
		lut_area_timeout_strategy(const logic_ntk& ntk, unsigned lut_size, unsigned backtrack_limit, timeout_behavior tb) {
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;
		vector<tt> timeoutfuncs;
	
		bool timeout_occurred = false;
		auto oldsize = ntk.nnodes();
		funcmap fm;
		const auto cut_map = enumerate_cuts_eval_funcs(ntk, cut_params.get(), fm);
		logic_ntk cntk(ntk);
		do {
			timeout_occurred = false;
			auto best_area = eval_matches_area_timeout(cntk, cut_map, fm, timeoutfuncs);
			auto area_cover = build_cover(cntk, best_area);
			it_exact_cover_timeout(cntk, area_cover, cut_map, best_area, fm, timeoutfuncs);
			auto lut_ntk = ntk_cover_to_logic_ntk(cntk, area_cover, best_area, fm);
			auto decomp_ntk = logic_ntk_from_luts(lut_ntk, timeoutfuncs, backtrack_limit, tb);
			if (decomp_ntk) {
				auto newsize = decomp_ntk.get_ptr()->nnodes();
				if (newsize < oldsize) {
					cntk = std::move(decomp_ntk.value());
					continue;
				} else {
					break;
				}
			} else {
				cerr << "timeout occurred" << endl;
				timeout_occurred = true;
			}
		} while (timeout_occurred);

		return std::move(cntk);
	}

	optional<xmg> 
		lut_area_timeout_strategy(const xmg& m, const xmg_params* frparams, unsigned lut_size, 
			unsigned timeout, timeout_behavior behavior) {
		xmg cmig(m, frparams);
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;
		vector<tt> timeoutfuncs;
	
		bool timeout_occurred = false;
		auto oldsize = cmig.nnodes();
		funcmap fm;
		const auto cut_map = enumerate_cuts_eval_funcs(cmig, cut_params.get(), fm);
		do {
			timeout_occurred = false;
			auto best_area = eval_matches_area_timeout(cmig, cut_map, fm, timeoutfuncs);
			auto area_cover = build_cover(cmig, best_area);
			it_exact_cover_timeout(cmig, area_cover, cut_map, best_area, fm, timeoutfuncs);
			auto lutxmg = xmg_from_luts(cmig, area_cover, best_area, fm, timeoutfuncs, timeout, behavior);
			if (lutxmg) {
				auto newsize = lutxmg.get_ptr()->nnodes();
				if (newsize < oldsize) {
					cmig = std::move(lutxmg.value());
					continue;
				} else {
					break;
				}
			} else {
				cerr << "timeout occurred" << endl;
				timeout_occurred = true;
			}
		} while (timeout_occurred);

		return std::move(cmig);
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

	inline vector<unsigned> inv(const vector<unsigned>& perm) {
		vector<unsigned> invperm(perm.size());
		for ( auto i = 0u; i < perm.size(); ++i ) { invperm[perm[i]] = i; }
		return invperm;
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
			assert(imap.find(expr[offset]) != imap.end());
			return imap.at(expr[offset]);
		}
	}

	optional<pair<nodeid,bool>> decompose_cut(
			xmg& xmg, const cut* cut, tt& cutfunction, strashmap& shmap, 
			nodemap& nodemap, function_store& fstore, vector<tt>& timeoutfuncs, unsigned timeout, timeout_behavior behavior) {
		const auto& cutnodes = cut->nodes();
		//const auto npn = exact_npn_canonization(cutfunction, phase, perm);
		//auto npn = fstore.npn_canon(cutfunction, phase, perm);
		auto num_vars = tt_num_vars(cutfunction);

        //npn.resize(cutfunction.size());
		
        vector<unsigned> perm; tt phase, npn;
        if (num_vars < 6) {
            npn = exact_npn_canonization(cutfunction, phase, perm);
        } else {
            npn = npn_canonization_lucky(cutfunction, phase, perm);
        }
        npn.resize(cutfunction.size());
        
		auto min_xmg = fstore.min_size_xmg(npn, timeout);
		if (!min_xmg) { // Exact synthesis may have timed out
			if (behavior == rebuild_cover) {
				timeoutfuncs.push_back(cutfunction);
				return boost::none;
			} 
			// Try to resynthesize, starting from the last size that failed
			auto olast_size = fstore.get_last_size(npn);
			assert(olast_size);
			auto start_size = olast_size.get() + 1;
			// We may have stored a heuristic version of this NPN class before,
			// so no need to compute it again
			min_xmg = fstore.heuristic_size_xmg(npn, timeout);
			if (!min_xmg) {
				min_xmg = exact_xmg_expression(npn, timeout, start_size);
				if (!min_xmg) { 
					if (!min_xmg) {
						// Try to compute a heuristic version
						min_xmg = heuristic_xmg_expression(npn, num_vars, timeout, start_size, behavior);
						if (!min_xmg) {
							return boost::none;
						} else {
							fstore.set_heuristic_size_xmg(npn, min_xmg.get(), timeout);
						}
					}
				} else {
					fstore.set_heuristic_size_xmg(npn, min_xmg.get(), timeout);
				}
			} else {
				cout << "Using cached heuristic result" << endl;
			}
		}
		//cout  << "got min: " << min_xmg << endl;
		input_map_t imap;

        auto invperm = inv(perm);
		for (auto i = 0u; i < num_vars; i++) {
			auto inode = nodemap[cutnodes[i]];
			imap['a' + invperm[i]] = make_pair(inode.first, phase.test(i) ^ inode.second);
		}
		const auto majbrackets = find_bracket_pairs(min_xmg.get(), '<', '>');
		const auto xorbrackets = find_bracket_pairs(min_xmg.get(), '[', ']');
		auto inv = false;
		auto offset = 0u;
		if (min_xmg.get()[0] == '!') {
			inv = true;
			offset = 1u;
		}
		auto res = frmaj3_from_string(min_xmg.get(), offset, majbrackets, xorbrackets, imap, xmg, shmap);
		res.second = (res.second != inv);
        res.second = (res.second != (phase.test(num_vars)));
		return res;
	}

	pair<nodeid, bool> decompose_cut(xmg& xmg, const cut* cut, tt& cutfunction, strashmap& shmap,
		nodemap& nodemap, function_store& fstore) {
		vector<tt> emptyset;
		return decompose_cut(xmg, cut, cutfunction, shmap, nodemap, fstore, emptyset, 0, rebuild_cover).get();
	}


	optional<xmg> xmg_from_luts(const xmg& m, const cover& cover,
		const bestmap& best, const funcmap& funcmap, vector<tt>& timeoutfuncs, unsigned timeout) {
		return xmg_from_luts(m, cover, best, funcmap, timeoutfuncs, timeout, rebuild_cover);
	}

	optional<xmg> xmg_from_luts(const xmg& m, const cover& cover,
			const bestmap& best, const funcmap& funcmap, vector<tt>& timeoutfuncs, unsigned timeout, timeout_behavior behavior) {
		bool timeout_occurred = false;

		xmg n;
		xmg_stats stats{
			0u, // Nr. strash hits
			0u, // nr_potentials
			0u, // nr_matches
			0u, // nr_misses
			0u, // nr_undefined
		};

		nodemap nodemap;
		strashmap shmap(m.nnodes() / 2, stats);
		function_store fstore;

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
				nodemap[i] = make_pair(n.create_input(innames[i - 1]), false);
				++progress;
				continue;
			}
			const auto& cut = best.at(i);
			auto& f = *funcmap.at(cut);
			auto decomp_cut = decompose_cut(n, cut, f, shmap, nodemap, fstore, timeoutfuncs, timeout, behavior);
			if (decomp_cut) {
				nodemap[i] = decomp_cut.get();
			} else {
				nodemap[i] = make_pair(0, false);
				timeout_occurred = true;
			}
			cout << "Progress: (" << ++progress << "/" << total_nodes << ")\r";
		}
		cout << endl;
		if (timeout_occurred) {
			cout << "Timeout occurred" << endl;
			return boost::none;
		}
		const auto& outputs = m.outputs();
		const auto& outcompls = m.outcompl();
		const auto& outnames = m.outnames();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto np = nodemap[outputs[i]];
			n.create_output(np.first, outcompls[i] != np.second, outnames[i]);
		}
		optional<xmg> res(move(n));
		return res;
	}

	xmg xmg_from_luts(const xmg& m, const cover& cover, 
			const bestmap& best, const funcmap& funcmap) {
		vector<tt> emptyset;
		return std::move(xmg_from_luts(m, cover, best, funcmap, emptyset, 0).value());
	}

	pair<nodeid, bool> 
		parse_into_logic_ntk(logic_ntk& ntk, const logic_ntk& opt_ntk, const nodemap& nodemap, const input_map_t& imap, bool invert_output) {
		const auto& opt_nodes = opt_ntk.nodes();
		vector<pair<nodeid,bool>> nids(opt_nodes.size());
		vector<nodeid> ntk_fanin(2);
		for (auto i = 0u; i < opt_nodes.size(); i++) {
			const auto& node = opt_nodes[i];
			if (node.pi) {
				nids[i] = imap.at('a' + i);// nodemap.at(cut_fanin[i]).first;
			} else {
				auto fanin1 = nids[node.fanin[0]];
				auto fanin2 = nids[node.fanin[1]];
				ntk_fanin[0] = fanin1.first;
				ntk_fanin[1] = fanin2.first;

				tt localfunc(4, 0);
				tt var0 = fanin1.second ? ~tt_nth_var(0) : tt_nth_var(0);
				tt var1 = fanin2.second ? ~tt_nth_var(1) : tt_nth_var(1);
				for (auto i = 0; i  < 4; i++) {
					localfunc[i] = node.function[var1[i] << 1 | var0[i]];
				}

				if ((i == opt_nodes.size() - 1) && invert_output) {
					nids[i] = make_pair(ntk.create_node(ntk_fanin, ~localfunc), false);
				} else {
					nids[i] = make_pair(ntk.create_node(ntk_fanin, localfunc), false);
				}
			}
		}
		return nids[opt_nodes.size() - 1];
	}

	optional<pair<nodeid, bool>> decompose_cut(logic_ntk& ntk, const ln_node& node, nodemap& nodemap, 
		function_store& fstore, vector<tt>& timeoutfuncs, const unsigned conflict_limit, timeout_behavior behavior) {
		
        vector<unsigned> perm; tt phase, npn;
        if (node.fanin.size() < 6) {
            npn = exact_npn_canonization(node.function, phase, perm);
        } else {
            npn = npn_canonization_lucky(node.function, phase, perm);
        }
        npn.resize(node.function.size());

		synth_spec spec;
		spec.nr_vars = node.fanin.size();
spec.verbose = true;
		//auto opt_ntk = size_optimum_ntk_ns<CMSat::SATSolver>(npn.to_ulong(), &spec);

		auto opt_ntk_str = fstore.size_optimum_ntk_ns(npn, &spec, conflict_limit);
		if (!opt_ntk_str) { // Exact synthesis may have timed out
		//	if (behavior == rebuild_cover) {
				timeoutfuncs.push_back(node.function);
				return boost::none;
		//	} 
			
				/*
			// Try to resynthesize, starting from the last size that failed
			auto olast_size = fstore.get_last_size(npn);
			assert(olast_size);
			auto start_size = olast_size.get() + 1;
			// We may have stored a heuristic version of this NPN class before,
			// so no need to compute it again
			min_xmg = fstore.heuristic_size_xmg(npn, timeout);
			if (!min_xmg) {
				min_xmg = exact_xmg_expression(npn, timeout, start_size);
				if (!min_xmg) { 
					if (!min_xmg) {
						// Try to compute a heuristic version
						min_xmg = heuristic_xmg_expression(npn, num_vars, timeout, start_size, behavior);
						if (!min_xmg) {
							return boost::none;
						} else {
							fstore.set_heuristic_size_xmg(npn, min_xmg.get(), timeout);
						}
					}
				} else {
					fstore.set_heuristic_size_xmg(npn, min_xmg.get(), timeout);
				}
			} else {
				cout << "Using cached heuristic result" << endl;
			}
			*/
		}
		auto opt_ntk = string_to_logic_ntk(opt_ntk_str.get());

		input_map_t imap;

        auto invperm = inv(perm);
		for (auto i = 0u; i < node.fanin.size(); i++) {
			auto inode = nodemap[node.fanin[i]];
			imap['a' + invperm[i]] = make_pair(inode.first, phase.test(i));
		}

		return parse_into_logic_ntk(ntk, opt_ntk, nodemap, imap, phase.test(node.fanin.size()));
	}
	
	logic_ntk logic_ntk_from_luts(const logic_ntk& lut_ntk) {
		vector<tt> timeoutfuncs;
		return logic_ntk_from_luts(lut_ntk, timeoutfuncs, 0, rebuild_cover).value();
	}

	optional<logic_ntk> logic_ntk_from_luts(const logic_ntk& lut_ntk, vector<tt>& timeoutfuncs, 
		const unsigned conflict_limit, timeout_behavior behavior) {
		bool timeout_occurred = false;

		logic_ntk ntk;

		nodemap nodemap;
		function_store fstore;

		const auto& nodes = lut_ntk.nodes();
		const auto total_nodes = lut_ntk.nnodes();
		auto progress = 0u;
		for (auto i = 0u; i < total_nodes; i++) {
			const auto& node = nodes[i];
			if (node.pi) {
				nodemap[i] = make_pair(ntk.create_input(), false);
				++progress;
				continue;
			}
			auto decomp_cut = decompose_cut(ntk, node, nodemap, fstore, timeoutfuncs, conflict_limit, behavior);
			if (decomp_cut) {
				nodemap[i] = decomp_cut.get();
			} else {
				nodemap[i] = make_pair(0, false);
				timeout_occurred = true;
			}
			cout << "Progress: (" << ++progress << "/" << total_nodes << ")\r";
		}
		cout << endl;
		if (timeout_occurred) {
			cout << "Timeout occurred" << endl;
			return boost::none;
		}

		const auto& outputs = lut_ntk.outputs();
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto np = nodemap[outputs[i]];
			ntk.create_output(np.first);
		}

		const auto& innames = lut_ntk.innames();
		for (const auto& name : innames) {
			ntk.add_inname(name);
		}

		const auto& outnames = lut_ntk.outnames();
		for (const auto& name : outnames) {
			ntk.add_outname(name);
		}

		optional<logic_ntk> res(move(ntk));
		return res;
	}



    /*
	// NPN canonization function from ABC
	tt jake_canon(const tt& ttf, unsigned* uCanonPhase, char* pCanonPerm, unsigned num_vars) {
		vector<word> pTruth( ttf.num_blocks() );
		boost::to_block_range( ttf, &pTruth[0] );

		resetPCanonPermArray( pCanonPerm, num_vars );
		*uCanonPhase = luckyCanonicizer_final_fast(&pTruth[0], num_vars, pCanonPerm);

		tt tt_npn( pTruth.begin(), pTruth.end() );
        tt_npn.resize(ttf.size());

		for ( auto i = 0u; i < num_vars; ++i ) {
			pCanonPerm[i] -= 'a';
		}
		Abc_TtImplementNpnConfig(&pTruth[0], num_vars, pCanonPerm, *uCanonPhase);
		tt ttf_new( pTruth.begin(), pTruth.end() );
		ttf_new.resize(ttf.size());
		if (ttf != ttf_new) {
			cout << to_string(ttf) << " != " << to_string(ttf_new) << endl;
			assert(ttf == ttf_new);
		}

		return tt_npn;
	}
    */
}
