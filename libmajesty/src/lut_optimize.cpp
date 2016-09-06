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
using boost::optional;

namespace majesty {

	xmg* ptr_lut_area_strategy(const xmg& m, unsigned lut_size, unsigned nr_backtracks) {
		auto frparams = default_xmg_params();
		frparams->nr_backtracks = nr_backtracks;
		return new xmg(lut_area_strategy(m, frparams.get(), lut_size));
	}

	xmg lut_area_strategy(const xmg& m, const xmg_params* frparams, unsigned lut_size) {
        return lut_area_timeout_strategy(m, frparams, lut_size, 0, 0).value();
	}

	xmg* ptr_lut_area_timeout_strategy(const xmg& m, unsigned lut_size, unsigned timeout, unsigned effort, unsigned nr_backtracks) {
		auto frparams = default_xmg_params();
		frparams->nr_backtracks = nr_backtracks;
		auto oxmg = lut_area_timeout_strategy(m, frparams.get(), lut_size, timeout, effort);
		if (oxmg) {
			return new xmg(std::move(oxmg.get()));
		} else {
			return NULL;
		}
	}

	optional<xmg> lut_area_timeout_strategy(const xmg& m, const xmg_params* frparams, unsigned lut_size, unsigned timeout, unsigned effort) {
		xmg cmig(m, frparams);
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;
		vector<tt> timeoutfuncs;
	
		unsigned count = 0;
		bool timeout_occurred = false;
		do {
			auto oldsize = cmig.nnodes();
			funcmap fm;
			const auto cut_map = filtered_enumerate_cuts(cmig, cut_params.get(), fm, timeoutfuncs);
			auto best_area = eval_matches_area(cmig, cut_map);
			auto area_cover = build_cover(cmig, best_area);
			it_exact_cover(cmig, area_cover, cut_map, best_area);
			auto lutxmg = xmg_from_luts(cmig, area_cover, best_area, fm, timeoutfuncs, timeout);
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
			++count;
			if (effort > 0 && count > effort) {
				cerr << "Timeout budget spent" << endl;
				return boost::none;
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

	inline void inv(const char* perm, char* invperm, unsigned num_vars) {
		for ( auto i = 0u; i < num_vars; ++i ) { invperm[(int)perm[i]] = i; }
	}

    // Maps variable indices between ABC's and Cirkit's internal representations
    static inline unsigned inv_var_idx(unsigned idx, unsigned num_vars) {
        return num_vars - idx - 1u;
    }

    static inline unsigned inv_phase(unsigned phase, unsigned num_vars) {
       unsigned res = 0u;
       for (auto i = 0u; i < num_vars; i++) {
           auto mask = (phase >> (num_vars - i - 1u)) & 1;
           res |= mask << i;
       }
       res |= (((phase >> num_vars) & 1u) << num_vars);
       return res;
    }


	static inline void jake_inv(const char* perm, char* invperm, unsigned num_vars) {
		for ( auto i = 0u; i < num_vars; ++i ) { invperm[(unsigned)perm[inv_var_idx(i, num_vars)]] = i; }
	}

    inline vector<unsigned> inv(const vector<unsigned>& perm) {
        vector<unsigned> invperm(perm.size());
		for ( auto i = 0u; i < perm.size(); ++i ) { invperm[perm[i]] = i; }
        return invperm;
	}

	pair<nodeid,bool> frmaj3_from_string(const string& expr, 
			unsigned offset, const bracket_map_t& majbrackets, 
			const bracket_map_t& xorbrackets,
			const input_map_t& imap, xmg& xmg, strashmap& shmap) {
        //cout << "expr: " << expr << endl;
        //cout << "offset: " << offset << endl;
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
			nodemap& nodemap, function_store& fstore, vector<tt>& timeoutfuncs, unsigned timeout) {
		const auto& cutnodes = cut->nodes();
		//const auto npn = exact_npn_canonization(cutfunction, phase, perm);
		//auto npn = fstore.npn_canon(cutfunction, phase, perm);
		auto num_vars = tt_num_vars(cutfunction);

		unsigned uCanonPhase; char pCanonPerm[num_vars]; char invperm[num_vars];
		auto npn = jake_canon(cutfunction, &uCanonPhase, pCanonPerm, num_vars);
        auto invphase = inv_phase(uCanonPhase, num_vars);
        //npn.resize(cutfunction.size());
		
        vector<unsigned> perm; tt phase;
		//const auto npn = exact_npn_canonization(cutfunction, phase, perm);
		const auto enpn = exact_npn_canonization(cutfunction, phase, perm);
        assert(npn == enpn);
        
		//cout  << "got npn: " << to_string(npn) << endl;
		const auto min_xmg = fstore.min_size_xmg(npn, timeout);
		const auto emin_xmg = fstore.min_size_xmg(enpn, timeout);
		if (!min_xmg) { // Exact synthesis may have timed out
			timeoutfuncs.push_back(cutfunction);
			return boost::none;
		}
		//cout  << "got min: " << min_xmg << endl;
		input_map_t imap;
		inv(pCanonPerm, invperm, num_vars);
		//jake_inv(pCanonPerm, invperm, num_vars);

        for (auto i = 0u; i < num_vars; i++) {
            if(perm[i] != pCanonPerm[i]) {
                cout << endl;
                for (auto j = 0u; j < num_vars; j++) {
                    cout << "perm[" << j << "] = " << perm[j] << ", ";
                    cout << "pCanonPerm[" << j << "] = " << (int)pCanonPerm[j] << endl;
                }
                for (auto j = 0u; j < num_vars; j++) {
                    //cout << "invperm[" << j << "] = " << (int)invperm[j] << endl;
                }
                for (auto j = 0u; j <= num_vars; j++) {
                    cout << "phase[" << j << "] = " << phase.test(j) << ", ";
                    cout << "uCanonPhase[" << j << "] = " << ((uCanonPhase >> j) & 1u) << endl;
                }
                for (auto j = 0u; j <= num_vars; j++) {
                    //cout << "invphase[" << j << "] = " << ((invphase >> j) & 1u) << endl;
                }
            }
            cout << "tt:        " << cutfunction << endl;
            cout << "npn:       " << to_string(npn) << endl;
            cout << "enpn:      " << to_string(enpn) << endl;
            cout << "expr:      " << min_xmg.get() << endl;
            cout << "enpn_expr: " << emin_xmg.get() << endl;
            assert(0);
        }
        if ((uCanonPhase >> num_vars) & 1) {
            cout << "Got inverted output!" << endl;
        }

        //auto invperm = inv(perm);
		for (auto i = 0u; i < num_vars; i++) {
			auto inode = nodemap[cutnodes[i]];
			imap['a' + invperm[i]] = make_pair(
					//inode.first, ((invphase >> i) & 1u) ^ inode.second);
					inode.first, ((uCanonPhase >> i) & 1u) ^ inode.second);
                    //inode.first, phase.test(i) ^ inode.second);
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
		//res.second = (res.second != ((invphase >> num_vars) & 1u ));
		res.second = (res.second != ((uCanonPhase >> num_vars) & 1u ));
        //res.second = (res.second != (phase.test(num_vars)));
        //cout << "nurpie!" << endl;
		return res;
	}
    // p[0] = 1, p[1] = 0 --> ip[0] = 1, p[1] = 0

	pair<nodeid, bool> decompose_cut(xmg& xmg, const cut* cut, tt& cutfunction, strashmap& shmap,
		nodemap& nodemap, function_store& fstore) {
		vector<tt> emptyset;
		return decompose_cut(xmg, cut, cutfunction, shmap, nodemap, fstore, emptyset, 0).get();
	}

	optional<xmg> xmg_from_luts(const xmg& m, const cover& cover,
			const bestmap& best, const funcmap& funcmap, vector<tt>& timeoutfuncs, unsigned timeout) {
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
			auto decomp_cut = decompose_cut(n, cut, f, shmap, nodemap, fstore, timeoutfuncs, timeout);
			if (decomp_cut) {
				nodemap[i] = decomp_cut.get();
			} else {
				//nodemap[i] = make_pair(0, false);
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
        /*
		Abc_TtImplementNpnConfig(&pTruth[0], num_vars, pCanonPerm, *uCanonPhase);
		tt ttf_new( pTruth.begin(), pTruth.end() );
		ttf_new.resize(ttf.size());
		if (ttf != ttf_new) {
			cout << to_string(ttf) << " != " << to_string(ttf_new) << endl;
			assert(ttf == ttf_new);
		}
        */

		return tt_npn;
	}
}
