#pragma once

#include "cut.h"
#include "lut_cover.h"
#include "function_store.h"
#include <logic_network.h>
#include <vector>
#include <boost/optional.hpp>
#include <convert.h>
#include "npn_canonization.hpp"
#include <maj_io.h>
#include <convert.h>
#include <stdexcept>

namespace majesty {

	class xmg;
	struct xmg_params;

	enum timeout_behavior { 
		rebuild_cover = 0, 
		optimize_heuristically, 
		combine 
	};

	xmg* ptr_lut_area_strategy(const xmg&, unsigned, unsigned);
   	xmg lut_area_strategy(const xmg&, const xmg_params*, unsigned);

	xmg* ptr_lut_area_timeout_strategy(const xmg&, unsigned lut_size, unsigned timeout, unsigned nr_backtracks);
	boost::optional<xmg> lut_area_timeout_strategy(const xmg&, unsigned, unsigned, timeout_behavior);
	boost::optional<xmg> lut_area_timeout_strategy(const xmg&, unsigned, unsigned);
	boost::optional<xmg> lut_area_timeout_strategy(const xmg&, const xmg_params*, unsigned, unsigned, timeout_behavior);

	xmg xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&);
	boost::optional<xmg> xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&, std::vector<cirkit::tt>&, unsigned);
	boost::optional<xmg> xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&, std::vector<cirkit::tt>&, unsigned, timeout_behavior);

	logic_ntk abc_heuristic_logic_ntk(const tt& func) {
		auto optcmd = (boost::format("abc -c \"read_truth %s; strash; resyn2; write_verilog tmp.v\"") % tt_to_hex(func)).str();
		auto success = system(optcmd.c_str());
		if (success != 0) {
			throw std::runtime_error("Heuristic optimization through ABC failed");
		}
		auto upperbound_xmg = read_verilog("tmp.v");
		auto ntk = xmg_to_logic_ntk(upperbound_xmg);
		
		auto xmg_func = simulate_xmg(upperbound_xmg);
		auto ntk_func = ntk.simulate();

		return std::move(ntk);
	}

	static inline std::vector<unsigned> inv(const std::vector<unsigned>& perm) {
		std::vector<unsigned> invperm(perm.size());
		for ( auto i = 0u; i < perm.size(); ++i ) { invperm[perm[i]] = i; }
		return invperm;
	}

	std::pair<nodeid, bool> 
		parse_into_logic_ntk(logic_ntk& ntk, const logic_ntk& opt_ntk, 
				const input_map_t& imap, bool invert_output) {
		const auto& opt_nodes = opt_ntk.nodes();
		std::vector<std::pair<nodeid,bool>> nids(opt_nodes.size());

		for (auto i = 0u; i < opt_nodes.size(); i++) {
			const auto& node = opt_nodes[i];
			if (node.pi) {
				nids[i] = imap.at(i);// nodemap.at(cut_fanin[i]).first;
			} else {
				const auto gate_size = node.fanin.size();
				const auto gate_tt_size = (1u << gate_size);
				std::vector<nodeid> ntk_fanin(gate_size);
				for (auto j = 0u; j < gate_size; j++) {
					ntk_fanin[j] = nids[node.fanin[j]].first;
				}
				tt localfunc(gate_tt_size, 0);
				for (auto j = 0u; j  < gate_tt_size; j++) {
					auto func_idx = 0u;
					for (auto k = 0u; k < gate_size; k++) {
						auto tbit = tt_nth_var(k).test(j);
						if (nids[node.fanin[k]].second) {
							tbit = !tbit;
						}
						func_idx += (tbit << k);
					}
					localfunc[j] = node.function[func_idx];
				}

				if ((i == opt_nodes.size() - 1) && invert_output) {
					nids[i] = std::make_pair(ntk.create_node(ntk_fanin, ~localfunc), false);
				} else {
					nids[i] = std::make_pair(ntk.create_node(ntk_fanin, localfunc), false);
				}
			}
		}
		return nids[opt_nodes.size() - 1];
	}

	template<typename S>
	boost::optional<std::pair<nodeid, bool>> decompose_cut(logic_ntk& ntk, const ln_node& node, nodemap& nodemap, 
		function_store& fstore, std::vector<tt>& timeoutfuncs, const unsigned conflict_limit, timeout_behavior behavior) {
		
		std::vector<unsigned> perm; tt phase, npn;
        if (node.fanin.size() < 6) {
            npn = exact_npn_canonization(node.function, phase, perm);
        } else {
            npn = npn_canonization_lucky(node.function, phase, perm);
        }
        npn.resize(node.function.size());

		synth_spec spec;
		spec.verbose = true;
		spec.nr_vars = node.fanin.size();
		if (spec.nr_vars >= 3) {
			// We have to manually set the gate size and gate_tt_size here,
			// as exact synthesis may not be called when we use stored results.
			// ES normally deduces these values, but if it's not called we
			// still need them to parse the results back into the logic_ntk. 
			spec.gate_size = 3;
			spec.gate_tt_size = 7;
		}

		logic_ntk opt_ntk;
		auto entry = fstore.get_entry(npn);
		if (!entry) { // This function hasn't been synthesized yet
			auto synth_ntk = size_optimum_ntk_ns<S>(npn, &spec, conflict_limit);
			if (synth_ntk) {
				opt_ntk = std::move(synth_ntk.get());
				fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), true, conflict_limit);
			} else { // Conflict limit was breached in exact synthesis
				if (behavior == rebuild_cover) {
					timeoutfuncs.push_back(node.function);
					return boost::none;
				}
				const auto last_size = spec.nr_gates;
				lbool exists = l_False;
				auto invert = npn.test(0);
				// Try to find a heuristic solution by increasing the nr of gate, starting
				// from the last attempted number.
				init_solver<S>(conflict_limit);
				for (auto nr_gates = last_size + 1; nr_gates < last_size + 4; nr_gates++) {
					restart_solver<S>();
					spec.nr_gates = nr_gates;
					if (spec.verbose) {
						std::cout << "trying to get heuristicic result with " 
							<< nr_gates << " gates" << std::endl;
					}
					exists = exists_fanin_3_ntk_ns<S>(invert ? ~npn : npn, &spec);
					if (exists == l_True) {
						opt_ntk = extract_fanin_3_ntk_ns<S>(&spec, invert);
						break;
					}
				}
				destroy_solver<S>();
				if (exists == l_True) { 
					fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), false, conflict_limit);
				} else { // Use ABC to find a heuristic decomposition
					opt_ntk = abc_heuristic_logic_ntk(npn);
					fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), false, conflict_limit);
				}
			}
		} else {
			opt_ntk = string_to_logic_ntk(entry.get().expr);
		}

		input_map_t imap;
        auto invperm = inv(perm);
		for (auto i = 0u; i < node.fanin.size(); i++) {
			auto inode = nodemap[node.fanin[i]];
			imap[invperm[i]] = std::make_pair(inode.first, phase.test(i));
		}

		return parse_into_logic_ntk(ntk, opt_ntk, imap, phase.test(node.fanin.size()));
	}

	template<typename S>
	boost::optional<logic_ntk> logic_ntk_from_luts(const logic_ntk& lut_ntk, std::vector<tt>& timeoutfuncs, 
		const unsigned conflict_limit, timeout_behavior behavior) {
		bool timeout_occurred = false;

		logic_ntk ntk;
		
		nodemap nodemap;
		function_store fstore;

		const std::vector<nodeid> emptyfanin;
		bool have_const1 = false, have_const0 = false;
		nodeid const1id = 0, const0id = 0;

		const auto& nodes = lut_ntk.nodes();
		const auto total_nodes = lut_ntk.nnodes();
		auto progress = 0u;
		for (auto i = 0u; i < total_nodes; i++) {
			const auto& node = nodes[i];
			if (node.pi) {
				nodemap[i] = std::make_pair(ntk.create_input(), false);
				++progress;
				continue;
			} else if (node.fanin.size() == 0) {
				// Const 0 or 1
				if (node.function == tt_const1()) {
					if (!have_const1) {
						have_const1 = true;
						const1id = ntk.create_node(emptyfanin, tt(1, 1));
						nodemap[i] = std::make_pair(const1id, false);
					} else {
						nodemap[i] = std::make_pair(const1id, false);
					}
				} else {
					if (!have_const0) {
						have_const0 = true;
						const0id = ntk.create_node(emptyfanin, tt(1, 0));
						nodemap[i] = std::make_pair(const0id, false);
					} else {
						nodemap[i] = std::make_pair(const0id, false);
					}
				}
				continue;
			} else if (node.fanin.size() == 1) {
				// Primary input
				nodemap[i] = nodemap[node.fanin[0]];
				continue;
			}
			auto decomp_cut = decompose_cut<S>(ntk, node, nodemap, fstore, timeoutfuncs, conflict_limit, behavior);
			if (decomp_cut) {
				nodemap[i] = decomp_cut.get();
			} else {
				nodemap[i] = std::make_pair(0, false);
				timeout_occurred = true;
			}
			std::cout << "Progress: (" << ++progress << "/" << total_nodes << ")\r";
		}

		std::cout << std::endl;
		if (timeout_occurred) {
			std::cout << "Timeout occurred" << std::endl;
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

		boost::optional<logic_ntk> res(std::move(ntk));
		return res;
	}

	void mine_functions(const xmg&, const xmg_params*, unsigned);
	cirkit::tt jake_canon(const cirkit::tt&, unsigned* uCanonPhase, char* pCanonPerm, unsigned num_vars);

	template<typename S>
	logic_ntk lut_area_strategy(const logic_ntk& ntk, unsigned lut_size, unsigned conflict_limit = 0, timeout_behavior tb = optimize_heuristically) {
		auto cut_params = default_cut_params();
		cut_params->klut_size = lut_size;
		std::vector<tt> timeoutfuncs;
	
		logic_ntk cntk(ntk);
		auto ctu = false;
		do {
			ctu = false;
			funcmap fm;
			auto oldsize = cntk.nnodes();
			const auto cut_map = enumerate_cuts_eval_funcs(cntk, cut_params.get(), fm);
			auto best_area = eval_matches_area_timeout(cntk, cut_map, fm, timeoutfuncs);
			auto area_cover = build_cover(cntk, best_area);
			it_exact_cover_timeout(cntk, area_cover, cut_map, best_area, fm, timeoutfuncs);
			auto lut_ntk = ntk_cover_to_logic_ntk(cntk, area_cover, best_area, fm);
			auto decomp_ntk = logic_ntk_from_luts<S>(lut_ntk, timeoutfuncs, conflict_limit, tb);
			if (decomp_ntk) {
				auto newsize = decomp_ntk.get_ptr()->nnodes();
				if (newsize < oldsize) {
					std::cout << "oldsize: " << oldsize << std::endl;
					std::cout << "newsize: " << newsize << std::endl;
					std::cout << "continuing" << std::endl;
					cntk = std::move(decomp_ntk.value());
					ctu = true;
				} else {
					std::cout << "oldsize: " << oldsize << std::endl;
					std::cout << "newsize: " << newsize << std::endl;
					std::cout << "not continuing" << std::endl;
				}
			} else {
				std::cerr << "timeout occurred" << std::endl;
				ctu = true;
			}
		} while (ctu);

		return cntk;
	}

	int recursive_deselect(const nodeid nid, const std::vector<ln_node>& nodes, std::unordered_map<nodeid,unsigned>& nref) {
		const auto& node = nodes[nid];
		if (node.pi) {
			return 0;
		}
		auto area = 1;
		for (auto nodeid : node.fanin) {
			nref[nodeid] -= 1;
			if (nref[nodeid] == 0) {
				area += recursive_deselect(nodeid, nodes, nref);
			}
		}
		return area;
	}

	int recursive_select(const nodeid nid, const std::vector<ln_node>& nodes, std::unordered_map<nodeid,unsigned>& nref) {
		const auto& node = nodes[nid];
		if (node.pi) {
			return 0;
		}
		auto area = 1;
		const auto& fanin = node.fanin;
		for (auto nodeid : fanin) {
			nref[nodeid] += 1;
			if (nref[nodeid] == 1) {
				area += recursive_select(nodeid, nodes, nref);
			}
		}
		return area;
	}

	int virtual_recursive_deselect(const logic_ntk& ntk, const logic_ntk& opt_ntk, 
			const std::vector<nodeid>& fanin, const std::unordered_map<nodeid, nodeid>& nodemap, 
			std::unordered_map<nodeid, unsigned>& nref) {
		auto area = 0;
		const auto& opt_ntk_nodes = opt_ntk.nodes();
		const auto& nr_opt_ntk_nodes = opt_ntk.nnodes();

		for (auto i = 0u; i < nr_opt_ntk_nodes; i++) {
			const auto& node = opt_ntk_nodes[i];
			if (node.pi) {
				auto faninid = nodemap.at(fanin[i]); 
				nref[faninid] -= 1;
				if (nref.at(faninid) == 0) {
					area += recursive_deselect(faninid, ntk.nodes(), nref);
				}
			} else {
				++area;
			}
		}

		return area;
	}

	int virtual_recursive_select(const logic_ntk& ntk, const logic_ntk& opt_ntk, 
			const std::vector<nodeid>& fanin, const std::unordered_map<nodeid,nodeid>& nodemap, 
			std::unordered_map<nodeid,unsigned>& nref) {
		auto area = 0;
		const auto& opt_ntk_nodes = opt_ntk.nodes();
		const auto& nr_opt_ntk_nodes = opt_ntk.nnodes();

		for (auto i = 0u; i < nr_opt_ntk_nodes; i++) {
			const auto& node = opt_ntk_nodes[i];
			if (node.pi) {
				assert(nodemap.find(fanin[i]) != nodemap.end());
				auto faninid = nodemap.at(fanin[i]); 
				nref[faninid] += 1;
				if (nref.at(faninid) == 1) {
					area += recursive_select(faninid, ntk.nodes(), nref);
				}
			} else {
				++area;
			}
		}

		return area;
	}

	inline nodeid select_opt_ntk(logic_ntk& ntk, const logic_ntk& opt_ntk, const input_map_t& imap, 
			bool invert_output, std::unordered_map<nodeid,unsigned>& nref) {
		const auto& opt_ntk_nodes = opt_ntk.nodes();
		const auto& nr_opt_ntk_nodes = opt_ntk.nnodes();

		std::vector<std::pair<nodeid,bool>> nids(nr_opt_ntk_nodes);

		for (auto i = 0u; i < nr_opt_ntk_nodes; i++) {
			const auto& node = opt_ntk_nodes[i];
			if (node.pi) {
				nids[i] = imap.at(i);
				const auto faninid = nids[i].first;
				nref[faninid] += 1;
				if (nref[faninid] == 1) {
					recursive_select(faninid, ntk.nodes(), nref);
				}
			} else {
				const auto gate_size = node.fanin.size();
				const auto gate_tt_size = (1u << gate_size);
				std::vector<nodeid> virtfanin;
				for (auto j = 0u; j < gate_size; j++) {
					virtfanin.push_back(nids[node.fanin[j]].first);
				}
				tt localfunc(gate_tt_size, 0);
				for (auto j = 0u; j < gate_tt_size; j++) {
					auto func_idx = 0u;
					for (auto k = 0u; k < gate_size; k++) {
						auto tbit = tt_nth_var(k).test(j);
						if (nids[node.fanin[k]].second) {
							tbit = !tbit;
						}
						func_idx += (tbit << k);
					}
					localfunc[j] = node.function[func_idx];
				}
				if ((i == nr_opt_ntk_nodes - 1) && invert_output) {
					localfunc = ~localfunc;
				}
				auto new_nodeid = ntk.create_node(virtfanin, localfunc);
				nids[i] = std::make_pair(new_nodeid, false);
				nref[new_nodeid] = 1;
			}
		}

		return nids[nr_opt_ntk_nodes - 1].first;
	}

	template<typename S>
	logic_ntk logic_ntk_from_cuts(const logic_ntk& cut_ntk, const cutmap& cut_map, unsigned conflict_limit) {
		logic_ntk tmp_ntk;
		
		std::unordered_map<nodeid,nodeid> nodemap;
		std::unordered_map<nodeid,unsigned> nref;
		function_store fstore;

		const auto& nodes = cut_ntk.nodes();
		const auto total_nodes = cut_ntk.nnodes();

		funcmap fm = compute_all_functions(cut_ntk, cut_map);

		auto progress = 0u;
		for (auto i = 0u; i < total_nodes; i++) {
			const auto& node = nodes[i];
			if (node.pi) {
				nodemap[i] = tmp_ntk.create_input();
				++progress;
				continue;
			}
			const auto& node_cuts = cut_map.at(i);
			auto smallest_add = std::numeric_limits<unsigned>::max();
			cut* best_cut = nullptr;
			auto found_const_cut = false;

			for (const auto& cut : node_cuts) {
				const auto& cutfunc = *fm.at(cut.get());
				if (cut->size() == 1 && cut->nodes()[0] == i) { // Trivial cut
					continue;
				} else if (cut->size() == 1) { // E.g. may map to a PI.
					nodemap[i] = nodemap[cut->nodes()[0]];
					continue;
				} else if (cut->size() == 0) { // Const 1 or 0
					found_const_cut = true;
					if (cutfunc == tt_const0()) {
						nodemap[i] = tmp_ntk.get_const0_node();
					} else {
						nodemap[i] = tmp_ntk.get_const1_node();
					}
					break;
				} else {
					std::vector<unsigned> perm; tt phase, npn;
					const auto& cutnodes = cut->nodes();
					if (cutnodes.size() < 6) {
						npn = exact_npn_canonization(cutfunc, phase, perm);
					} else {
						npn = npn_canonization_lucky(cutfunc, phase, perm);
					}
					npn.resize(cutfunc.size());

					synth_spec spec;
					spec.verbose = true;
					spec.nr_vars = cut->size();
					if (spec.nr_vars >= 3) {
						// We have to manually set the gate size and gate_tt_size here,
						// as exact synthesis may not be called when we use stored results.
						// ES normally deduces these values, but if it's not called we
						// still need them to parse the results back into the logic_ntk. 
						spec.gate_size = 3;
						spec.gate_tt_size = 7;
					}

					logic_ntk opt_ntk;
					auto entry = fstore.get_entry(npn);
					if (!entry) { // This function hasn't been synthesized yet
						auto synth_ntk = size_optimum_ntk_ns<S>(npn, &spec, conflict_limit);
						if (synth_ntk) {
							opt_ntk = std::move(synth_ntk.get());
							fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), true, conflict_limit);
						} else { // Conflict limit was breached in exact synthesis
							const auto last_size = spec.nr_gates;
							lbool exists = l_False;
							auto invert = npn.test(0);
							// Try to find a heuristic solution by increasing the nr of gate, starting
							// from the last attempted number.
							init_solver<S>(conflict_limit);
							for (auto nr_gates = last_size + 1; nr_gates < last_size + 4; nr_gates++) {
								restart_solver<S>();
								spec.nr_gates = nr_gates;
								if (spec.verbose) {
									std::cout << "trying to get heuristicic result with " 
										<< nr_gates << " gates" << std::endl;
								}
								exists = exists_fanin_3_ntk_ns<S>(invert ? ~npn : npn, &spec);
								if (exists == l_True) {
									opt_ntk = extract_fanin_3_ntk_ns<S>(&spec, invert);
									break;
								}
							}
							destroy_solver<S>();
							if (exists == l_True) { 
								fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), false, conflict_limit);
							} else { // Use ABC to find a heuristic decomposition
								opt_ntk = abc_heuristic_logic_ntk(npn);
								fstore.set_entry(npn, logic_ntk_to_string(opt_ntk), false, conflict_limit);
							}
						}
					} else {
						opt_ntk = string_to_logic_ntk(entry.get().expr);
					}

					auto virt_nodes_added = virtual_recursive_select(tmp_ntk, opt_ntk, cut->nodes(), nodemap, nref);
					if (virt_nodes_added < smallest_add) {
						best_cut = cut.get();
						smallest_add = virt_nodes_added;
					}
					auto virt_nodes_saved = virtual_recursive_deselect(tmp_ntk, opt_ntk, cut->nodes(), nodemap, nref);
					assert(virt_nodes_added == virt_nodes_saved);
				}
			}
			if (found_const_cut) {
				continue;
			}
			assert(best_cut != nullptr);
			const auto& cutfunc = *fm.at(best_cut);
			std::vector<unsigned> perm; tt phase, npn;
			const auto& cutnodes = best_cut->nodes();
			if (cutnodes.size() < 6) {
				npn = exact_npn_canonization(cutfunc, phase, perm);
			} else {
				npn = npn_canonization_lucky(cutfunc, phase, perm);
			}
			npn.resize(cutfunc.size());

			// NOTE: optimized network must exist in function store at this point
			auto entry = fstore.get_entry(npn);
			auto opt_ntk = string_to_logic_ntk(entry.get().expr);

			input_map_t imap;
			auto invperm = inv(perm);
			for (auto i = 0u; i < cutnodes.size(); i++) {
				const auto inode = nodemap[cutnodes[i]];
				imap[invperm[i]] = std::make_pair(inode, phase.test(i));
			}

			nodemap[i] = select_opt_ntk(tmp_ntk, opt_ntk, imap, phase.test(cutnodes.size()), nref);
			std::cout << "Progress: (" << ++progress << "/" << total_nodes << ")\r";
		}
		std::cout << std::endl;

		nref.clear();
		{
			const auto& outputs = cut_ntk.outputs();
			for (auto i = 0u; i < outputs.size(); i++) {
				const auto nid = nodemap[outputs[i]];
				tmp_ntk.create_output(nid);
				nref[nid] = 1;
				recursive_select(nid, tmp_ntk.nodes(), nref);
			}
		}

		logic_ntk ntk;
		auto tmp_nodes = tmp_ntk.nodes();
		nodemap.clear();
		for (auto i = 0u; i < tmp_ntk.nnodes(); i++) {
			const auto& node = tmp_nodes[i];
			if (node.pi) {
				nodemap[i] = ntk.create_input();
				++progress;
			} else {
				auto refs = nref[i];
				if (refs > 0) {
					std::vector<nodeid> nfanin;
					for (auto id : node.fanin) {
						nfanin.push_back(nodemap[id]);
					}
					nodemap[i] = ntk.create_node(nfanin, node.function);
				}
			}
		}
		
		const auto& tmp_outputs = tmp_ntk.outputs();
		for (auto i = 0u; i < tmp_outputs.size(); i++) {
			const auto nid = nodemap[tmp_outputs[i]];
			ntk.create_output(nid);
		}

		const auto& innames = cut_ntk.innames();
		for (const auto& name : innames) {
			ntk.add_inname(name);
		}

		const auto& outnames = cut_ntk.outnames();
		for (const auto& name : outnames) {
			ntk.add_outname(name);
		}

		return ntk;
	}

	template<typename S>
	logic_ntk size_rewrite_strategy(const logic_ntk& ntk, unsigned cut_size, unsigned conflict_limit) {
		auto cut_params = default_cut_params();
		cut_params->klut_size = cut_size;
	
		logic_ntk cntk(ntk);
		auto ctu = false;
		do {
			ctu = false;
			funcmap fm;
			auto oldsize = cntk.nnodes();
			const auto cut_map = enumerate_cuts_eval_funcs(cntk, cut_params.get(), fm);
			auto decomp_ntk = logic_ntk_from_cuts<S>(cntk, cut_map, conflict_limit);
			auto newsize = decomp_ntk.nnodes();
			std::cout << "oldsize: " << oldsize << std::endl;
			std::cout << "newsize: " << newsize << std::endl;
			if (newsize < oldsize) {
				std::cout << "continuing" << std::endl;
				cntk = std::move(decomp_ntk);
				ctu = true;
			} else {
				std::cout << "not continuing" << std::endl;
			}
		} while (ctu);

		return std::move(cntk);
	}

}
