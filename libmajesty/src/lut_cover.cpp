#include "lut_cover.h"
#include "cut.h"
#include "truth_table_utils.hpp"
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;
using namespace cirkit;

namespace majesty {

	static inline unsigned int 
		cut_arrival_time(const cut* c, const nintmap& arrival_times) {
			auto max_node_arrival = 0u;
			auto nodes = c->nodes();
			for (auto node : nodes) {
				auto atime = arrival_times.at(node);
				if (atime > max_node_arrival) {
					max_node_arrival = atime;
				}
			}
			return max_node_arrival+1;
		}

	static inline unsigned int 
		cut_area_flow(const cut* c, const nintmap& area_flows) {
			auto res = 1;
			auto cutnodes = c->nodes();
			for (auto node : cutnodes) {
				res += area_flows.at(node);
			}
			return res;
		}

	static inline const tuple<cut*,unsigned> 
		best_cut_area(const vector<cut*>& cuts, const nintmap& area_flows) {
			cut* best_cut = NULL;
			auto best_aflow = numeric_limits<unsigned>::max();
			auto best_size = numeric_limits<unsigned>::max();
			for (const auto& cut : cuts) {
				if (cut->size() == 1) { // Trivial cut
					continue;
				}
				auto aflow = cut_area_flow(cut, area_flows);
				if ((aflow < best_aflow) ||
						(aflow == best_aflow && cut->size() < best_size)) {
					best_cut = cut;
					best_size = cut->size();
					best_aflow = aflow;
				}
			}
			return make_tuple(best_cut, best_aflow);
		}

	static inline const tuple<cut*,unsigned,unsigned> 
		best_cut_area_recover(const vector<cut*>& cuts, const nintmap& area_flows,
				const nintmap& atimes, unsigned int req) {
			cut* best_cut = NULL;
			auto best_atime = numeric_limits<unsigned>::max();
			auto best_size = numeric_limits<unsigned>::max();
			auto best_aflow = numeric_limits<unsigned>::max();
			for (const auto& cut : cuts) {
				if (cut->size() == 1) { // Trivial cut
					continue;
				}
				auto atime = cut_arrival_time(cut, atimes);
				auto aflow = cut_area_flow(cut, area_flows);
				if ((atime <= req && aflow < best_aflow) ||
						(atime <= req && aflow == best_aflow && cut->size() < best_size)) {
					best_cut = cut;
					best_atime = atime;
					best_size = cut->size();
					best_aflow = aflow;
				}
			}
			return make_tuple(best_cut, best_aflow, best_atime);
		}

	static inline const tuple<cut*,unsigned,unsigned> 
		best_cut_depth(const vector<cut*>& cuts, const nintmap& arrival_times, 
				const nintmap& area_flows) {
			cut* best_cut = NULL;
			auto best_atime = numeric_limits<unsigned>::max();
			auto best_size = numeric_limits<unsigned>::max();
			auto best_aflow = numeric_limits<unsigned>::max();
			for (const auto& cut : cuts) {
				if (cut->size() == 1) { // Trivial cut
					continue;
				}
				auto atime = cut_arrival_time(cut, arrival_times);
				auto aflow = cut_area_flow(cut, area_flows);
				if ((atime < best_atime) ||
						(atime == best_atime && cut->size() < best_size) || 
						(atime == best_atime && cut->size() == best_size &&
						 aflow < best_aflow)) {
					best_cut = cut;
					best_atime = atime;
					best_size = cut->size();
					best_aflow = aflow;
				}
			}
			return make_tuple(best_cut, best_atime, best_aflow);
		}

	bestmap eval_matches_depth(const xmg& m, const cutmap& cut_map, 
				nintmap& atimes) {
		cout << "Evaluating matches depth..." << endl;
		bestmap best(m.nnodes());
		nintmap area_flows(m.nnodes());
		auto nodes = m.nodes();
		vector<cut*> eval_cuts;
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes.at(i);
			if (is_pi(node)) {
				best[i] = cut_map.at(i)[0].get();
				atimes[i] = 0;
				area_flows[i] = 0;
			} else if (node.ecrep == static_cast<nodeid>(i)) {
				// Eq. class representative
				const auto& cuts = cut_map.at(i);
				eval_cuts.clear();
				for (const auto& cut : cuts) {
					eval_cuts.push_back(cut.get());
				}
				auto depth_eval = best_cut_depth(eval_cuts, atimes, area_flows);
				best[i] = get<0>(depth_eval);
				atimes[i] = get<1>(depth_eval);
				area_flows[i] = get<2>(depth_eval);
			}
		}
		return best;
	}

	bestmap eval_matches_area(const xmg& m, const cutmap& cut_map) {
		cout << "Evaluating matches area..." << endl;
		bestmap best(m.nnodes());
		nintmap area_flows(m.nnodes());
		auto nodes = m.nodes();
		vector<cut*> eval_cuts;
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes.at(i);
			if (is_pi(node)) {
				best[i] = cut_map.at(i)[0].get();
				area_flows[i] = 0;
			} else if (node.ecrep == static_cast<nodeid>(i)) {
				// Eq. class representative
				const auto& cuts = cut_map.at(i);
				eval_cuts.clear();
				for (const auto& cut : cuts) {
					eval_cuts.push_back(cut.get());
				}
				auto area_eval = best_cut_area(eval_cuts, area_flows);
				best[i] = get<0>(area_eval);
				area_flows[i] = get<1>(area_eval);
			}
		}

		return best;
	}

	bestmap eval_matches_area_timeout(const xmg& m, const cutmap& cut_map, const funcmap& fm, const vector<tt>& timeoutfuncs) {
		cout << "Evaluating matches area..." << endl;
		bestmap best(m.nnodes());
		nintmap area_flows(m.nnodes());
		auto nodes = m.nodes();
		vector<cut*> eval_cuts;
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes.at(i);
			if (is_pi(node)) {
				best[i] = cut_map.at(i)[0].get();
				area_flows[i] = 0;
			} else if (node.ecrep == static_cast<nodeid>(i)) {
				// Eq. class representative
				const auto& cuts = cut_map.at(i);
				eval_cuts.clear();
				for (const auto& cut : cuts) {
					const auto& cutfunction = *fm.at(cut.get());
					if (find(timeoutfuncs.begin(), timeoutfuncs.end(), cutfunction) == timeoutfuncs.end()) {
						eval_cuts.push_back(cut.get());
					}
				}
				auto area_eval = best_cut_area(eval_cuts, area_flows);
				best[i] = get<0>(area_eval);
				area_flows[i] = get<1>(area_eval);
			}
		}

		return best;
	}

	bestmap eval_matches_af_recover(const xmg& m, const cutmap& cut_map, 
			const reqmap& req, const nintmap& nref) {
		cout << "Evaluating matches area flow (recover)..." << endl;
		bestmap best(m.nnodes());
		nintmap arrival_times(m.nnodes());
		nintmap area_flows(m.nnodes());
		auto nodes = m.nodes();
		vector<cut*> eval_cuts;
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes.at(i);
			if (is_pi(node)) {
				best[i] = cut_map.at(i)[0].get();
				arrival_times[i] = 0;
				area_flows[i] = 0;
			} else if (node.ecrep == static_cast<nodeid>(i)) {
				const auto& cuts = cut_map.at(i);
				eval_cuts.clear();
				for (const auto& cut : cuts) {
					eval_cuts.push_back(cut.get());
				}
				const auto eval_area_rec = best_cut_area_recover(eval_cuts, area_flows, arrival_times, req.at(i)); 
				auto numfanouts = nref.at(i) == 0 ? 1 : nref.at(i);
				best[i] = get<0>(eval_area_rec);
				area_flows[i] = get<1>(eval_area_rec) / numfanouts;
				arrival_times[i] = get<2>(eval_area_rec);
			}
		}

		return best;
	}

	cover build_cover(const xmg& m, bestmap& best) {
		cout << "Building cover..." << endl;
		cover nref(m.nnodes());

		auto nodes = m.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			nref[i] = 0;
		}

		auto outputs = m.outputs();
		for (auto outid : outputs) {
			nref[outid] += 1;
		}

		int i = nodes.size();
		while (i-- > 0) {
			const auto& node = nodes[i];
			if (node.ecrep != static_cast<nodeid>(i)) {
				// Not an eq. class representative
				continue;
			}
			if (nref[i] > 0) {
				auto cutnodes = best[i]->nodes();
				for (auto nodeid : cutnodes) {
					nref[nodeid] += 1;
				}
			}
		}

		return nref;
	}

	reqmap compute_required_times(const xmg& m, const cover& nref, 
			const nintmap& arrival_times, const bestmap& best) {
		cout << "Computing required times..." << endl;
		reqmap req(m.nnodes());

		for (auto i = 0u; i < m.nnodes(); i++) {
			req[i] = numeric_limits<unsigned>::max();
		}

		auto outputs = m.outputs();
		auto maxoutcost = 0u;
		for (auto output : outputs) {
			if (arrival_times.at(output) > maxoutcost) {
				maxoutcost = arrival_times.at(output);
			}
		}
		for (auto output : outputs) {
			req[output] = maxoutcost;
		}
		int i = m.nnodes();
		while (i-- > 0) {
			if (nref.at(i) == 0) {
				continue;
			}
			auto best_nodes = best.at(i)->nodes();
			auto child_req = req[i] - 1;
			for (auto cut_node : best_nodes) {
				if (req[cut_node] > child_req) {
					req[cut_node] = child_req;
				}
			}
		}
		return req;
	}

	unsigned int recursive_deselect(const node& n, 
			const vector<node>& nodes, bestmap& best, nintmap& nref) {
		if (is_pi(n)) {
			return 0;
		}
		auto area = 1;
		auto cutnodes = best[n.ecrep]->nodes();
		for (auto nodeid : cutnodes) {
			nref[nodeid] -= 1;
			if (nref[nodeid] == 0) {
				const auto& node = nodes[nodeid];
				area += recursive_deselect(node, nodes, best, nref);
			}
		}
		return area;
	}

	unsigned int recursive_select(const node& n, 
			const vector<node>& nodes, bestmap& best, nintmap& nref) {
		if (is_pi(n)) {
			return 0;
		}
		auto area = 1;
		auto cutnodes = best[n.ecrep]->nodes();
		for (auto nodeid : cutnodes) {
			nref[nodeid] += 1;
			if (nref[nodeid] == 1) {
				const auto& node = nodes[nodeid];
				area += recursive_select(node, nodes, best, nref);
			}
		}
		return area;
	}

	void improve_cover_exact_area_timeout(const xmg& m, const cutmap& cm,
			bestmap& best, nintmap& nref, const funcmap& fm, const vector<tt>& timeoutfuncs) {
		cout << "Running exact area improvements..." << endl;
		auto nodes = m.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			} else if (node.ecrep != static_cast<nodeid>(i)) {
				// Not an eq. class representative
				continue;
			}
			if (nref.at(i) != 0) {
				recursive_deselect(node, nodes, best, nref);
			}
			auto best_area = numeric_limits<unsigned int>::max();
			cut* best_cut = NULL;
			const auto& cuts = cm.at(i);
			for (const auto& cut : cuts) {
				if (cut->size() == 1) { // Trivial cut
					continue;
				}
				const auto& cutfunction = *fm.at(cut.get());
				if (find(timeoutfuncs.begin(), timeoutfuncs.end(), cutfunction) != timeoutfuncs.end()) {
					// Don't select this function as exact synthesis will timeout on it
					continue;
				}
				best[i] = cut.get();
				auto area1 = recursive_select(node, nodes, best, nref);
				auto area2 = recursive_deselect(node, nodes, best, nref);
				assert(area1 == area2);
				if (area1 < best_area) {
					best_area = area1;
					best_cut = cut.get();
				}
			}
			best[i] = best_cut;
			if (nref.at(i) != 0) {
				recursive_select(node, nodes, best, nref);
			}
		}
	}

	void improve_cover_exact_area(const xmg& m, const cutmap& cm, 
			bestmap& best, nintmap& nref) {
		cout << "Running exact area improvements..." << endl;
		auto nodes = m.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			} else if (node.ecrep != static_cast<nodeid>(i)) {
				// Not an eq. class representative
				continue;
			}
			if (nref.at(i) != 0) {
				recursive_deselect(node, nodes, best, nref);
			}
			auto best_area = numeric_limits<unsigned int>::max();
			cut* best_cut = NULL;
			const auto& cuts = cm.at(i);
			for (const auto& cut : cuts) {
				if (cut->size() == 1) { // Trivial cut
					continue;
				}
				best[i] = cut.get();
				auto area1 = recursive_select(node, nodes, best, nref);
				auto area2 = recursive_deselect(node, nodes, best, nref);
				assert(area1 == area2);
				if (area1 < best_area) {
					best_area = area1;
					best_cut = cut.get();
				}
			}
			best[i] = best_cut;
			if (nref.at(i) != 0) {
				recursive_select(node, nodes, best, nref);
			}
		}
	}

	void 
		improve_cover_exact_area(const xmg& m, const cutmap& cm, bestmap& best, 
				cover& nref, nintmap& atimes, const reqmap& req) {
			cout << "Running exact area improvements (recover)..." << endl;
			auto nodes = m.nodes();
			for (auto i = 0u; i < nodes.size(); i++) {
				const auto& node = nodes[i];
				if (is_pi(node)) {
					continue;
				} else if (node.ecrep != static_cast<nodeid>(i)) {
					// Not an eq. class representative
					continue;
				}
				if (nref.at(i) != 0) {
					recursive_deselect(node, nodes, best, nref);
				}
				auto best_area = numeric_limits<unsigned int>::max();
				cut* best_cut = NULL;
				auto& cuts = cm.at(i);
				for (const auto& cut : cuts) {
					if (cut->size() == 1) { // Trivial cut
						continue;
					}
					best[i] = cut.get();
					auto atime = cut_arrival_time(cut.get(), atimes);
					auto area1 = recursive_select(node, nodes, best, nref);
					auto area2 = recursive_deselect(node, nodes, best, nref);
					assert(area1 == area2);
					if (area1 < best_area && atime <= req.at(i)) {
						best_area = area1;
						best_cut = cut.get();
					}
				}
				best[i] = best_cut;
				atimes[i] = cut_arrival_time(best_cut, atimes);
				if (nref.at(i) != 0) {
					recursive_select(node, nodes, best, nref);
				}
			}
			cout << endl;
		}

	unsigned int cover_size(const xmg& m, const cover& cover) {
		auto size = 0u;
		auto nodes = m.nodes();
		auto nnodes = m.nnodes();
		for (auto i = 0u; i < nnodes; i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			}
			if (cover.at(i) != 0) {
				++size;
			}
		}
		return size;
	}

	void it_exact_cover(const xmg& m, cover& cover, 
			const cutmap& cm, bestmap& best) {
		auto csize = cover_size(m, cover);
		while (true) {
			improve_cover_exact_area(m, cm, best, cover);
			auto nsize = cover_size(m, cover);
			if (nsize < csize) {
				csize = nsize;
			} else {
				break;
			}
		}
	}
	
	void it_exact_cover_timeout(const xmg& m, cover& cover, const cutmap& cm, bestmap& best, 
		const funcmap& fm, const vector<tt>& timeoutfuncs) {
		auto csize = cover_size(m, cover);
		while (true) {
			improve_cover_exact_area_timeout(m, cm, best, cover, fm, timeoutfuncs);
			auto nsize = cover_size(m, cover);
			if (nsize < csize) {
				csize = nsize;
			} else {
				break;
			}
		}
	}

	void it_exact_cover(const xmg& m, cover& cover, const cutmap& cm, 
			bestmap& best, nintmap& atimes, const reqmap& req) {
		auto csize = cover_size(m, cover);
		while (true) {
			improve_cover_exact_area(m, cm, best, cover, atimes, req);
			auto nsize = cover_size(m, cover);
			if (nsize < csize) {
				csize = nsize;
			} else {
				break;
			}
		}
	}

	unsigned int cover_depth(const xmg& m, const nintmap& atimes) {
		auto outputs = m.outputs();
		auto maxdepth = 0u;
		for (auto nodeid : outputs) {
			auto atime = atimes.at(nodeid);
			if (atime > maxdepth) {
				maxdepth = atime;
			}
		}
		return maxdepth;
	}

	funcmap compute_functions(const xmg& xmg, const cover& cover, const bestmap& best, const cutmap& cutmap) {
		funcmap fm;

		const auto& nodes = xmg.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			if (!contains(cover, i)) {
				continue;
			}
			auto& cut = best[i];
			cut->set_required();
		}

		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (i == 0) {
				const auto& cut = best[i];
				unique_ptr<tt> f(new tt(tt_const1()));
				fm[cut] = std::move(f);
				continue;
			} else if (is_pi(node)) {
				const auto& cut = best[i];
				unique_ptr<tt> f(new tt(tt_nth_var(0)));
				fm[cut] = std::move(f);
				continue;
			}
			const auto& cuts = cutmap[i];
			for (auto& cut : cuts) {
				if (cut->is_required()) {
					cut->computefunction(node, fm);
				}
			}
		}
		
		return fm;
	}

	extern void write_blif(const string&, const xmg&, const cover&, const bestmap&, const funcmap&);
	
	void lut_map(const xmg& m, const cut_params* params, const string& f) {
		const auto cut_map = enumerate_cuts(m, params);
		auto best_area = eval_matches_area(m, cut_map);
		auto area_cover = build_cover(m, best_area);
		it_exact_cover(m, area_cover, cut_map, best_area);
		auto csize = cover_size(m, area_cover);
		cout << "Cover size: " << csize << endl;
		auto functionmap = compute_functions(m, area_cover, best_area, cut_map);
		write_blif(f, m, area_cover, best_area, functionmap);
	}

}
