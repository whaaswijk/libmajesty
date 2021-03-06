#ifndef LUT_COVER_H
#define LUT_COVER_H

#include <logic_network.h>
#include <cut.h>
#include <map>
#include <unordered_map>
#include <vector>

namespace majesty {
	using bestmap = std::vector<cut*>;
	using nintmap = std::vector<unsigned int>;
	using cover = std::vector<unsigned int>;

	class xmg;

	cover build_cover(const xmg& xmg, bestmap& best);
	cover build_cover(const logic_ntk& ntk, bestmap& best);

	bestmap eval_matches_depth(const xmg&, const cutmap&, nintmap&);
	bestmap eval_matches_area(const xmg&, const cutmap&);
	bestmap eval_matches_area_timeout(const xmg&, const cutmap&, const funcmap&, const std::vector<cirkit::tt>&);
	bestmap eval_matches_area_timeout(const logic_ntk&, const cutmap&, const funcmap&, const std::vector<cirkit::tt>&);
	bestmap eval_matches_af_recover(const xmg&, const cutmap&, const reqmap&, const nintmap&);

	static inline bool contains(const cover& cover, nodeid id) {
		return cover.at(id) != 0;
	}

	reqmap compute_required_times(
			const xmg&, const nintmap&, const nintmap&, const bestmap&);

	void improve_cover_exact_area(const xmg&, const cutmap&, bestmap&, nintmap&);
	void improve_cover_exact_area_timeout( const xmg&, const cutmap&, bestmap&, nintmap&, const funcmap&, const std::vector<cirkit::tt>&);
	void improve_cover_exact_area( const xmg&, const cutmap&, bestmap&, nintmap&, nintmap&, const reqmap&);

	void improve_cover_exact_area(const logic_ntk&, const cutmap&, bestmap&, nintmap&);
	void improve_cover_exact_area_timeout(const logic_ntk&, const cutmap&, bestmap&, nintmap&, const funcmap&, const std::vector<cirkit::tt>&);

	void it_exact_cover(const xmg&, cover&, const cutmap&, bestmap&);
	void it_exact_cover_timeout(const xmg&, cover&, const cutmap&, bestmap&, const funcmap&, const std::vector<cirkit::tt>&);
	void it_exact_cover(const xmg&, cover&, const cutmap&, bestmap&, 
			nintmap&, const reqmap&);

	void it_exact_cover(const logic_ntk& m, cover& cover, const cutmap& cm, bestmap& best);
	void it_exact_cover_timeout(const logic_ntk&, cover&, const cutmap&, bestmap&, const funcmap&, const std::vector<cirkit::tt>&);

	unsigned int cover_depth(const xmg& m, const nintmap& atimes);
	
	funcmap compute_functions(const xmg&, const bestmap&, const cutmap&);
	funcmap compute_functions(const xmg&, const cover&, const bestmap&, const cutmap&);
	
	funcmap compute_all_functions(const logic_ntk& ntk, const cutmap& cutmap);
	funcmap compute_functions(const logic_ntk& ntk, const cover& cover, const bestmap& best, const cutmap& cutmap);

	logic_ntk lut_map_area(const xmg&, const cut_params*);
	logic_ntk lut_map_area(const logic_ntk&, const cut_params*);
}

#endif
