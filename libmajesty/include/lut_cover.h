#ifndef LUT_COVER_H
#define LUT_COVER_H

#include "xmg.h"
#include "cut.h"
#include <map>
#include <unordered_map>
#include <vector>

namespace majesty {
	using bestmap = std::vector<cut*>;
	using nintmap = std::vector<unsigned int>;
	using cover = std::vector<unsigned int>;

	cover build_cover(const xmg&, bestmap&);

	bestmap eval_matches_depth(const xmg&, const cutmap&, nintmap&);
	bestmap eval_matches_area(const xmg&, const cutmap&);
	bestmap eval_matches_af_recover(
			const xmg&, const cutmap&, const reqmap&, const nintmap&);

	static inline bool contains(const cover& cover, nodeid id) {
		return cover.at(id) != 0;
	}

	reqmap compute_required_times(
			const xmg&, const nintmap&, const nintmap&, const bestmap&);

	unsigned int cover_size(const xmg&, const cover&);
	void improve_cover_exact_area(
			const xmg&, const cutmap&, bestmap&, nintmap&);
	void improve_cover_exact_area(
			const xmg&, const cutmap&, bestmap&, nintmap&, 
			nintmap&, const reqmap&);

	void it_exact_cover(const xmg&, cover&, const cutmap&, bestmap&);
	void it_exact_cover(const xmg&, cover&, const cutmap&, bestmap&, 
			nintmap&, const reqmap&);

	unsigned int cover_depth(const xmg& m, const nintmap& atimes);
	
	funcmap compute_functions(const xmg&, const bestmap&, const cutmap&);
	funcmap compute_functions(const xmg&, const cover&, const bestmap&, const cutmap&);

	void lut_map(const xmg&, const cut_params*, const std::string&);

	void abc_lut_map(const xmg&, const cut_params*, const std::string&);
}

#endif
