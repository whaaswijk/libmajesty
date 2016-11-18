#pragma once

#include "cut.h"
#include "lut_cover.h"
#include "function_store.h"
#include <logic_network.h>

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

	template<typename S>
	logic_ntk lut_area_strategy(const logic_ntk&, unsigned lut_size, unsigned conflict_limit = 0, timeout_behavior = optimize_heuristically);

	xmg xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&);
	boost::optional<xmg> xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&, std::vector<cirkit::tt>&, unsigned);
	boost::optional<xmg> xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&, std::vector<cirkit::tt>&, unsigned, timeout_behavior);

	boost::optional<logic_ntk> logic_ntk_from_luts(const logic_ntk& lut_ntk, std::vector<cirkit::tt>&, unsigned, timeout_behavior);

	void mine_functions(const xmg&, const xmg_params*, unsigned);
	cirkit::tt jake_canon(const cirkit::tt&, unsigned* uCanonPhase, char* pCanonPerm, unsigned num_vars);

	template<typename S>
	logic_ntk size_rewrite_strategy(const logic_ntk&, unsigned cut_size, unsigned conflict_limit, timeout_behavior = optimize_heuristically);

	logic_ntk abc_heuristic_logic_ntk(const cirkit::tt&);
}
