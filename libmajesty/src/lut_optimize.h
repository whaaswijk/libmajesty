#ifndef LUT_OPTIMIZE_H
#define LUT_OPTIMIZE_H

#include "xmg.h"
#include "cut.h"
#include "lut_cover.h"

namespace majesty {

	enum NODE_TYPE { ALL, MAJ };
	

   	xmg lut_area_strategy(const xmg&, const xmg_params*, unsigned);
   	xmg mig_lut_area_strategy(const xmg&, const xmg_params*, unsigned);
   	xmg gen_lut_area_strategy(const xmg&, const xmg_params*, unsigned, NODE_TYPE);
	xmg xmg_from_luts(const xmg&, const cover&, 
			const bestmap&, const funcmap&, NODE_TYPE);
	void mine_functions(const xmg&, const xmg_params*, unsigned);
	std::string min_size_depth_xmg(cirkit::tt npn, NODE_TYPE);
	cirkit::tt jake_canon(const cirkit::tt&, 
			unsigned* uCanonPhase, char* pCanonPerm);
}

#endif
