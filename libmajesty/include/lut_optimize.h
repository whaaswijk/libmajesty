#ifndef LUT_OPTIMIZE_H
#define LUT_OPTIMIZE_H

#include "xmg.h"
#include "cut.h"
#include "lut_cover.h"
#include "function_store.h"

namespace majesty {

	xmg* ptr_lut_area_strategy(const xmg&, unsigned, unsigned);
   	xmg lut_area_strategy(const xmg&, const xmg_params*, unsigned);
   	xmg mig_lut_area_strategy(const xmg&, const xmg_params*, unsigned);
   	xmg gen_lut_area_strategy(const xmg&, const xmg_params*, unsigned);


	xmg* ptr_lut_area_timeout_strategy(const xmg&, unsigned, unsigned);
	xmg lut_area_timeout_strategy(const xmg& m, const xmg_params* frparams, unsigned lut_size);

	xmg xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&);
	xmg xmg_from_luts(const xmg&, const cover&, const bestmap&, const funcmap&, std::unordered_set<unsigned long>&);
	void mine_functions(const xmg&, const xmg_params*, unsigned);
	cirkit::tt jake_canon(const cirkit::tt&, 
			unsigned* uCanonPhase, char* pCanonPerm);
}

#endif
