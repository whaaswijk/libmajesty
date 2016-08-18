#pragma once

#include <xmg.h>
#include <vector>
#include <truth_table_utils.hpp>

namespace majesty {
	xmg xmg_from_string(const std::string&, unsigned);
	xmg xmg_from_string(const std::string&, unsigned, const cirkit::tt&,  const std::vector<unsigned>&);
	
	xmg mig_shannon_decompose(unsigned ninputs, const cirkit::tt& func);

	xmg mig_decompose(unsigned ninputs, unsigned function);

	xmg mig_decompose(unsigned ninputs, const std::string& function);
}