#pragma once

#include <xmg.h>
#include <vector>
#include <truth_table_utils.hpp>
#include <map>
#include <unordered_map>

namespace majesty {
	xmg xmg_from_string(unsigned, const std::string&);
	xmg xmg_from_string(unsigned, const std::string&, const cirkit::tt&,  const std::vector<unsigned>&);
	std::pair<nodeid,bool> xmg_parse_string(xmg&, nodemap&, const std::string&);
	
	cirkit::tt tt_from_long(unsigned ninputs, unsigned function);
	xmg mig_shannon_decompose(unsigned ninputs, const cirkit::tt& func);

	xmg mig_decompose(unsigned ninputs, unsigned function);

	xmg mig_decompose(unsigned ninputs, const std::string& function);
	
	std::string exact_mig_expression(const cirkit::tt&);
	std::string exact_xmg_expression(const cirkit::tt&);

	xmg exact_mig(const cirkit::tt&);
	xmg exact_xmg(const cirkit::tt&);
	
	std::string xmg_expression_from_file(const std::string& filename);
	
	using bracket_map_t = std::unordered_map<unsigned, unsigned>;
	using input_map_t = std::unordered_map<char, std::pair<nodeid, bool>>;
	bracket_map_t find_bracket_pairs(const std::string&, char, char);
}