#pragma once

#include <logic_network.h>
#include <vector>
#include <truth_table_utils.hpp>
#include <map>
#include <unordered_map>
#include <boost/optional.hpp>
#include <sstream>
#include <lut_cover.h>

namespace majesty {

	class xmg;
	using nodemap = std::unordered_map<nodeid, std::pair<nodeid, bool>>;

	xmg xmg_from_string(unsigned, const std::string&);
	xmg xmg_from_string(unsigned, const std::string&, const cirkit::tt&,  const std::vector<unsigned>&);
	std::pair<nodeid,bool> xmg_parse_string(xmg&, nodemap&, const std::string&);

	std::string logic_ntk_to_string(const logic_ntk&);
	logic_ntk string_to_logic_ntk(const std::string&);
	
	cirkit::tt tt_from_long(unsigned ninputs, unsigned function);
	xmg mig_shannon_decompose(unsigned ninputs, const cirkit::tt& func);

	xmg mig_decompose(unsigned ninputs, unsigned function);

	xmg mig_decompose(unsigned ninputs, const std::string& function);
	
	std::string exact_mig_expression(const cirkit::tt&);
	boost::optional<std::string> exact_mig_expression(const cirkit::tt&, unsigned);
	std::string exact_xmg_expression(const cirkit::tt&);
	boost::optional<std::string> exact_xmg_expression(const cirkit::tt&, unsigned);
	boost::optional<std::string> exact_xmg_expression(const cirkit::tt&, unsigned, unsigned);

	boost::optional<unsigned> last_size_from_file(const std::string&);

	xmg exact_mig(const cirkit::tt&);
	xmg exact_depth_mig(const cirkit::tt&);
	xmg exact_xmg(const cirkit::tt&);
	
	boost::optional<std::string> xmg_expression_from_file(const std::string& filename);
	
	using bracket_map_t = std::unordered_map<unsigned, unsigned>;
	using input_map_t = std::unordered_map<char, std::pair<nodeid, bool>>;
	bracket_map_t find_bracket_pairs(const std::string&, char, char);

	std::string xmg_to_expr(const xmg&);
	logic_ntk xmg_to_logic_ntk(const xmg&);

	logic_ntk xmg_cover_to_logic_ntk(const xmg&, const cover&, const bestmap&, const funcmap&);
	logic_ntk ntk_cover_to_logic_ntk(const logic_ntk&, const cover&, const bestmap&, const funcmap&);

	xmg verilog_to_xmg(const std::string&);

	static inline void split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss;
		ss.str(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
	}

	static inline std::vector<std::string> split(const std::string &s, char delim) {
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	}
	
	boost::optional<std::string> min_size_expression(const cirkit::tt& func, unsigned timeout, unsigned start_size, const std::string& synth_type);
}
