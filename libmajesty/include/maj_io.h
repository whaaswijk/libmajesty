#pragma once

#include <xmg.h>
#include <fstream>

namespace majesty {
	xmg read_bench(const std::string filename);
	xmg* ptr_read_bench(const std::string filename);
	xmg read_bench(std::ifstream& instream);
	void write_verilog(const std::string&, const majesty::xmg&);
	void write_verilog(std::ofstream&, const majesty::xmg&);
}