#pragma once

#include <xmg.h>
#include <fstream>

namespace majesty {
	xmg read_bench(const std::string);
	xmg* ptr_read_bench(const std::string);
	xmg read_bench(std::ifstream&);

	xmg read_verilog(const std::string);
	xmg* ptr_read_verilog(const std::string);
	xmg read_verilog(std::ifstream&);

	void write_verilog(const std::string&, const majesty::xmg&);
	void write_verilog(std::ofstream&, const majesty::xmg&);
}