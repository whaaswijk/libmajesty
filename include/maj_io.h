#pragma once

#include <xmg.h>
#include <logic_network.h>
#include <fstream>

namespace majesty {
	xmg read_bench(const std::string);
	xmg* ptr_read_bench(const std::string);
	xmg read_bench(std::ifstream&);

	xmg read_verilog(const std::string& filename);
	xmg read_verilog(const std::string& filename, const xmg_params* frparams);
	xmg* ptr_read_verilog(const std::string&);
	xmg read_verilog(std::istream&);

	void write_verilog(const majesty::xmg&, const std::string& filename);
	void write_verilog(const majesty::xmg&, std::ostream&);

	void lut_map_area(const majesty::xmg&, const std::string&);

	void write_blif(const logic_ntk&, const std::string&);
	
	unsigned int writeVerilogMIGreduced(FILE *file, MIG *net);
}
