#include <maj_io.h>
#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <map>
#include <unordered_map>
#include <convert.h>
#include <truth_table_utils.hpp>

using namespace std;

namespace majesty {

	xmg read_bench(const string filename) {
		ifstream f(filename);
		xmg res = read_bench(f);
		f.close();
		return res;
	}
	
	xmg* ptr_read_bench(const std::string filename) {
		return new xmg(read_bench(filename));
	}

	/*
	xmg read_verilog(const std::string filename) {
		ifstream f(filename);
		xmg res = read_verilog(f);
		f.close();
		return res;
	}

	xmg* ptr_read_verilog(const std::string filename) {
		return new xmg(read_verilog(filename));
	}
	*/

	void parse_bench_file(xmg& xmg, ifstream& instream) {
		string line;
		vector<string> outputnames;
		unordered_map<string, pair<nodeid,bool>> varnames;
		nodemap nodemap;
		while (!instream.eof()) {
			getline(instream, line);
			if (line.size() == 0) {
				continue;
			}
			boost::trim(line);
			if (line[0] == '#') { // Comment line
				continue;
			} else if (line.substr(0, 5) == "INPUT") {
				auto inputname = line.substr(6, line.size() - 7);
				varnames[inputname] = make_pair(xmg.create_input(inputname), false);
			} else if (line.substr(0, 6) == "OUTPUT") {
				auto outputname = line.substr(7, line.size() - 8);
				outputnames.push_back(outputname);
			} else {
				vector<string> split_strs;
				vector<string> param_names;
				boost::split(split_strs, line, boost::is_any_of("\t "), boost::token_compress_on);
				// First entry is the name
				auto name = split_strs[0];
				auto funcstr = split_strs[3];
				auto funclong = stoul(funcstr, nullptr, 16);
				// Between the brackets are the input parameters
				assert(split_strs[4] == "(");
				assert(split_strs[split_strs.size()-1] == ")");
				for (auto i = 5u; i < split_strs.size() - 1; i++) {
					auto param_name = split_strs[i];
					if (param_name == ",") {
						continue;
					}
					param_names.push_back(param_name);
				}
				auto func = tt_from_long(param_names.size(), funclong);
				auto min_xmg = exact_xmg_expression(func);
				
				nodemap[0] = make_pair(0, false);
				for (auto i = 0u; i < param_names.size(); i++) {
					auto param_name = param_names[i];
					nodemap[i + 1] = varnames[param_name];
				}
				varnames[name] = xmg_parse_string(xmg, nodemap, min_xmg);
			}
		}
		for (const auto& outname : outputnames) {
			auto outp = varnames[outname];
			xmg.create_output(outp.first, outp.second, outname);
		}
	}

	xmg read_bench(ifstream& instream) {
		assert(instream.good());
		xmg res;
		res.create_input();

		parse_bench_file(res, instream);

		return res;
	}

	void write_verilog(const string& filename, const majesty::xmg& xmg) {
		ofstream outfile(filename);
		write_verilog(outfile, xmg);
		outfile.close();
	}

	static inline string node_name(const nodeid idx, const vector<node>& nodes, const vector<string>& innames, bool c) {
			const auto& node = nodes[idx];
			if (idx == 0) {
				if (c) {
					return "1'b0";
				} else {
					return "1'b1";
				}
			} else if (is_pi(node)) {
				if (c) {
					return "~" + innames[node.ecrep - 1];
				} else {
					return innames[node.ecrep - 1];
				}
			} if (c) {
				return "~w" + to_string(node.ecrep);
			} else {
				return "w" + to_string(node.ecrep);
			}
		}

	void write_verilog(ofstream& file, const xmg& xmg) {
		time_t now;
		time(&now);
		file << "// Written by the Majesty Logic Package " << ctime(&now); 

		file << "module top (" << endl;
		file << "\t\t\t";

		const auto& innames = xmg.innames();
		const auto& outnames = xmg.outnames();
		const auto& nodes = xmg.nodes();
		const auto nin = xmg.nin();
		for(auto i = 0u; i < nin; i++){
			file << innames[i] << ", ";
		}
		file << endl << "\t\t\t";
		for(auto i = 0u; i < (xmg.nout()-1); i++) {
			file << outnames[i] << ", ";
		}
		file << outnames[xmg.nout()-1] << ");" << endl;

		if (nin > 0) {
			file << "input ";
			for (auto i = 0u; i < nin - 1; i++) {
				file << innames[i] << ", ";
			}
			file << innames[nin - 1] << ";" << endl;
		}

		file << "output ";
		for(auto i = 0u; i< (xmg.nout()-1); i++){
			file << outnames[i] << ", ";
		}	
		file << outnames[xmg.nout()-1] << ";" << endl;

		if (nin > 0) {
			file << "wire ";
			bool begun = false;
			for (auto i = 0u; i < xmg.nnodes(); i++) {
				const auto& node = nodes[i];
				if (is_pi(node)) {
					continue;
				} else if (static_cast<nodeid>(i) != node.ecrep) {
					continue;
				}
				if (begun) {
					file << ", ";
				} else {
					begun = true;
				}
				file << "w" << node.ecrep;
			}
			file << ";" << endl;
		}

		for (auto i = 0u; i < xmg.nnodes(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			} else if (static_cast<nodeid>(i) != node.ecrep) {
				continue;
			}
			file << "assign w" << node.ecrep << " = ";
			if (is_xor(node)) {
				auto c1 = is_c1(node);
				auto in1 = 
					node_name(node.in1, nodes, innames, c1);
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				file << in1 << " ^ " << in2 << ";" << endl;
			} else if (is_and(node)) {
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				auto c3 = is_c3(node);
				auto in3 = 
					node_name(node.in3, nodes, innames, c3);
				file << in2 << " & " << in3 << ";" << endl;
			} else if (is_or(node)) {
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				auto c3 = is_c3(node);
				auto in3 = 
					node_name(node.in3, nodes, innames, c3);
				file << in2 << " | " << in3 << ";" << endl;
			} else {
				auto c1 = is_c1(node);
				auto in1 = 
					node_name(node.in1, nodes, innames, c1);
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				auto c3 = is_c3(node);
				auto in3 = 
					node_name(node.in3, nodes, innames, c3);
				file << "(" << in1 << " & " << in2 << ") | "
					<< "(" << in1 << " & " << in3 << ") | "
					<< "(" << in2 << " & " << in3 << ");" << endl;
			}
		}

		const auto& outputs = xmg.outputs();
		const auto& outcompl = xmg.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			auto outid = outputs[i];
			auto c = outcompl[i];
			const auto noderef = 
				node_name(outid, nodes, innames, c);
			const auto& outname = outnames[i];
			file << "assign " << outname << " = " << noderef << ";" << endl;
		}

		file << "endmodule" << endl;
	}

}