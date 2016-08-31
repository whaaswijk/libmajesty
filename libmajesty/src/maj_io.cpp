#include <maj_io.h>
#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <map>
#include <unordered_map>
#include <convert.h>
#include <truth_table_utils.hpp>
#include <cut.h>
#include <lut_cover.h>
#include <truth_table_utils.hpp>

using namespace std;

namespace majesty {

	xmg read_bench(const string filename) {
		ifstream f(filename);
		xmg res = read_bench(f);
		return res;
	}
	
	xmg* ptr_read_bench(const std::string filename) {
		return new xmg(read_bench(filename));
	}

	xmg read_verilog(const std::string filename) {
		ifstream f(filename);
		xmg res = read_verilog(f);
		return res;
	}

	xmg* ptr_read_verilog(const std::string filename) {
		return new xmg(read_verilog(filename));
	}

	void parse_bench_file(xmg& xmg, ifstream& instream) {
		string line;
		vector<string> outputnames;
		unordered_map<string, pair<nodeid,bool>> varnames;
		nodemap nodemap;
		while (getline(instream, line)) {
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
			file << innames[i] << " , ";
		}
		file << endl << "\t\t\t";
		for(auto i = 0u; i < (xmg.nout()-1); i++) {
			file << outnames[i] << " , ";
		}
		file << outnames[xmg.nout()-1] << " ) ;" << endl;

		if (nin > 0) {
			file << "input ";
			for (auto i = 0u; i < nin - 1; i++) {
				file << innames[i] << " , ";
			}
			file << innames[nin - 1] << " ;" << endl;
		}

		file << "output ";
		for(auto i = 0u; i< (xmg.nout()-1); i++){
			file << outnames[i] << " , ";
		}	
		file << outnames[xmg.nout()-1] << " ;" << endl;

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
					file << " , ";
				} else {
					begun = true;
				}
				file << "w" << node.ecrep;
			}
			file << " ;" << endl;
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
				file << in1 << " ^ " << in2 << " ;" << endl;
			} else if (is_and(node)) {
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				auto c3 = is_c3(node);
				auto in3 = 
					node_name(node.in3, nodes, innames, c3);
				file << in2 << " & " << in3 << " ;" << endl;
			} else if (is_or(node)) {
				auto c2 = is_c2(node);
				auto in2 = 
					node_name(node.in2, nodes, innames, c2);
				auto c3 = is_c3(node);
				auto in3 = 
					node_name(node.in3, nodes, innames, c3);
				file << in2 << " | " << in3 << " ;" << endl;
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
				file << "( " << in1 << " & " << in2 << " ) | "
					<< "( " << in1 << " & " << in3 << " ) | "
					<< "( " << in2 << " & " << in3 << " ) ;" << endl;
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
			file << "assign " << outname << " = " << noderef << " ;" << endl;
		}

		file << "endmodule" << endl;
	}

	static inline void parse_verilog_inputs(xmg& xmg, const string& line, unordered_map<string,pair<nodeid,bool>>& varnames) {
		vector<string> split_strs;
		boost::split(split_strs, line, boost::is_any_of("\t ,;"), boost::token_compress_on);
		for (auto& str : split_strs) {
			boost::trim(str);
			if (str == "input" || str == "," || str == ";" || str.length() == 0) { // Ignore the input keyword and commas
				continue;
			}
			varnames[str] = make_pair(xmg.create_input(str), false);
		}
	}
	
	static inline void parse_verilog_wires(xmg& xmg, const string& line, unordered_map<string, pair<nodeid, bool>>& varnames) {
		vector<string> split_strs;
		boost::split(split_strs, line, boost::is_any_of("\t ,;"), boost::token_compress_on);
		for (auto& str : split_strs) {
			boost::trim(str);
			if (str == "wire" || str == "one" || str == "," || str == ";" || str.length() == 0) { // Ignore the input keyword and commas
				continue;
			}
			varnames[str] = xmg.create(0, false, 0, false, 0, false);
		}
	}

	static inline void parse_verilog_outputs(const string& line, vector<string>& outputnames) {
		vector<string> split_strs;
		boost::split(split_strs, line, boost::is_any_of("\t ,;"), boost::token_compress_on);
		for (auto& str : split_strs) {
			boost::trim(str);
			if (str == "output" || str == "," || str == ";" || str.size() == 0) { // Ignore the output keyword and commas
				continue;
			}
			outputnames.push_back(str);
		}
	}

	xmg read_verilog(ifstream& instream) {
		assert(instream.good());
		xmg xmg;
		unordered_map<string, pair<nodeid,bool>> varnames;
		varnames["one"] = make_pair(xmg.create_input(), false);

		string line;
		vector<string> outputnames;
		nodemap nodemap;
		bool parse_module = false, parse_inputs = false, parse_outputs = false, parse_wires = false, parse_assignments = false, endmodule = false;
		auto linenr = 0u;
		while (getline(instream, line)) {
			if (!instream) {
				break;
			}
			++linenr;
			if (line.size() == 0) {
				continue;
			}
			boost::trim(line);
			if (line.substr(0, 2) == "//") { // Comment line
				continue;
			} else if (line.substr(0, 6) == "module") { // Top declaration
				parse_module = true;
				continue;
			} else if (line.substr(0, 5) == "input") {
				parse_module = false;
				parse_inputs = true;
				parse_verilog_inputs(xmg, line, varnames);
			} else if (line.substr(0, 6) == "output") {
				parse_module = parse_inputs = false;
				parse_outputs = true;
				parse_verilog_outputs(line, outputnames);
			} else if (line.substr(0, 4) == "wire") {
				parse_module = parse_inputs = parse_outputs = false;
				parse_wires = true;
				parse_verilog_wires(xmg, line, varnames);
			} else if (line.substr(0, 6) == "assign") {
				parse_module = parse_inputs = parse_outputs = parse_wires = false;
				parse_assignments = true;

				// Remove the semicolon at the end
				auto semidx = line.find_first_of(';');
				line = line.substr(0, semidx);
				vector<string> split_strs;
				vector<string> param_names;
				boost::split(split_strs, line, boost::is_any_of("\t "), boost::token_compress_on);
				// Second entry is the name
				auto name = split_strs[1];
				if (name == "one") {
					continue;
				}

				if (split_strs.size() == 6 || split_strs.size() == 7) { // AND/OR/XOR
					auto param1 = split_strs[3];
					auto compl2 = false;
					if (param1[0] == '~') {
						compl2 = true;
						param1 = param1.substr(1, param1.size() - 1);
					}
					auto symbol = split_strs[4];
					auto param2 = split_strs[5];
					auto compl3 = false;
					if (param2[0] == '~') {
						compl3 = true;
						param2 = param2.substr(1, param2.size() - 1);
					}

					auto np1 = varnames[param1];
					auto np2 = varnames[param2];

					auto nodepfound = varnames.find(name);
					if (nodepfound == varnames.end()) { // An output is being assigned in a non-unary way. Create a new node for it.
						varnames[name] = xmg.create(0, false, 0, false, 0, false);
					}
					const auto& nodep = varnames[name];
					auto& node = xmg.get_node(nodep.first);
					auto in1 = node.in1; auto c1 = false;
					auto in2 = np1.first; 
					auto in3 = np2.first; 
					sort_inputs(in1, c1, in2, compl2, in3, compl3);

					node.in2 = in2;
					if (compl2) {
						set_c2(node);
					}
					node.in3 = in3;
					if (compl3) {
						set_c3(node);
					}
					if (symbol == "^") {
						set_xor(node);
					} else if (symbol == "&") {
						set_c1(node);
					} else {
						reset_c1(node);
					}
				} else if (split_strs.size() == 14 || split_strs.size() == 15) { // Majority
					auto param1 = split_strs[3].substr(1, split_strs[3].size() - 1);
					auto compl1 = false;
					if (param1[0] == '~') {
						compl1 = true;
						param1 = param1.substr(1, param1.size() - 1);
					}
					auto param2 = split_strs[5].substr(0, split_strs[5].size() - 1);
					auto compl2 = false;
					if (param2[0] == '~') {
						compl2 = true;
						param2 = param2.substr(1, param2.size() - 1);
					}
					auto param3 = split_strs[9].substr(0, split_strs[9].size() - 1);
					auto compl3 = false;
					if (param3[0] == '~') {
						compl3 = true;
						param3 = param3.substr(1, param3.size() - 1);
					}

					auto np1 = varnames[param1];
					auto np2 = varnames[param2];
					auto np3 = varnames[param3];
					auto in1 = np1.first; auto in2 = np2.first; auto in3 = np3.first;
					sort_inputs(in1, compl1, in2, compl2, in3, compl3);

					auto nodepfound = varnames.find(name);
					if (nodepfound == varnames.end()) { // An output is being assigned in a non-unary way. Create a new node for it.
						varnames[name] = xmg.create(0, false, 0, false, 0, false);
					}
					const auto& nodep = varnames[name];
					auto& node = xmg.get_node(nodep.first);
					node.in1 = in1;
					if (compl1) {
						set_c1(node);
					}
					node.in2 = in2;
					if (compl2) {
						set_c2(node);
					}
					node.in3 = in3;
					if (compl3) {
						set_c3(node);
					}
				} else { // Unary assignment
					assert(split_strs.size() == 4 || split_strs.size() == 5);
					auto param1 = split_strs[3];
					auto compl1 = false;
					if (param1[0] == '~') {
						compl1 = true;
						param1 = param1.substr(1, param1.size() - 1);
					}
					auto np1 = varnames[param1];
					varnames[name] = make_pair(np1.first, np1.second != compl1);
				}
			} else if (line.substr(0, 9) == "endmodule") {
				parse_module = parse_inputs = parse_outputs = parse_wires = parse_assignments = false;
				endmodule = true;
				break;
			} else if (parse_module) {
				continue;
			} else if (parse_inputs) {
				parse_verilog_inputs(xmg, line, varnames);
			} else if (parse_outputs) {
				parse_verilog_outputs(line, outputnames);
			} else if (parse_wires) { // We ignore the wire definitions
				parse_verilog_wires(xmg, line, varnames);
			} else { // Some unknown state was reached
				assert(false);
			}
		}
		for (auto& outname : outputnames) {
			auto outp = varnames[outname];
			xmg.create_output(outp.first, outp.second, outname);
		}


		return xmg;
	}

	static bool is_one(const cirkit::tt& function) {
		for (auto i = 0u; i < function.size(); i++) {
			if (!function.test(i))
				return false;
		}
		return true;
	}

	static bool is_zero(const cirkit::tt& function) {
		for (auto i = 0u; i < function.size(); i++) {
			if (function.test(i))
				return false;
		}
		return true;
	}
	
	inline string lut_name(const node& n) {
		assert(!is_pi(n) && !is_po(n));
		return "n" + to_string(n.ecrep);
	}

	// For BLIF mapping we don't require the escape slash that we use in Verilog
	string blif_name(const string& name) {
		if (name[0] == '\\') {
			return name.substr(1, name.size() - 1);
		} else {
			return string(name);
		}
	}

	string internal_name(const xmg& m, nodeid nodeid) {
		const auto& nodes = m.nodes();
		const auto& n = nodes[nodeid];
		assert(n.ecrep != 0);
		assert(n.ecrep == nodeid);
		if (is_pi(n)) {
			auto innames = m.innames();
			return blif_name(innames[n.ecrep-1]);
		} else if (is_po(n)) {
			const auto& outputs = m.outputs();
			const auto& outcompl = m.outcompl();
			const auto& outnames = m.outnames();
			bool c = false;
			auto idx = -1;
			for (auto i = 0u; i < outputs.size(); i++) {
				if (outputs[i] == n.ecrep) {
					idx = i;
					if (outcompl[i]) {
						c = true;
					}
				}
			}
			if (c) {
				return blif_name(outnames[idx] + "_bar");
			} else {
				return blif_name(outnames[idx]);
			}
		} else {
			return "n" + to_string(nodeid);
		}
	}

	// Generates a LUT for the PO and it's complement. 
	void po_lut(ofstream& f, const bestmap& best, const xmg& m,
			const cover& cover, nodeid nodeid , bool c, string lutname,
			const funcmap& funcmap) {
		const auto& mnode = m.nodes().at(nodeid);
		assert(is_po(mnode));

		auto cut = best.at(mnode.ecrep);
		const auto& inputs = cut->nodes();
		const auto& function = *funcmap.at(cut);
		if (is_zero(function)) {
			// This node never evaluates to true,
			// so we set it to constant 0
			if (c) {
				f << ".names " << lutname << endl;
				f << "1" << endl;
				if (cover.at(mnode.ecrep) > 1) {
					f << ".names " << lutname << "_bar" << endl;
				}
			} else {
				f << ".names " << lutname << endl;
				f << ".names " << lutname << "_bar" << endl;
				f << "1" << endl;
			}
			return;
		} else if (is_one(function)) {
			if (c) {
				f << ".names " << lutname << endl;
				if (cover.at(mnode.ecrep) > 1) {
					f << ".names " << lutname << "_bar" << endl;
					f << "1" << endl;
				}
			} else {
				f << ".names " << lutname << endl;
				f << "1" << endl;
				f << ".names " << lutname << "_bar" << endl;
			}
			return;
		}

		//  First print the function header
		f << ".names ";
		for (auto j = 0u; j < inputs.size(); j++) {
			f << internal_name(m, inputs[j]) << " ";
		}
		f << lutname << endl;

		// Simply generate all input vectors
		for (auto i = 0u; i < function.size(); i++) {
			auto eval = function.test(i) ^ c;
			if (eval) {
				for (auto j = 0u; j < inputs.size(); j++) {
					f << ((i >> j) & 1);
				}
				f << " " << 1 << endl;
			}
		}

		if (cover.at(mnode.ecrep) > 1 && c) {
			// Generate the complement
			f << ".names ";
			for (auto j = 0u; j < inputs.size(); j++) {
				f << internal_name(m, inputs[j]) << " ";
			}
			f << lutname << "_bar" << endl;

			for (auto i = 0u; i < function.size(); i++) {
				auto eval = function.test(i) ^ c;
				if (!eval) {
					for (auto j = 0u; j < inputs.size(); j++) {
						f << ((i >> j) & 1);
					}
					f << " " << 1 << endl;
				}
			}
		}
	}

	void in_lut(ofstream& f, const xmg& m, const node& mnode, 
			const bestmap& best, const funcmap& funcmap) {
		assert(!is_po(mnode));

		auto cut = best.at(mnode.ecrep);
		const auto& inputs = cut->nodes();
		const auto& function = *funcmap.at(cut);

		if (is_zero(function)) {
			// This node never evaluates to true,
			// so we set it to constant 0
			f << ".names " << lut_name(mnode) << endl;
			return;
		} else if (is_one(function)) {
			f << ".names " << lut_name(mnode) << endl;
			f << "1" << endl;
			return;
		}

		// Print the function header
		f << ".names ";
		for (auto j = 0u; j < inputs.size(); j++) {
			f << internal_name(m, inputs[j]) << " ";
		}
		f << lut_name(mnode) << endl;

		// Generate K-LUT

		// Simply generate all input vectors
		for (auto i = 0u; i < function.size(); i++) {
			auto eval = function.test(i);
			if (eval) {
				for (auto j = 0u; j < inputs.size(); j++) {
					f << ((i >> j) & 1);
				}
				f << " " << 1 << endl;
			}
		}
	}

	void write_blif(const string& fname, const xmg& m, const cover& cover,
			const bestmap& best, const funcmap& fm) {
		ofstream f;
		time_t now;
		time(&now);
		f.open(fname, ios::out | ios::trunc);

		f << "# Written by Majesty " << ctime(&now);

		f << ".model mapping" << endl;
		f << ".inputs ";
		const auto& innames = m.innames();
		for (auto i = 0u; i < m.nin(); i++) {
			f << blif_name(innames[i]) << " ";
		}
		f << endl;
		f << ".outputs ";
		const auto& outnames = m.outnames();
		for (auto i = 0u; i < m.nout(); i++) {
			f << blif_name(outnames[i]) << " ";
		}
		f << endl;

		const auto& outputs = m.outputs();
		const auto& outcompl = m.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			po_lut(f, best, m, cover, outputs[i], outcompl[i], blif_name(outnames[i]), fm);
		}	
	
		const auto& nodes = m.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			if (!contains(cover, i)) {
				continue;
			}
			const auto& node = nodes[i];
			if (!is_pi(node) && !is_po(node)) {
				in_lut(f, m, node, best, fm);
			}
		}

		f << ".end" << endl;
		f.close();
	}
	
	void lut_map_area(const majesty::xmg& xmg, const std::string& filename) {
		auto cparams = default_cut_params();
		const auto cut_map = enumerate_cuts(xmg, cparams.get());
		auto best_area = eval_matches_area(xmg, cut_map);
		auto area_cover = build_cover(xmg, best_area);
		it_exact_cover(xmg, area_cover, cut_map, best_area);
		auto functionmap = compute_functions(xmg, area_cover, best_area, cut_map);
		write_blif(filename, xmg, area_cover, best_area, functionmap);
	}

}