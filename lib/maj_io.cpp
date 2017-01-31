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
#include <cstdio>

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

	extern "C" void parse_verilog(FILE *file, MIG **m);

	xmg read_verilog(const std::string& filename) {
		MIG* mig;
		auto fp = fopen(filename.c_str(), "r");
		if (fp == NULL) {
			throw runtime_error("Unable to open input file");
		}
		parse_verilog(fp, &mig);
		fclose(fp);
		xmg res(mig);
		freemig(mig);
		return res;
	}

	xmg read_verilog(const std::string& filename, const xmg_params* frparams) {
		MIG* mig;
		auto fp = fopen(filename.c_str(), "r");
		if (fp == NULL) {
			throw runtime_error("Unable to open input file");
		}
		parse_verilog(fp, &mig);
		fclose(fp);
		xmg res(mig, frparams);
		freemig(mig);
		return res;
	}

	xmg* ptr_read_verilog(const std::string& filename) {
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

	void write_verilog(const majesty::xmg& xmg, const string& filename) {
		ofstream outfile(filename);
		write_verilog(xmg, outfile);
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

	void write_verilog(const xmg& xmg, ostream& file) {
		file << "// Written by the Majesty Logic Package\n";

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

	// For BLIF mapping we don't require the escape slash that we use in Verilog
	string blif_name(const string& name) {
		if (name[0] == '\\') {
			return name.substr(1, name.size() - 1);
		} else {
			return string(name);
		}
	}

	static inline string nodename(nodeid id, const vector<nodeid>& outputs, const vector<string>& outnames) {
		for (auto i = 0u; i < outputs.size(); i++) {
			if (outputs[i] == id) {
				return outnames[i];
			}
		}
		return "n_" + to_string(id);
	}

	void write_blif(const logic_ntk& ntk, const std::string& fname) {
		ofstream f;
		f.open(fname, ios::out | ios::trunc);

		f << "# Written by the Majesty Logic Package\n";

		f << ".model mapping" << endl;
		f << ".inputs ";

		const auto& innames = ntk.innames();
		assert(innames.size() > 0);
		for (auto i = 0u; i < ntk.nin(); i++) {
			f << blif_name(innames[i]) << " ";
		}
		f << endl;
		f << ".outputs ";
		const auto& outnames = ntk.outnames();
		assert(outnames.size() > 0);
		for (auto i = 0u; i < ntk.nout(); i++) {
			f << blif_name(outnames[i]) << " ";
		}
		f << endl;

		const auto& nodes = ntk.nodes();
		const auto& outputs = ntk.outputs();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (node.pi) {
				continue;
			}
			f << ".names ";
			for (auto nodeid : node.fanin) {
				string name;
				const auto& fanin_node = nodes[nodeid];
				if (fanin_node.pi) {
					name = innames[nodeid];
				} else {
					name = nodename(nodeid, outputs, outnames);
				}
				f << name << " ";
			}
			f << nodename(i, outputs, outnames) << endl;
			// Check if it's a constant node
			if (ntk.has_const0_node() && i == ntk.const0_id()) {
				// Const 0: no line necessary 
			} else if (ntk.has_const1_node() && i == ntk.const1_id()) {
				// Const 1
				f << "1" << endl;
			} else {
				// Simply print the truth table implemented by this node
				const auto& tt = node.function;
				for (auto j = 0u; j < tt.size(); j++) {
					if (tt.test(j)) {
						for (auto k = 0u; k < node.fanin.size(); k++) {
							f << ((j >> k) & 1u);
						}
						f << " " << tt.test(j) << endl;
					}
				}
			}
		}

		// Check if any output was mapped to a (constant) input
		for (auto i = 0u; i < outputs.size(); i++) {
			const auto outid = outputs[i];
			const auto& node = nodes[outid];
			if (node.pi) {
				f << ".names " << innames[outid] << " " << outnames[i] << endl;
				f << "1 1" << endl;
			}
		}

		f << ".end";
	}

	unsigned int writeVerilogMIGreduced(FILE *file, MIG *net) {
		unsigned int ret, i;
		ret = 1;

		time_t now;

		time(&now);

		fprintf(file, "//Written by the Majority Logic Package %s", ctime(&now));
		fprintf(file, "module top (\n");
		fprintf(file, "            ");
		for (i = 0;i < net->Nin;i++) {
			fprintf(file, "%s , ", net->innames[i]);
		}
		fprintf(file, "\n            ");
		for (i = 0;i < net->Nout - 1;i++) {
			fprintf(file, "%s , ", net->outnames[i]);
		}
		fprintf(file, "%s ) ;\n", net->outnames[net->Nout - 1]);

		fprintf(file, "input ");
		for (i = 0;i < net->Nin - 1;i++) {//---
			fprintf(file, "%s , ", net->innames[i]);
		}
		fprintf(file, "%s ;\n", net->innames[net->Nin - 1]);

		fprintf(file, "output ");
		for (i = 0;i < net->Nout - 1;i++) {
			fprintf(file, "%s , ", net->outnames[i]);
		}
		fprintf(file, "%s ;\n", net->outnames[net->Nout - 1]);

		fprintf(file, "wire one");
		for (i = 0u;i < net->Nnodes;i++) {
			fprintf(file, " , w%u", net->nodes[i]->label);
		}
		fprintf(file, " ;\n");

		for (i = 0;i < net->Nnodes;i++) {
			if ((net->nodes[i]->in1 == net->one) || (net->nodes[i]->in2 == net->one) || (net->nodes[i]->in3 == net->one)) {
				if (net->nodes[i]->in1 == net->one) {
					if (net->nodes[i]->compl1 == 1) {//AND
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl2 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in2 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in2->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in2->label);
						}
						fprintf(file, " & ");
						if (net->nodes[i]->compl3 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in3 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in3->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in3->label);
						}
						fprintf(file, " ;\n");
					} else {//OR
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl2 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in2 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in2->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in2->label);
						}
						fprintf(file, " | ");
						if (net->nodes[i]->compl3 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in3 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in3->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in3->label);
						}
						fprintf(file, " ;\n");
					}
				} else if (net->nodes[i]->in2 == net->one) {
					if (net->nodes[i]->compl2 == 1) {//AND
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl1 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in1 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in1->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in1->label);
						}
						fprintf(file, " & ");
						if (net->nodes[i]->compl3 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in3 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in3->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in3->label);
						}
						fprintf(file, " ;\n");
					} else {//OR
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl1 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in1 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in1->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in1->label);
						}
						fprintf(file, " | ");
						if (net->nodes[i]->compl3 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in3 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in3->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in3->label);
						}
						fprintf(file, " ;\n");
					}
				} else if (net->nodes[i]->in3 == net->one) {
					if (net->nodes[i]->compl3 == 1) {//AND
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl2 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in2 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in2->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in2->label);
						}
						fprintf(file, " & ");
						if (net->nodes[i]->compl1 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in1 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in1->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in1->label);
						}
						fprintf(file, " ;\n");
					} else {//OR
						fprintf(file, "assign w%u = ", net->nodes[i]->label);
						if (net->nodes[i]->compl2 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in2 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in2->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in2->label);
						}
						fprintf(file, " | ");
						if (net->nodes[i]->compl1 == 1) {
							fprintf(file, "~");
						}
						if (net->nodes[i]->in1 == net->one) {
							fprintf(file, "one");
						} else if (net->nodes[i]->in1->PI == 1) {
							fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
						} else {
							fprintf(file, "w%u", net->nodes[i]->in1->label);
						}
						fprintf(file, " ;\n");
					}
				} else {
					printf("error in MIG writing");
				}
			} else {
				fprintf(file, "assign w%u = ", net->nodes[i]->label);
				fprintf(file, "(");
				if (net->nodes[i]->compl1 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in1 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in1->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in1->label);
				}
				fprintf(file, " & ");
				if (net->nodes[i]->compl2 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in2 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in2->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in2->label);
				}
				fprintf(file, ") | ");
				fprintf(file, "(");
				if (net->nodes[i]->compl1 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in1 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in1->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in1->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in1->label);
				}
				fprintf(file, " & ");
				if (net->nodes[i]->compl3 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in3 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in3->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in3->label);
				}
				fprintf(file, ") | ");
				fprintf(file, "(");
				if (net->nodes[i]->compl2 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in2 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in2->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in2->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in2->label);
				}
				fprintf(file, " & ");
				if (net->nodes[i]->compl3 == 1) {
					fprintf(file, "~");
				}
				if (net->nodes[i]->in3 == net->one) {
					fprintf(file, "one");
				} else if (net->nodes[i]->in3->PI == 1) {
					fprintf(file, "%s", net->innames[net->nodes[i]->in3->label]);
				} else {
					fprintf(file, "w%u", net->nodes[i]->in3->label);
				}
				fprintf(file, ") ;\n");
			}
		}

		fprintf(file, "assign one = 1;\n");

		for (i = 0;i < net->Nout;i++) {
			if (net->outcompl[i] == 1) {
				if (net->out[i] == net->one) {
					fprintf(file, "assign %s = ~one ;// level 0\n", net->outnames[i]);
				} else if (net->out[i]->PI == 1) {
					fprintf(file, "assign %s = ~%s ;// level 0\n", net->outnames[i], net->innames[net->out[i]->label]);
				} else {
					fprintf(file, "assign %s = ~w%u ;// level %u\n", net->outnames[i], net->out[i]->label, net->out[i]->level);
				}
			} else {
				if (net->out[i] == net->one) {
					fprintf(file, "assign %s = one ;// level 0\n", net->outnames[i]);
				} else if (net->out[i]->PI == 1) {
					fprintf(file, "assign %s = %s ;// level 0\n", net->outnames[i], net->innames[net->out[i]->label]);
				} else {
					fprintf(file, "assign %s = w%u ;// level %u\n", net->outnames[i], net->out[i]->label, net->out[i]->level);
				}
			}
		}


		fprintf(file, "endmodule\n");


		return ret;
	}

}
