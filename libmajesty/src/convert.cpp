#include <convert.h>
#include <stack>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/format.hpp>
#include <maj_io.h>
#include <sstream>

using namespace std;
using namespace cirkit;
using boost::optional;
using boost::property_tree::ptree;


namespace majesty {

	pair<nodeid, bool> xmg_parse_string(const string& expr, unsigned offset, const bracket_map_t& majbrackets,
		const bracket_map_t& xorbrackets, nodemap& nodemap, xmg& xmg) {
		assert(expr[offset] != '!');
		if (expr[offset] == '<') {
			pair<nodeid, bool> children[3];

			auto child_pos = offset + 1u;
			auto to = 0u;
			auto inv = false;
			for (auto i = 0u; i < 3; i++) {
				child_pos += to;
				inv = expr[child_pos] == '!';
				if (inv)
					++child_pos;
				if (expr[child_pos] == '<')
					to = majbrackets.at(child_pos) - child_pos + 1u;
				else if (expr[child_pos] == '[')
					to = xorbrackets.at(child_pos) - child_pos + 1u;
				else
					to = 1u;
				children[i] = xmg_parse_string(expr, child_pos, majbrackets, xorbrackets, nodemap, xmg);
				children[i].second ^= inv;
			}
			return xmg.create(
				children[0].first, children[0].second,
				children[1].first, children[1].second,
				children[2].first, children[2].second
			);
		} else if (expr[offset] == '[') {
			pair<nodeid, bool> children[2];

			auto child_pos = offset + 1u;
			auto to = 0u;
			auto inv = false;
			for (auto i = 0u; i < 2; i++) {
				child_pos += to;
				inv = expr[child_pos] == '!';
				if (inv)
					++child_pos;
				if (expr[child_pos] == '<')
					to = majbrackets.at(child_pos) - child_pos + 1u;
				else if (expr[child_pos] == '[')
					to = xorbrackets.at(child_pos) - child_pos + 1u;
				else
					to = 1u;
				children[i] = xmg_parse_string(expr, child_pos, majbrackets, xorbrackets, nodemap, xmg);
				children[i].second ^= inv;
			}
			return xmg.create(children[0].first, children[0].second, children[1].first, children[1].second);
		} else if (expr[offset] == '0') {
			return make_pair(0, true);
		} else if (expr[offset] == '1') {
			return make_pair(0, false);
		} else {
			return nodemap.at((expr[offset] - 'a') + 1);
		}
	}
	
	bracket_map_t find_bracket_pairs(const string& s, char open_bracket, char closed_bracket) {
		stack<unsigned> open;
		bracket_map_t pairs;

		for (auto i = 0u; i < s.size(); i++) {
			auto c = s[i];
			if (c == open_bracket) {
				open.push(i);
			} else if (c == closed_bracket) {
				pairs.insert({ open.top(), i });
				open.pop();
			}
		}

		assert(open.empty());

		return pairs;
	}

	xmg xmg_from_string(unsigned ninputs, const string& str) {
		xmg res;

		nodemap nodemap;
		nodemap[0] = make_pair(res.create_input(), false);
		for (auto i = 0u; i < ninputs; i++) {
			nodemap[i + 1] = make_pair(res.create_input(), false);
		}

		const auto majbrackets = find_bracket_pairs(str, '<', '>');
		const auto xorbrackets = find_bracket_pairs(str, '[', ']');

		bool complout = false;
		pair<nodeid, bool> outp;
		if (str[0] == '!') {
			complout = true;
			outp = xmg_parse_string(str, 1, majbrackets, xorbrackets, nodemap, res);
		} else {
			outp = xmg_parse_string(str, 0, majbrackets, xorbrackets, nodemap, res);
		}
		res.create_output(outp.first, outp.second != complout);

		return res;
	}

	pair<nodeid,bool> xmg_parse_string(xmg& xmg, nodemap& nodemap, const string& str) {
		const auto majbrackets = find_bracket_pairs(str, '<', '>');
		const auto xorbrackets = find_bracket_pairs(str, '[', ']');

		bool complout = false;
		pair<nodeid, bool> outp;
		if (str[0] == '!') {
			complout = true;
			outp = xmg_parse_string(str, 1, majbrackets, xorbrackets, nodemap, xmg);
		} else {
			outp = xmg_parse_string(str, 0, majbrackets, xorbrackets, nodemap, xmg);
		}
		outp.second = (outp.second != complout);
		return outp;
	}

	inline vector<unsigned> inv(const vector<unsigned>& perm) {
		vector<unsigned> invperm(perm.size());
		for ( auto i = 0u; i < perm.size(); ++i ) { invperm[perm[i]] = i; }
		return invperm;
	}
	
	xmg xmg_from_string(unsigned ninputs, const string& str, const tt& phase, const vector<unsigned>& perm) {
		xmg res;

		nodemap nodemap;
		nodemap[0] = make_pair(res.create_input(), false);
		for (auto i = 0u; i < ninputs; i++) {
			res.create_input();
		}
		auto invperm = inv(perm);
		for (auto i = 0u; i < ninputs; i++) {
			nodemap[invperm[i] + 1] = make_pair(i+1, phase.test(i));
		}

		const auto majbrackets = find_bracket_pairs(str, '<', '>');
		const auto xorbrackets = find_bracket_pairs(str, '[', ']');

		bool complout = phase.test(ninputs);
		pair<nodeid, bool> outp;
		if (str[0] == '!') {
			complout = !complout;
			outp = xmg_parse_string(str, 1, majbrackets, xorbrackets, nodemap, res);
		} else {
			outp = xmg_parse_string(str, 0, majbrackets, xorbrackets, nodemap, res);
		}
		res.create_output(outp.first, outp.second != complout);

		return res;
	}

	tt maj_tt_const0()
	{
		return boost::dynamic_bitset<>(1u, 0u);
	}

	tt maj_tt_const1()
	{
		return ~maj_tt_const0();
	}

	tt maj_tt_cof0(const tt& t, unsigned i, unsigned ninputs) {
		auto n = tt_num_vars(t);
		assert(i < n);
		auto tv = ~tt_nth_var(i);
		tv.resize(1u << ninputs);
		// tt_extend(tv, n);

		auto tc = t;
		// if (n < tt_store::i().width) { tt_extend(tc, tt_store::i().width); }

		return (tc & tv) | ((tc & tv) << (1 << i));
	}

	tt maj_tt_cof1(const tt& t, unsigned i, unsigned ninputs) {
		auto n = tt_num_vars(t);
		assert(i < n);
		auto tv = tt_nth_var(i);
		tv.resize(1u << ninputs);

		//tt_extend(tv, n);

		auto tc = t;
		//if (n < tt_store::i().width) { tt_extend(tc, tt_store::i().width); }

		return (tc & tv) | ((tc & tv) >> (1 << i));
	}

	
	pair<nodeid, bool> mig_shannon_decompose_node(xmg& mig, const nodemap& nodemap, const tt& func, const unsigned ninputs, unsigned var_idx) {
		assert(tt_num_vars(func) == ninputs);
		auto const0 = maj_tt_const0();
		tt_extend(const0, 6u);
		auto const1 = maj_tt_const1();
		tt_extend(const1, 6u);

		tt copy(func);
		tt_to_minbase(copy);
		tt_extend(copy, 6u);
		auto copystring = to_string(copy);
		if (copy == const0) {
			return make_pair(0, true);
		} else if (copy == const1) {
			return make_pair(0, false);
		}

		auto var_node = nodemap.at(var_idx+1); // +1 because idx 0 is the "one" node
		auto compl_var_node = var_node;
		compl_var_node.second = true;
		
		auto cof1 = maj_tt_cof1(func, var_idx, ninputs);
		auto cof1_node = mig_shannon_decompose_node(mig, nodemap, cof1, ninputs, var_idx + 1);
		auto cof1_and = mig.create(make_pair(0, true), cof1_node, var_node);

		auto cof0 = maj_tt_cof0(func, var_idx, ninputs);
		auto cof0_node = mig_shannon_decompose_node(mig, nodemap, cof0, ninputs, var_idx + 1);
		auto cof0_and = mig.create(make_pair(0, true), cof0_node, compl_var_node);

		return mig.create(make_pair(0, false), cof1_and, cof0_and);
	}

	xmg mig_shannon_decompose(unsigned ninputs, const tt& func) {
		assert(ninputs <= 6);
		xmg res;

		nodemap nodemap;
		nodemap[0] = make_pair(res.create_input(), false);

		for (auto i = 0u; i < ninputs; i++) {
			nodemap[i + 1] = make_pair(res.create_input(), false);
		}
		assert(tt_num_vars(func) == ninputs);
		auto output = mig_shannon_decompose_node(res, nodemap, func, ninputs, 0);
		res.create_output(output.first, output.second);

		return res;
	}

	tt tt_from_long(unsigned ninputs, unsigned function) {
		tt func(1u << ninputs, function);
		return func;
	}

	xmg mig_decompose(unsigned ninputs, unsigned function) {
		tt func(1u << ninputs, function);
		return mig_shannon_decompose(ninputs, func);
	}

	xmg mig_decompose(unsigned ninputs, const string& function) {
		tt func(function);
		return mig_shannon_decompose(ninputs, func);
	}

	optional<string> xmg_expression_from_file(const string& filename) {
		ifstream infile(filename);

		boost::property_tree::ptree pt;
		boost::property_tree::read_json(infile, pt);
		optional<string> expression;
		for (ptree::const_iterator it = pt.begin(); it != pt.end(); ++it) {
			auto obj = it->second;
			for (auto it2 = obj.begin(); it2 != obj.end(); it2++) {
				if (it2->first == "expression") {
					expression = it2->second.get_value<string>();
				}
			}
		}
		return expression;
	}

	optional<unsigned> last_size_from_file(const string& filename) {
		ifstream infile(filename);
		//stringstream infile("[{ \"command\": \"tt 1000\", \"time\" : \"2016-08-23 15:00:20\", \"tt\" : \"1000\" }, { \"command\": \"exact_mig\", \"time\" : \"2016-08-23 15:00:20\", \"min_depth\" : false, \"all_solutions\" : false, \"start\" : 1, \"runtime\" : 0.02 }, { \"command\": \"convert --mig_to_expr\", \"time\" : \"2016-08-23 15:00:20\" }, { \"command\": \"ps -e\", \"time\" : \"2016-08-23 15:00:29\", \"expression\" : \"<0ab>\" }]");

		boost::property_tree::ptree pt;
		boost::property_tree::read_json(infile, pt);
		optional<unsigned> last_size;
		for (ptree::const_iterator it = pt.begin(); it != pt.end(); ++it) {
			auto obj = it->second;
			for (auto it2 = obj.begin(); it2 != obj.end(); it2++) {
				if (it2->first == "last_size") {
					last_size = it2->second.get_value<unsigned>();
				}
			}
		}
		return last_size;
	}


	optional<string> min_size_expression(const tt& func, unsigned timeout, unsigned start_size, const string& synth_type) {
		auto cmdstr = "cirkit -l cirkit.log -c \"tt " + to_string(func) + "; exact_" + synth_type;
		if (timeout > 0) {
			cmdstr += " --timeout " + to_string(timeout);
		}
		if (start_size > 0) {
			cmdstr += " --start " + to_string(start_size);
		}
		cmdstr += "; convert --" + synth_type + "_to_expr; ps -e; quit\" > /dev/null";
		auto success = system(cmdstr.c_str());
		if (success != 0) {
			throw runtime_error("Exact synthesis through Cirkit failed");
		}
		return xmg_expression_from_file("cirkit.log");
	}


	boost::optional<std::string> exact_xmg_expression(const cirkit::tt& func, unsigned timeout, unsigned start_size) {
		return min_size_expression(func, timeout, start_size, "xmg");
	}
	
	optional<string> exact_xmg_expression(const tt& func, unsigned timeout) {
		return min_size_expression(func, timeout, 0, "xmg");
	}

	string exact_xmg_expression(const tt& func) {
		return exact_xmg_expression(func, 0).get();
	}

	optional<string> exact_mig_expression(const tt& func, unsigned timeout) {
		return min_size_expression(func, timeout, 0, "mig");
	}

	string exact_mig_expression(const tt& func) {
		return exact_mig_expression(func, 0).get();
	}

	xmg exact_mig(const tt& func) {
		auto expression = exact_mig_expression(func);
		auto ninputs = tt_num_vars(func);
		auto exact_parsed = xmg_from_string(ninputs, expression);
		return strash(exact_parsed);
	}

	xmg exact_xmg(const tt& func) {
		auto expression = exact_xmg_expression(func);
		auto ninputs = tt_num_vars(func);
		auto exact_parsed = xmg_from_string(ninputs, expression);
		return strash(exact_parsed);
	}

	static const auto heuristic_threshold = 5u;
	optional<string> heuristic_xmg_expression(const tt& func, unsigned ninputs, unsigned timeout, 
		unsigned last_size, timeout_behavior behavior) {
		// Get an upper bound for optimization by using ABC's size optimization.
		auto optcmd = (boost::format("abc -c \"read_truth %s; strash; resyn2; write_verilog tmp.v\"") % tt_to_hex(func)).str();
		auto success = system(optcmd.c_str());
		if (success != 0) {
			throw runtime_error("Heuristic optimization through ABC failed");
		}
		auto upperbound_xmg = read_verilog("tmp.v");
		auto start = upperbound_xmg.nnodes() - upperbound_xmg.nin() - 1;
		// What we do now depends on the specified behavior. If the difference between the heuristic result and the last
		// is too great it is unlikely that we will find an optimum result by going down from the starting point. We
		// We also may not want to invoke exact synthesis too often.
		if (behavior == combine) {
			if (start - last_size > heuristic_threshold) {
				return boost::none;
			}
		}
		optional<string> expr = xmg_to_expr(upperbound_xmg);
		while (true) {
			auto sat_xmg_expr = min_size_expression(func, timeout, start - 1, "xmg");
			if (sat_xmg_expr) {
				auto sat_xmg = xmg_from_string(ninputs, sat_xmg_expr.get());
				if (sat_xmg.nnodes() == upperbound_xmg.nnodes()) {
					// The upperbound found by ABC was already optimum (otherwise we would've found and XMG with size start - 1)
					expr = sat_xmg_expr;
					break;
				} else {
					--start;
				}
			} else {
				// Exact synthesis times out before finding a solution better than the heuristic one.
				break;
			}
		}

		return expr;
	}

	void node_to_expr(const vector<node>& nodes, const vector<string> innames, nodeid nodeid, stringstream& os) {
		const auto& node = nodes[nodeid];
		if (nodeid == 0) {
			os << "1";
		} else if (is_pi(node)) {
			os << innames[nodeid - 1];
		} else {
			if (is_xor(node)) {
				os << "[";
				if (is_c1(node)) {
					os << "!";
				}
				node_to_expr(nodes, innames, node.in1, os);
				if (is_c2(node)) {
					os << "!";
				}
				node_to_expr(nodes, innames, node.in2, os);
				os << "]";
			} else {
				os << "<";
				if (is_c1(node)) {
					os << "!";
				}
				node_to_expr(nodes, innames, node.in1, os);
				if (is_c2(node)) {
					os << "!";
				}
				node_to_expr(nodes, innames, node.in2, os);
				if (is_c3(node)) {
					os << "!";
				}
				node_to_expr(nodes, innames, node.in3, os);
				os << ">";
			}
		}
	}


	// NOTE: only works for single-output XMGs!
	string xmg_to_expr(const xmg& xmg) {
		stringstream os;
		const auto& nodes = xmg.nodes();
		const auto& innames = xmg.innames();
		const auto outnodeid = xmg.outputs()[0];
		const auto outnodec = xmg.outcompl()[0];
		if (outnodec) {
			os << "!";
		}
		node_to_expr(nodes, innames, outnodeid, os);
		return os.str();
	}

	static inline tt tt_from_node(const node& node) {
		tt res;

		const auto ttlen = (is_xor(node) || is_and(node) || is_or(node)) ? 4u : 8u;
		for (auto i = 0u; i < ttlen; i++) {
			bool i1pol = ((i >> 0) & 1);
			bool i2pol = ((i >> 1) & 1);
			bool i3pol = ((i >> 2) & 1);
			if (is_xor(node)) {
				const auto i1val = is_c1(node) ? !i1pol : i1pol;
				const auto i2val = is_c2(node) ? !i2pol : i2pol;
				res.push_back(i1val ^ i2val);
			} else if (is_and(node)) {
				const auto i1val = is_c2(node) ? !i1pol : i1pol;
				const auto i2val = is_c3(node) ? !i2pol : i2pol;
				res.push_back(i1val & i2val);
			} else if (is_or(node)) {
				const auto i1val = is_c2(node) ? !i1pol : i1pol;
				const auto i2val = is_c3(node) ? !i2pol : i2pol;
				res.push_back(i1val | i2val);
			} else {
				const auto i1val = is_c1(node) ? !i1pol : i1pol;
				const auto i2val = is_c2(node) ? !i2pol : i2pol;
				const auto i3val = is_c3(node) ? !i3pol : i3pol;
				res.push_back((i1val & i2val) | (i1val & i3val) | (i2val & i3val));
			}
		}

		return res;
	}
	
	logic_ntk xmg_to_logic_ntk(const xmg& xmg) {
		logic_ntk ntk;

		unordered_map<nodeid, nodeid> nodemap;

		const auto& innames = xmg.innames();
		for (auto i = 0u; i < innames.size(); i++) {
            const auto& name = innames[i];
            // Add one because we don't use the constant input
            nodemap[i+1] = ntk.create_input();
		}
		const auto& nodes = xmg.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			}
            vector<nodeid> fanin;
            if (is_xor(node)) {
                fanin.push_back(nodemap[node.in1]);
                fanin.push_back(nodemap[node.in2]);
            } else if (is_and(node) || is_or(node)) {
                fanin.push_back(nodemap[node.in2]);
                fanin.push_back(nodemap[node.in3]);
            } else {
                fanin.push_back(nodemap[node.in1]);
                fanin.push_back(nodemap[node.in2]);
                fanin.push_back(nodemap[node.in3]);
            }
			const auto node_tt = tt_from_node(node);
			nodemap[i] = ntk.create_node(fanin, node_tt);
		}

		const auto& outputs = xmg.outputs();
		const auto& outcompl = xmg.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			ntk.create_output(nodemap[outputs[i]], outcompl[i]);
		}
		
		return ntk;
	}
	
	logic_ntk xmg_cover_to_logic_ntk(const xmg& xmg, const cover& cover, const bestmap& best, const funcmap& fm) {
		logic_ntk ntk;
		unordered_map<nodeid, nodeid> nodemap;

		for (auto i = 0u; i < xmg.nin(); i++) {
            // Add one because we don't use the constant input
            nodemap[i+1] = ntk.create_input();
		}
		const auto& nodes = xmg.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (is_pi(node)) {
				continue;
			}
			if (!contains(cover, i)) {
				continue;
			}
			auto cut = best.at(node.ecrep);
			const auto& inputs = cut->nodes();
			const auto& function = *fm.at(cut);
			vector<nodeid> fanin;
			for (auto input : inputs) {
				fanin.push_back(nodemap[input]);
			}
			nodemap[i] = ntk.create_node(fanin, function);
		}

		const auto& outputs = xmg.outputs();
		const auto& outcompl = xmg.outcompl();
		for (auto i = 0u; i < outputs.size(); i++) {
			ntk.create_output(nodemap[outputs[i]], outcompl[i]);
		}

		return ntk;
	}

	logic_ntk ntk_cover_to_logic_ntk(const logic_ntk& net, const cover& cover, const bestmap& best, const funcmap& fm) {
		logic_ntk ntk;
		unordered_map<nodeid, nodeid> nodemap;

		const auto& nodes = net.nodes();
		for (auto i = 0u; i < nodes.size(); i++) {
			const auto& node = nodes[i];
			if (node.pi) {
				ntk.create_input();
				continue;
			}
			if (!contains(cover, i)) {
				continue;
			}
			auto cut = best.at(i);
			const auto& inputs = cut->nodes();
			const auto& function = *fm.at(cut);
			vector<nodeid> fanin;
			for (auto input : inputs) {
				fanin.push_back(nodemap[input]);
			}
			nodemap[i] = ntk.create_node(fanin, function);
		}

		const auto& outputs = net.outputs();
		for (auto i = 0u; i < outputs.size(); i++) {
			ntk.create_output(nodemap[outputs[i]]);
		}

		const auto& innames = net.innames();
		for (const auto& name : innames) {
			ntk.add_inname(name);
		}
		
		const auto& outnames = net.outnames();
		for (const auto& name : outnames) {
			ntk.add_inname(name);
		}

		return ntk;
	}

#ifndef _WIN32
    extern "C" void parse_verilog_string(const char *verilog_string, MIG **m);
#else
	void parse_verilog_string(const char *verilog_string, MIG **m) {

	}
#endif

	xmg verilog_to_xmg(const std::string& verilog_string) {
		MIG* mig;
		parse_verilog_string(verilog_string.c_str(), &mig);
		xmg res(mig);
		freemig(mig);
		return res;
	}
}
