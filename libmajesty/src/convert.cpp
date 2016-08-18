#include <convert.h>
#include <stack>
#include <map>
#include <unordered_map>

using namespace std;
using namespace cirkit;

using bracket_map_t = unordered_map<unsigned,unsigned>;

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

	xmg xmg_from_string(const string& str, unsigned ninputs) {
		xmg res;

		nodemap nodemap;
		nodemap[0] = make_pair(res.create_input(), false);
		for (auto i = 0u; i < ninputs; i++) {
			nodemap[i + 1] = make_pair(res.create_input(string(1, 'a'+(char)i)), false);
		}

		const auto majbrackets = find_bracket_pairs(str, '<', '>');
		const auto xorbrackets = find_bracket_pairs(str, '[', ']');

		bool complout = false;
		if (str[0] == '!') {
			complout = true;
			xmg_parse_string(str, 1, majbrackets, xorbrackets, nodemap, res);
		} else {
			xmg_parse_string(str, 0, majbrackets, xorbrackets, nodemap, res);
		}
		// The output is always the outer node of the expression. This corresponds to the last node in the xmg.
		auto last_idx = res.nnodes() - 1;
		res.create_output(last_idx, complout, "F");

		return res;
	}
	
	xmg xmg_from_string(const string& str, unsigned ninputs, const tt& phase, const vector<unsigned>& perm) {
		xmg res;

		nodemap nodemap;
		nodemap[0] = make_pair(res.create_input(), false);
		for (auto i = 0u; i < ninputs; i++) {
			res.create_input(string(1, 'a'+(char)i));
		}
		for (auto i = 0u; i < ninputs; i++) {
			nodemap[perm[i] + 1] = make_pair(i+1, phase.test(i));
		}

		const auto majbrackets = find_bracket_pairs(str, '<', '>');
		const auto xorbrackets = find_bracket_pairs(str, '[', ']');

		bool complout = phase.test(ninputs);
		if (str[0] == '!') {
			complout = !complout;
			xmg_parse_string(str, 1, majbrackets, xorbrackets, nodemap, res);
		} else {
			xmg_parse_string(str, 0, majbrackets, xorbrackets, nodemap, res);
		}
		// The output is always the outer node of the expression. This corresponds to the last node in the xmg.
		auto last_idx = res.nnodes() - 1;
		res.create_output(last_idx, complout, "F");

		return res;
	}
}