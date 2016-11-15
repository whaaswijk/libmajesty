#include <logic_network.h>
#include <string>
#include <unordered_map>
#include <boost/functional/hash.hpp>

using namespace cirkit;

namespace majesty {

	logic_ntk::logic_ntk(const logic_ntk& ntk) {
		_nodes = ntk._nodes;
		_outputs = ntk._outputs;
		_innames = ntk._innames;
		_outnames = ntk._outnames;
	}

	logic_ntk::logic_ntk(logic_ntk&& ntk) {
		_nodes = std::move(ntk._nodes);
		_outputs = std::move(ntk._outputs);
		_innames = std::move(ntk._innames);
		_outnames = std::move(ntk._outnames);
	}

	logic_ntk& logic_ntk::operator=(logic_ntk&& ntk) {
		_nodes = std::move(ntk._nodes);
		_outputs = std::move(ntk._outputs);
		_innames = std::move(ntk._innames);
		_outnames = std::move(ntk._outnames);
		return *this;
	}

	unsigned logic_ntk::nin() const {
		auto res = 0u;
		for (const auto& node : _nodes) {
			if (node.pi) {
				++res;
			}
		}
		return res;
	}

	nodeid logic_ntk::create_input() {
		ln_node input;
		input.pi = true;
		_nodes.push_back(input);
		return _nodes.size() - 1;
	}

	nodeid logic_ntk::create_input(const std::string& name) {
		_innames.push_back(name);
        return create_input();
	}

	nodeid logic_ntk::create_output(nodeid id, bool c) {
		auto& old_node = _nodes[id];
        if (c) {
			if (old_node.fanout.size() > 0) {
				id = create_node(old_node.fanin, ~old_node.function);
			} else {
				old_node.function = ~old_node.function;
			}
        }
        _outputs.push_back(id);
        return _outputs.size() - 1;
	}
	
	nodeid logic_ntk::get_const0_node() {
		if (_has_const0_node) {
			return _const0_id;
		}
		std::vector<nodeid> zerofanin;
		_const0_id = create_node(zerofanin, tt_const0());
		_has_const0_node = true;
		return _const0_id;
	}
	
	nodeid logic_ntk::get_const1_node() {
		if (_has_const1_node) {
			return _const1_id;
		}
		std::vector<nodeid> zerofanin;
		_const1_id = create_node(zerofanin, tt_const1());
		_has_const1_node = true;
		return _const1_id;
	}

	nodeid logic_ntk::create_output(nodeid id) {
		return create_output(id, false);
	}

	nodeid logic_ntk::create_output(nodeid id, const std::string& name) {
		return create_output(id, name, false);
	}

	nodeid logic_ntk::create_output(nodeid id, const std::string& name, bool c) {
        _outnames.push_back(name);
        return create_output(id, c);
	}

	nodeid logic_ntk::create_node(const std::vector<nodeid>& fanin, const cirkit::tt& function) {
		auto hash_key = boost::hash_range(fanin.begin(), fanin.end());
		auto sh_lookup = _shmap.equal_range(hash_key);
		for (auto it = sh_lookup.first; it != sh_lookup.second; it++) {
			const auto& lookup_fanin = it->second.first;
			if (lookup_fanin == fanin) {
				// The fanin is the same, now we have to check if the function is the same too
				const auto& lookup_func = it->second.second.first;
				if (lookup_func == function) {
					// We have a node with the same fanin and the same function, strash hit!
					//std::cout << "logic_ntk strash hit!" << std::endl;
					return it->second.second.second;
				}
			}
		}

		ln_node node;
		node.pi = false;
		node.fanin = fanin;
		for (auto nodeid : fanin) {
			auto& faninnode = _nodes[nodeid];
			faninnode.fanout.push_back(_nodes.size());
		}
		node.function = function;
		_nodes.push_back(node);
		auto id = _nodes.size() - 1;
		// Insert into the strash map
		_shmap.emplace(hash_key, std::make_pair(fanin, std::make_pair(function, id)));
		return id;
	}

	void logic_ntk::create_dummy_innames() {
        _innames.clear();
		auto count = 0u;
		for (const auto& node : _nodes) {
			if (!node.pi) {
				break;
			}
			_innames.push_back("x" + std::to_string(count));
			++count;
		}
	}

	void logic_ntk::create_dummy_outnames() {
        _outnames.clear();
		for (auto i = 0u; i < _outputs.size(); i++) {
			_outnames.push_back("f[" + std::to_string(i) + "]");
		}
	}

	void logic_ntk::create_dummy_names() {
		create_dummy_innames();
		create_dummy_outnames();
	}

	static inline bool simulate_node(const ln_node& node, const std::unordered_map<nodeid, bool>& simval) {
		uint64_t func_idx = 0u;
		for (auto i = 0u; i < node.fanin.size(); i++) {
			auto sval = simval.at(node.fanin[i]);
			func_idx += (sval << i);
		}
		return node.function[func_idx];
	}

	tt logic_ntk::simulate() {
		tt func;

		const auto nsimvectors = (1u << nin());
		const auto _nnodes = nnodes();
		std::unordered_map<nodeid, bool> simval;
		for (auto j = 0u; j < nsimvectors; j++) {
			for (auto i = 0u; i < _nnodes; i++) {
				const auto& node = _nodes[i];
				if (node.pi) {
					simval[i] = (j >> i) & 1;
				} else {
					simval[i] = simulate_node(node, simval);
				}
			}
			auto outval = simval[_outputs[0]];
			func.push_back(outval);
		}
		assert(func.size() == nsimvectors);

		return func;
	}

}
