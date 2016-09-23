#include <logic_network.h>
#include <string>

using namespace cirkit;

namespace majesty {

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
			if (!node.pi) {
				break;
			}
			++res;
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
		ln_node node;
		node.pi = false;
		node.fanin = fanin;
		for (auto nodeid : fanin) {
			auto& faninnode = _nodes[nodeid];
			faninnode.fanout.push_back(_nodes.size());
		}
		node.function = function;
		_nodes.push_back(node);
		return _nodes.size() - 1;
	}

	void logic_ntk::create_dummy_innames() {
        _innames.clear();
		auto count = 0u;
		for (const auto& node : _nodes) {
			if (!node.pi) {
				break;
			}
			_innames.push_back("x_" + std::to_string(count));
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

}