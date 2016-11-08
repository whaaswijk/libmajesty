#pragma once

#include <vector>
#include <truth_table_utils.hpp>

/**
 * A class and note type to represent arbitrary homogenous logic network. Functions of nodes in the
 * network are represented by truth tables.
 */

namespace majesty {
	using nodeid = uint32_t;

	struct ln_node {
		std::vector<nodeid> fanin;
		std::vector<nodeid> fanout;
		cirkit::tt function;
		bool pi;
	};

	class logic_ntk {
		private:
			std::vector<ln_node> _nodes;
			std::vector<nodeid> _outputs;
			std::vector<std::string> _innames;
			std::vector<std::string> _outnames;
			
		public:
			logic_ntk() {  }
			logic_ntk(const logic_ntk&);
			logic_ntk(logic_ntk&&);
			logic_ntk& operator=(logic_ntk&) = delete;
			logic_ntk& operator=(logic_ntk&&);

			// Counts the nr. of PI nodes
			unsigned nin() const;
			// Counts the total number of nodes
			unsigned nnodes() const { return _nodes.size(); }
			// Counts the nr. of PO nodes
			unsigned nout() const { return _outputs.size(); }
			// Counts the nr. of non-PI nodes
			unsigned ninternal() const { return _nodes.size() - nin(); }

			const std::vector<ln_node>& nodes() const { return _nodes; }
			ln_node& get_node(nodeid id) { return _nodes[id]; }
			const std::vector<nodeid>& outputs() const { return _outputs; }

			const std::vector<std::string>& innames() const { 
				return _innames;
		   	}
			
			const std::vector<std::string>& outnames() const { 
				return _outnames;
		   	}

			void add_inname(const std::string& name) {
                _innames.push_back(name);
            }

            void add_outname(const std::string& name) {
                _outnames.push_back(name);
            }

			nodeid create_input();
			nodeid create_input(const std::string&);

			nodeid create_output(nodeid);
			nodeid create_output(nodeid, bool);
			nodeid create_output(nodeid, const std::string&);
			nodeid create_output(nodeid, const std::string&, bool);

			nodeid create_node(const std::vector<nodeid>&, const cirkit::tt&);

			void create_dummy_innames();
			void create_dummy_outnames();
			void create_dummy_names();

			cirkit::tt simulate();
	};


}
