#ifndef xmg_H
#define xmg_H

#include <bitset>
#include <vector>
#include <limits>
#include <boost/dynamic_bitset.hpp>
#include "../src/minisat/Solver.h"
#include "../src/minisat/SolverTypes.h"
#include "../src/mlp.h"
#include <map>
#include <unordered_map>
#include <functional>

using nodeid = uint32_t;
using input = std::pair<nodeid,bool>;
#define EC_NULL std::numeric_limits<uint32_t>::max()
using bv = std::vector<unsigned int>;
using hashmap = std::unordered_map<unsigned int,MAJ3*>;
using xhashmap = std::unordered_map<unsigned int,nodeid>;
using bvmap = std::unordered_map<MAJ3*,std::shared_ptr<bv>>;
using xbvmap = std::vector<std::shared_ptr<bv>>;
using simmap = std::unordered_map<MAJ3*,MAJ3*>;
using xsimmap = std::vector<nodeid>;
using varmap = std::unordered_map<nodeid,Minisat::Var>;
using fanoutmap = std::vector<std::vector<nodeid>>;
using nodemap = std::unordered_map<nodeid,std::pair<nodeid,bool>>;


#define maj3signature nodeid, bool, nodeid, bool,\
   	nodeid, bool

#define xorsignature nodeid, bool, nodeid, bool

#define maj3inputs nodeid in1, bool c1, nodeid in2, bool c2,\
   	nodeid in3, bool c3

#define xorinputs nodeid in1, bool c1, nodeid in2, bool c2

#define BITVECTOR_SIZE 8000
#define DEFAULT_BACKTRACKS 4096

struct MIG;

/* function pair
Synopsis [NxN->N pairing function]

Description [Cantor pairing function to uniquely map 2 integer into another
integer]

Side effects [-]
*/ 
static inline unsigned cantor_pair(unsigned i, unsigned j) { 
	unsigned int p;

	p=((i+j)*(i+j+1)/2)+i;

	return p;
}


#define SWAP(x, xc, y, yc) if (y < x) {\
	tn = x; tc = xc; x = y; xc = yc; y = tn; yc = tc; }

static inline void sort_inputs(
		nodeid& in1, bool& c1, 
		nodeid& in2, bool& c2,
		nodeid& in3, bool& c3) {
		nodeid tn; bool tc;
		SWAP(in2, c2, in3, c3);
		SWAP(in1, c1, in3, c3);
		SWAP(in1, c1, in2, c2);
}

static inline void 
sort_inputs(nodeid& in1, bool& c1, nodeid& in2, bool& c2) {
		nodeid tn; bool tc;
		SWAP(in1, c1, in2, c2);
}

namespace majesty {
	class strashmap;

	struct xmg_params {
		unsigned nr_backtracks;
	};	
	std::unique_ptr<xmg_params> default_xmg_params(); 

	struct node {
		nodeid in1, in2, in3;
		/**
		 * Bit 0: is primary input
		 * Bit 1: is primary output
		 * Bit 2: in1 is complemented
		 * Bit 3: in2 is complemented
		 * Bit 4: in3 is complemented
		 * Bit 5: is complement of eq. class representative
		 * Bit 6: flag for misc purposes
		 * Bit 7: if set this is an XOR node
		 */
		uint8_t flag;
		nodeid ecrep, ecnext;
	};

	struct xmg_stats {
		unsigned strash_hits;
		unsigned nr_potentials;
		unsigned nr_matches;
		unsigned nr_misses;
		unsigned nr_undefined;
	};

	static inline bool is_set(const uint8_t& flag, uint8_t idx) {
		return (flag >> idx) & 1;
	}

	static inline void set_bit(uint8_t& flag, uint8_t idx) {
		flag |= (1 << idx);
	}

	static inline void clear_bit(uint8_t& flag, uint8_t idx) {
		flag &= ~(1 << idx);
	}

	static inline bool is_pi(const node& n) {
		return is_set(n.flag, 0);
	}

	static inline bool is_po(const node& n) {
		return is_set(n.flag, 1);
	}

	static inline bool is_c1(const node& n) {
		return is_set(n.flag, 2);
	}

	static inline bool is_c2(const node& n) {
		return is_set(n.flag, 3);
	}

	static inline bool is_c3(const node& n) {
		return is_set(n.flag, 4);
	}

	static inline bool is_c(const node& n) {
		return is_set(n.flag, 5);
	}

	static inline bool is_pi_c(const node& n) {
		return is_pi(n) && is_c1(n);
	}

	static inline bool is_flag_set(const node& n) {
		return is_set(n.flag, 6);
	}

	static inline bool is_xor(const node& n) {
		return is_set(n.flag, 7);
	}

	static inline bool is_and(const node& n) {
		return n.in1 == 0 && is_c1(n);
	}

	static inline bool is_or(const node& n) {
		return n.in1 == 0 && !is_c1(n);
	}

	static inline bool is_maj(const node& n) {
		return (!is_set(n.flag, 7));
	}

	static inline void set_pi(node& n) {
		set_bit(n.flag, 0);
	}

	static inline void set_po(node& n) {
		set_bit(n.flag, 1);
	}

	static inline void set_c1(node& n) {
		set_bit(n.flag, 2);
	}

	static inline void set_c2(node& n) {
		set_bit(n.flag, 3);
	}

	static inline void set_c3(node& n) {
		set_bit(n.flag, 4);
	}

	static inline void set_c(node& n) {
		set_bit(n.flag, 5);
	}
	
	static inline void set_flag(node& n) {
		set_bit(n.flag, 6);
	}

	static inline void set_xor(node& n) {
		set_bit(n.flag, 7);
	}
	
	static inline void reset_c1(node& n) {
		clear_bit(n.flag, 2);
	}
	
	static inline void reset_c2(node& n) {
		clear_bit(n.flag, 3);
	}

	static inline void reset_c3(node& n) {
		clear_bit(n.flag, 4);
	}

	static inline void reset_flag(node& n) {
		clear_bit(n.flag, 6);
	}

	static inline std::vector<nodeid> fanin(const node& node) {
		if (is_xor(node) || is_and(node) || is_or(node)) {
			return{ node.in2, node.in3 };
		} else {
			return{ node.in1, node.in2, node.in3 };
		}
	}

	static inline std::vector<std::pair<nodeid, bool>> get_children(const node& node) {
		return { std::make_pair(node.in1, is_c1(node)), std::make_pair(node.in2, is_c2(node)), std::make_pair(node.in3, is_c3(node)) };
	}

	// Drops the first occurrences of a child from a vector of nodes.
	static inline std::vector<std::pair<nodeid, bool>> drop_child(const std::vector<std::pair<nodeid,bool>> children, std::pair<nodeid,bool> child) {
		std::vector<std::pair<nodeid, bool>> res;
		bool dropped = false;
		for (const auto& c : children) {
			if (!dropped && c.first == child.first && c.second == child.second) {
				dropped = true;
				continue;
			}
			res.push_back(c);
		}
		return res;
	}

	static inline std::vector<std::pair<nodeid, bool>> filter_nodes(const std::vector<std::pair<nodeid, bool>>& nodes, 
		std::function<bool(std::pair<nodeid,bool>)> filter) {
		std::vector<std::pair<nodeid, bool>> res;
		for (const auto& np : nodes) {
			if (filter(np)) {
				res.push_back(np);
			}
		}
		return res;
	}

	class xmg {
		friend class strashmap;
		private:
			std::vector<node> _nodes;
			std::vector<nodeid> _outputs;
			std::vector<bool> _outcompl;
			std::vector<std::string> _innames;
			std::vector<std::string> _outnames;
			
			nodeid create_node(maj3signature, 
					varmap&, Minisat::Solver&, fanoutmap&);
			nodeid create_node(maj3signature);
			nodeid create_node(xorsignature,
					varmap&, Minisat::Solver&, fanoutmap&);
			nodeid create_node(xorsignature);
							
			void set_choice(nodeid, nodeid, bool, fanoutmap&);
			bool majrecinfan(nodeid, nodeid);
			bool xorrecinfan(nodeid, nodeid);
			bool ecrecinfan(nodeid, nodeid);

		public:
			xmg() {  }
			xmg(const xmg&) = delete;
			xmg(xmg&&);
			xmg& operator=(xmg&&);
			xmg(MIG* mig, const xmg_params*);
			xmg(const xmg&, const xmg_params*);
			// Counts the nr. of PI nodes
			unsigned nin() const;
			// Counts the nr. of non-PI nodes
			unsigned nnodes() const { return _nodes.size(); }
			// Counts the nr. of PO nodes
			unsigned nout() const { return _outputs.size(); }
			const std::vector<node>& nodes() const { return _nodes; }
			node& get_node(nodeid id) { return _nodes[id]; }
			const std::vector<nodeid>& outputs() const { return _outputs; }
			const std::vector<std::string>& innames() const { 
				return _innames;
		   	}
			const std::vector<std::string>& outnames() const { 
				return _outnames;
		   	}
			const std::vector<bool>& outcompl() const { 
				return _outcompl;
		   	}

			nodeid create_input();
			nodeid create_input(bool c);
			nodeid create_input(const std::string&);
			nodeid create_input(varmap&, Minisat::Solver&);
			nodeid create_input(const std::string&, varmap&, Minisat::Solver&);
			void create_output(nodeid, bool, const std::string&);
			// Creates a raw node, without strashing or propagation
			std::pair<nodeid,bool> create(maj3signature);
			std::pair<nodeid,bool> create(input, input, input);
			// Creates a node with propagation but without strashing
			std::pair<nodeid,bool> prop_create(maj3signature);
			std::pair<nodeid,bool> prop_create(input, input, input);
			// Creates a node using both strashing and propagation
			std::pair<nodeid,bool> find_or_create(maj3signature, strashmap&);
			std::pair<nodeid,bool> find_or_create(std::pair<nodeid,bool>, std::pair<nodeid,bool>, std::pair<nodeid,bool>, strashmap&);
			std::pair<nodeid,bool> 
				find_or_create(maj3signature, strashmap&, 
						varmap&, Minisat::Solver&, fanoutmap&);
			// Creates a raw node, without strashing or propagation
			std::pair<nodeid,bool> create(xorsignature);
			std::pair<nodeid, bool> create(input, input);
			std::pair<nodeid,bool> find_or_create(xorsignature, strashmap&);
			std::pair<nodeid,bool> 
				find_or_create(xorsignature, strashmap&, 
						varmap&, Minisat::Solver&, fanoutmap&);
			bool infan(nodeid, nodeid);
			void resetflag();

			bool is_mig() const;
			MIG* extractmig() const;

			// Check for equality using combinational equivalence checking
			bool equals(const xmg&) const;
	};

	xmg strash(const xmg&);

	

	// Simulates every possible input vector on an xmg and returns the function it computes
	boost::dynamic_bitset<> simulate_xmg(const xmg&);
}

#endif
