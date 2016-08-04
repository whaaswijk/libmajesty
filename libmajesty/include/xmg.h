#ifndef xmg_H
#define xmg_H

#include <bitset>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "../src/minisat/Solver.h"
#include "../src/minisat/SolverTypes.h"
#include "../src/mlp.h"
#include <map>
#include <unordered_map>

#define maj3signature int32_t, bool, int32_t, bool,\
   	int32_t, bool

#define xorsignature int32_t, bool, int32_t, bool

#define maj3inputs int32_t in1, bool c1, int32_t in2, bool c2,\
   	int32_t in3, bool c3

#define xorinputs int32_t in1, bool c1, int32_t in2, bool c2

#define BITVECTOR_SIZE 8000
#define DEFAULT_BACKTRACKS 4096

struct MIG;

using nodeid = int32_t;
using bv = std::vector<unsigned int>;
using hashmap = std::unordered_map<unsigned int,MAJ3*>;
using xhashmap = std::unordered_map<unsigned int,int32_t>;
using bvmap = std::unordered_map<MAJ3*,std::shared_ptr<bv>>;
using xbvmap = std::vector<std::shared_ptr<bv>>;
using simmap = std::unordered_map<MAJ3*,MAJ3*>;
using xsimmap = std::vector<int32_t>;
using varmap = std::unordered_map<int32_t,Minisat::Var>;
using fanoutmap = std::vector<std::vector<int32_t>>;
using nodemap = std::unordered_map<int32_t,std::pair<int32_t,bool>>;

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
		int32_t ecrep, ecnext;
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

	class xmg {
		friend class strashmap;
		private:
			std::vector<node> _nodes;
			std::vector<int32_t> _outputs;
			boost::dynamic_bitset<> _outcompl;
			std::vector<std::string> _innames;
			std::vector<std::string> _outnames;
			
			int32_t create_node(maj3signature, 
					varmap&, Minisat::Solver&, fanoutmap&);
			int32_t create_node(maj3signature);
			int32_t create_node(xorsignature,
					varmap&, Minisat::Solver&, fanoutmap&);
			int32_t create_node(xorsignature);
							
			void set_choice(int32_t, int32_t, bool, fanoutmap&);
			bool majrecinfan(int32_t, int32_t);
			bool xorrecinfan(int32_t, int32_t);
			bool ecrecinfan(int32_t, int32_t);
			bool infan(int32_t, int32_t);
			void resetflag();

		public:
			xmg() { }
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
			const std::vector<int32_t>& outputs() const { return _outputs; }
			const std::vector<std::string>& innames() const { 
				return _innames;
		   	}
			const std::vector<std::string>& outnames() const { 
				return _outnames;
		   	}
			const boost::dynamic_bitset<>& outcompl() const { 
				return _outcompl;
		   	}

			int32_t create_input();
			int32_t create_input(const std::string&);
			int32_t create_input(varmap&, Minisat::Solver&);
			void create_output(int32_t, bool, const std::string&);
			std::pair<int32_t,bool> find_or_create(maj3signature, strashmap&);
			std::pair<int32_t,bool> rfind_or_create(maj3signature, strashmap&);
			std::pair<int32_t,bool> 
				find_or_create(maj3signature, strashmap&, 
						varmap&, Minisat::Solver&, fanoutmap&);
			std::pair<int32_t,bool> find_or_create(xorsignature, strashmap&);
			std::pair<int32_t,bool> 
				find_or_create(xorsignature, strashmap&, 
						varmap&, Minisat::Solver&, fanoutmap&);

			bool is_mig() const;
			MIG* extractmig() const;
	};
}

#endif
