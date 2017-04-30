#pragma once

#include <random>
#include <xmg.h>

namespace majesty {

	enum MoveType {
		MAJ3_PROP = 0,
		IDENTITY,
		INVERTER_PROP,
		DIST_RIGHT_LEFT,
		SWAP_TERNARY,
		DIST_LEFT_RIGHT,
		SUBSTITUTION,
		RELEVANCE,
		MAJ3_XXY,
		MAJ3_XYY
	};

	struct move {
		MoveType type;
		nodeid nodeid1;
		nodeid nodeid2;
		nodeid nodeid3;

		bool operator==(const move &o) const
		{
			return type == o.type && nodeid1 == o.nodeid1 && nodeid2 == o.nodeid2 && nodeid3 == o.nodeid3;
		}
		
		bool operator<(const move &o) const
		{
			return (type < o.type) ||
				(type == o.type && nodeid1 < o.nodeid1) ||
				(type == o.type && nodeid1 == o.nodeid1 && nodeid2 < o.nodeid2) ||
				(type == o.type && nodeid1 == o.nodeid1 && nodeid2 == o.nodeid2 && nodeid3 < o.nodeid3);
		}

	};

	struct partial_move {
		move move;
		int filled;
	};

	const unsigned NR_UNARY_MOVES = 3;
	const unsigned NR_BINARY_MOVES = 1;
	const unsigned NR_TERNARY_MOVES = 4;
	const unsigned NR_EDGE_TYPES = 2;
	const unsigned DEFAULT_SEED = 100;

	inline unsigned get_nr_unary_moves() {
		return NR_UNARY_MOVES;
	}

	inline unsigned get_nr_binary_moves() {
		return NR_BINARY_MOVES;
	}

	inline unsigned get_nr_ternary_moves() {
		return NR_TERNARY_MOVES;
	}

	inline unsigned get_nr_edge_types() {
		return NR_EDGE_TYPES;
	}

	class mig_manager {
	private:
		unsigned seed = 0;
		std::mt19937_64 rng;

	public:
		mig_manager() : mig_manager(DEFAULT_SEED) {
		}

		mig_manager(unsigned seed) : seed(seed) {
			rng.seed(seed);
		}

		unsigned get_seed() {
			return seed;
		}

		void set_seed(unsigned seed) {
			rng.seed(seed);
		}

		xmg* create_random_graph(unsigned ninputs, unsigned nnodes);
		xmg* create_random_graph(unsigned ninputs, unsigned nnodes, unsigned noutputs);
		xmg* random_mig_decomposition(unsigned ninputs);
	};

	xmg* apply_move(const xmg&, const move&);
	float compute_reward(const xmg&, const xmg&);
	std::vector<move> compute_moves(const xmg&);
	std::vector<move> compute_moves_fast(const xmg&, const unsigned max_nr_moves = 0);
	bool partial_move_applies(const std::vector<node>& nodes, const nodeid nid, const partial_move& pm);
	xmg* mig_string_decompose(const std::string& truth_table);
	xmg* mig_expression_decompose(unsigned ninputs, const std::string& expr);
	xmg* mig_int_decompose(unsigned ninputs, unsigned truth_table);
	xmg* get_optimum_mig(const xmg&);
	xmg* get_depth_optimum_mig(const xmg&);
	xmg* get_npn_representative(const xmg&);
	unsigned long get_truth_table(const xmg&);
	xmg* get_optimum_xmg(const xmg&);
	xmg* strash_xmg(const xmg&, bool no_compl=true);
	xmg* remove_duplicates(const xmg&);
    xmg* verilog_to_xmg_ptr(const std::string&);

	xmg* resyn2(const xmg&, const std::string& tmpfilename = "tmp.v");

	// Counts the number of shared inputs (with same polarity)
	inline char count_shared_inputs(const node& n1, const node& n2) {
		char share1_counter = 0;
		char share2_counter = 0;
		char share3_counter = 0;

		share1_counter += (n1.in1 == n2.in1 && is_c1(n1) == is_c1(n2));
		share1_counter += (n1.in1 == n2.in2 && is_c1(n1) == is_c2(n2));
		share1_counter += (n1.in1 == n2.in3 && is_c1(n1) == is_c3(n2));
		if (share1_counter > 1) share1_counter = 1;

		share2_counter += (n1.in2 == n2.in1 && is_c2(n1) == is_c1(n2));
		share2_counter += (n1.in2 == n2.in2 && is_c2(n1) == is_c2(n2));
		share2_counter += (n1.in2 == n2.in3 && is_c2(n1) == is_c3(n2));
		if (share2_counter > 1) share2_counter = 1;

		share3_counter += (n1.in3 == n2.in1 && is_c3(n1) == is_c1(n2));
		share3_counter += (n1.in3 == n2.in2 && is_c3(n1) == is_c2(n2));
		share3_counter += (n1.in3 == n2.in3 && is_c3(n1) == is_c3(n2));
		if (share3_counter > 1) share3_counter = 1;
		
		return share1_counter + share2_counter + share3_counter;
	}

	// Checks if two nodes share two inputs with the same polarity
	inline bool share_input_polarity(const node& n1, const node& n2) {
		auto nshared = count_shared_inputs(n1, n2);
		return nshared > 0;
	}

	// Checks if two nodes share two inputs with the same polarity
	inline bool share_two_input_polarity(const node& n1, const node& n2) {
		auto nshared = count_shared_inputs(n1, n2);
		return nshared > 1;
	}

	inline bool maj3_applies(const std::vector<node>& nodes, const node& n) {
		return !is_pi(n) && ((n.in1 == n.in2) || (n.in1 == n.in3) || (n.in2 == n.in3));
	}

	inline bool pm_start_inv_prop(const std::vector<node>& nodes, const node& n) {
		return !is_pi(n);
	}

	// Applies if n has 2 child nodes that have 2 nodes in common
	inline bool pm_start_dist_right_left(const std::vector<node>& nodes, const node& n) {
		if (is_pi(n))
			return false;
		auto in1 = nodes[n.in1];
		auto in2 = nodes[n.in2];
		auto in3 = nodes[n.in3];
		return (!is_pi(in1) && !is_pi(in2) && share_two_input_polarity(in1, in2)) || 
				(!is_pi(in1) && !is_pi(in3) && share_two_input_polarity(in1, in3)) ||
				(!is_pi(in2) && !is_pi(in3) && share_two_input_polarity(in2, in3));
	}
                
	// Applies if n has a child in common with one of its children
	inline bool pm_start_ternary_swap(const std::vector<node>& nodes, const node& n) {
		if (is_pi(n))
			return false;
		auto in1 = nodes[n.in1];
		auto in2 = nodes[n.in2];
		auto in3 = nodes[n.in3];
		return (!is_pi(in1) && share_input_polarity(n, in1)) || 
				(!is_pi(in2) && share_input_polarity(n, in2)) || 
				(!is_pi(in3) && share_input_polarity(n, in3));
	}
    
	// Applies if the node has a child that is not a PI
	inline bool pm_start_dist_left_right(const std::vector<node>& nodes, const node& n) {
		if (is_pi(n))
			return false;
		auto in1 = nodes[n.in1];
		auto in2 = nodes[n.in2];
		auto in3 = nodes[n.in3];
		return !is_pi(in1) || !is_pi(in2) || !is_pi(in3);
	}
	
	inline bool pm_start_substitution(const std::vector<node>& nodes, const node& n) {
		return !is_pi(n);
	}
        
	inline bool pm_start_relevance(const std::vector<node>& nodes, const node& n) {
		return !is_pi(n);
	}

};
