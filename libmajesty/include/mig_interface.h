#pragma once

#include <random>
#include "xmg.h"

namespace majesty {

	enum MoveType {
		MAJ3_PROP = 0,
		INVERTER_PROP,
		SWAP,
		DIST_RIGHT_LEFT,
		SWAP_TERNARY,
		DIST_LEFT_RIGHT,
		MAJ3_XXY,
		MAJ3_XYY
	};

	struct move {
		MoveType type;
		nodeid nodeid1;
		nodeid nodeid2;
		nodeid nodeid3;
	};

	const unsigned NR_UNARY_MOVES = 2;
	const unsigned NR_BINARY_MOVES = 2;
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

	xmg* apply_move(const xmg&, move&);
	float compute_reward(const xmg&, const xmg&);
	std::vector<move> compute_moves(const xmg&);

	xmg* mig_string_decompose(const std::string& truth_table);
	xmg* mig_expression_decompose(unsigned ninputs, const std::string& expr);
	xmg* mig_int_decompose(unsigned ninputs, unsigned truth_table);
	xmg* get_optimum_mig(const xmg&);
	xmg* strash_xmg(const xmg&);
};