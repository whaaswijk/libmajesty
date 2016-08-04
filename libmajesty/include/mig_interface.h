#pragma once

#include <random>
#include "xmg.h"

namespace majesty {

	enum UNARY_MOVE {
		MAJ3_LEFT_RIGHT = 0,
		INVERTER_PROP
	};

	enum BINARY_MOVE {
		MAJ3_XXY = 0,
		MAJ3_XYY,
		SWAP
	};

	const unsigned NR_UNARY_MOVES = 2;
	const unsigned NR_BINARY_MOVES = 3;
	const unsigned NR_EDGE_TYPES = 2;
	const unsigned DEFAULT_SEED = 100;


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

		unsigned get_nr_unary_moves() const {
			return NR_UNARY_MOVES;
		}

		unsigned get_nr_binary_moves() const {
			return NR_BINARY_MOVES;
		}

		unsigned get_nr_edge_types() const {
			return NR_EDGE_TYPES;
		}

		unsigned get_seed() {
			return seed;
		}

		void set_seed(unsigned seed) {
			rng.seed(seed);
		}

		xmg create_random_graph(unsigned ninputs, unsigned nnodes);
	};

	xmg* apply_unary_move(const xmg&, UNARY_MOVE, nodeid);
	xmg* apply_binary_move(const xmg&, BINARY_MOVE, nodeid, nodeid);
	void mig_to_array(const xmg&, nodeid*);
	float compute_reward(const xmg&, const xmg&);
};