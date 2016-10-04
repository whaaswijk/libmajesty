#pragma once

#include <logic_network.h>
#include <truth_table_utils.hpp>
#include <boost/optional.hpp>
#include "minisat/Solver.h"
#include "minisat/SolverTypes.h"

namespace majesty {

	logic_ntk size_optimum_ntk(cirkit::tt& function, unsigned gate_size);
	Minisat::lbool exists_fanin_2_ntk(const cirkit::tt& func, Minisat::Solver&, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_gates);
	logic_ntk extract_fanin_3_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_gates);
	boost::optional<logic_ntk> find_fanin_2_ntk(const cirkit::tt& function, const unsigned nr_gates);
	boost::optional<logic_ntk> find_fanin_3_ntk(const cirkit::tt& function, const unsigned nr_gates);
	
	inline unsigned nr_ordered_tuples(unsigned tuple_size, unsigned bound) {
		if (bound <= tuple_size) {
			return 0;
		} else if (tuple_size == bound + 1) {
			return 1;
		}
		return nr_ordered_tuples(tuple_size, bound - 1) + nr_ordered_tuples(tuple_size - 1, bound - 1);
	}

	static inline Minisat::Var gate_variable(int tt_size, int gate_i, int t) {
		return tt_size * gate_i + t;
	}

	static inline Minisat::Var function_variable(int gate_func_size, int f, int idx, int nr_gate_vars) {
		return f * gate_func_size + idx + nr_gate_vars;
	}
	
	// Tests a primary input's truth table (nonzero) at the specified index
	static inline bool pi_value(const int pi_num, const int tt_idx) {
		const auto num_zeros = (1u << pi_num);
		const auto proj_idx = tt_idx % (2 * num_zeros);
		return proj_idx >= num_zeros;
	}

	void print_fanin_2_solution(const cirkit::tt& func, Minisat::Solver& solver, const unsigned nr_gates);
}