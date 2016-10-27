#pragma once

#include <logic_network.h>
#include <truth_table_utils.hpp>
#include <boost/optional.hpp>
#include "minisat/Solver.h"
#include "minisat/SolverTypes.h"

namespace majesty {

	struct synth_options {
		bool verbose;
		bool use_cegar;
		bool colex_order;
		bool no_triv_ops;
		bool use_all_gates;
		bool exact_nr_svars;
		bool no_reapplication;
	};

	logic_ntk size_optimum_ntk(cirkit::tt& func, synth_options*, const unsigned nr_vars, const unsigned gate_size);
	logic_ntk size_optimum_ntk(uint64_t func, synth_options*, const unsigned nr_vars, const unsigned gate_size);
	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_options*, const unsigned nr_vars, const unsigned gate_size);

	Minisat::lbool exists_fanin_2_ntk(const cirkit::tt& func, Minisat::Solver&, synth_options*, const unsigned nr_vars, const unsigned nr_gates);
	Minisat::lbool exists_fanin_2_ntk(const uint64_t func, Minisat::Solver& solver, synth_options*, const unsigned nr_vars, const unsigned nr_gates);
	Minisat::lbool exists_fanin_2_ntk_ns(const uint64_t func, Minisat::Solver& solver, synth_options*, const unsigned nr_vars, const unsigned nr_gates);

	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates, bool invert);
	logic_ntk extract_fanin_2_ntk(const uint64_t func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk(const uint64_t func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates, bool invert);

	logic_ntk extract_fanin_2_ntk_ns(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk_ns(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates, bool invert);
	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, const Minisat::Solver&, const unsigned nr_vars, const unsigned nr_gates, bool invert);
	
	void print_fanin_2_solution(const cirkit::tt& func, Minisat::Solver& solver, const unsigned nr_vars, const unsigned nr_gates);
	void print_fanin_2_solution(const uint64_t func, Minisat::Solver& solver, const unsigned nr_vars, const unsigned nr_gates);
	void print_fanin_2_solution_ns(const cirkit::tt& func, Minisat::Solver& solver, const unsigned nr_vars, const unsigned nr_gates);
	void print_fanin_2_solution_ns(const uint64_t func, Minisat::Solver& solver, const unsigned nr_vars, const unsigned nr_gates);

	static inline int lbool_to_int(Minisat::lbool b) {
		if (b == l_False) {
			return 0;
		} else if (b == l_True) {
			return 1;
		} else {
			return 2;
		}
	}
	
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

	static inline int gate_variable(int tt_size, int gate_i, int t) {
		return tt_size * gate_i + t;
	}

	static inline int function_variable(int gate_func_size, int f, int idx, int nr_gate_vars) {
		return f * gate_func_size + idx + nr_gate_vars;
	}
	
	// Tests a primary input's truth table (nonzero) at the specified index
	static inline bool pi_value(const int pi_num, const int tt_idx) {
		const auto num_zeros = (1u << pi_num);
		const auto proj_idx = tt_idx % (2 * num_zeros);
		return proj_idx >= num_zeros;
	}
}
