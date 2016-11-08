#pragma once

#include <logic_network.h>
#include <truth_table_utils.hpp>
#include <boost/optional.hpp>
#include <tuple>
//#include "minisat/Solver.h"
//#include "minisat/SolverTypes.h"
extern "C" {
#include <base/abc/abc.h>
#include <misc/vec/vecInt.h>
#include <misc/vec/vecPtr.h>
#include <sat/bsat/satSolver.h>
}
	
namespace majesty {

	struct synth_spec {
		bool verbose;
		bool use_cegar;
		bool use_colex_order;
		bool use_no_triv_ops;
		bool use_all_gates;
		bool use_exact_nr_svars;
		bool use_no_reapplication;
		bool use_colex_functions;
		unsigned nr_vars;
		unsigned tt_size;
		unsigned nr_gates;
		unsigned selection_var_offset;
		unsigned nr_selection_vars;
		unsigned gate_var_offset;
		unsigned nr_gate_vars;
		unsigned simulation_var_offset;
		unsigned nr_simulation_vars;
	};
	
	logic_ntk size_optimum_ntk(cirkit::tt& func, synth_spec*);
	logic_ntk size_optimum_ntk(uint64_t func, synth_spec*);
	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_spec*);

	unsigned optimum_ntk_size(uint64_t func, synth_spec*);
	unsigned optimum_ntk_size_ns(uint64_t func, synth_spec*);

	lbool exists_fanin_2_ntk(const cirkit::tt& func, sat_solver*, synth_spec*);
	lbool exists_fanin_2_ntk(const uint64_t func, sat_solver*, synth_spec*);
	lbool exists_fanin_2_ntk_ns(const uint64_t func, sat_solver*, synth_spec*);

	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, sat_solver*, synth_spec*);
	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, sat_solver*, synth_spec*, bool invert);
	logic_ntk extract_fanin_2_ntk(const uint64_t func, sat_solver*, synth_spec*);
	logic_ntk extract_fanin_2_ntk(const uint64_t func, sat_solver*, synth_spec*, bool invert);
	
	logic_ntk extract_fanin_2_ntk_ns(const cirkit::tt& func, sat_solver*, synth_spec*);
	logic_ntk extract_fanin_2_ntk_ns(const cirkit::tt& func, sat_solver*, synth_spec*, bool invert);
	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, sat_solver*, synth_spec*);
	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, sat_solver*, synth_spec*, bool invert);
	
	void print_fanin_2_solution(const cirkit::tt& func, sat_solver*, synth_spec*);
	void print_fanin_2_solution(const uint64_t func, sat_solver*, synth_spec*);
	void print_fanin_2_solution_ns(const cirkit::tt& func, sat_solver*, synth_spec*);
	void print_fanin_2_solution_ns(const uint64_t func, sat_solver*, synth_spec*);

	/*
	Minisat::lbool exists_fanin_2_ntk(const cirkit::tt& func, Minisat::Solver&, synth_spec*, const unsigned nr_vars, const unsigned nr_gates);
	Minisat::lbool exists_fanin_2_ntk(const uint64_t func, Minisat::Solver& solver, synth_spec*, const unsigned nr_vars, const unsigned nr_gates);
	Minisat::lbool exists_fanin_2_ntk_ns(const uint64_t func, Minisat::Solver& solver, synth_spec*, const unsigned nr_vars, const unsigned nr_gates);

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
	*/

	static inline int selection_variable(const synth_spec* spec, unsigned gate_i, unsigned fanin_j, unsigned fanin_k) {
		auto ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					if (i == gate_i && j == fanin_j && k == fanin_k) {
						break;
					}
					++ctr;
				}
			}
		}
		return spec->selection_var_offset + ctr;
	}

	static inline int selection_variable_ns(const synth_spec* spec, unsigned gate_i, unsigned fanin_j) {
		return spec->selection_var_offset + fanin_j + gate_i * spec->nr_vars + (((gate_i - 1) * gate_i) / 2);
	}
	
	static inline int simulation_variable(int gate_i, int t, const synth_spec* spec) {
		return spec->tt_size * gate_i + t + spec->simulation_var_offset;
	}

	static inline int gate_variable(int gate, int idx, const synth_spec* spec) {
		return gate * 3 + idx + spec->gate_var_offset;
	}
	
	// Tests a primary input's truth table (nonzero) at the specified index
	static inline bool pi_value(const int pi_num, const int tt_idx) {
		const auto num_zeros = (1u << pi_num);
		const auto proj_idx = tt_idx % (2 * num_zeros);
		return proj_idx >= num_zeros;
	}
}
