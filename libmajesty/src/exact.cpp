#include <exact.h>
#include <convert.h>
#include <vector>

using Minisat::Solver;
using Minisat::Var;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::lbool;
using std::string;
using std::vector;

using namespace cirkit;
using namespace boost;

namespace majesty {


	logic_ntk size_optimum_xmg(tt& func, unsigned gate_size) {
		// TODO: Check if function is constant or variable. If so, return early.


		// Make the function normal if it isn't already
		bool invert = func.test(0);
		if (invert) {
			func = ~func;
		}


		// Else start by checking for networks with increasing nrs. of gates.
		optional<logic_ntk> ntk;
		auto nr_gates = 1u;
		while (true) {
			Solver solver;
			auto network_exists = exists_fanin_2_ntk(func, solver, nr_gates);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk(func, solver, nr_gates);
				break;
			}
			++nr_gates;
		}
		
		if (invert) { // Invert the output node back to the original specification
			func = ~func;
			(*ntk).get_node((*ntk).nnodes() - 1).function = func;
		}

		return move(*ntk);
	}


	static inline Var gate_variable(int tt_size, int gate_i, int t) {
		return tt_size * gate_i + t;
	}

	static inline Var function_variable(int gate_func_size, int f, int idx, int nr_gate_vars) {
		return f * gate_func_size + idx + nr_gate_vars;
	}

	// Tests a primary input's truth table (nonzero) at the specified index
	static inline bool pi_value(const int pi_num, const int tt_idx) {
		const auto num_zeros = (1u << pi_num);
		const auto proj_idx = tt_idx % (2 * num_zeros);
		return proj_idx > num_zeros;
	}
	
	static inline void add_selection_clause(Solver& solver, int num_vars, Var sel_var, int nr_gate_vars,
		int t, int i, int j, int k, bool a, bool b, bool c) {
		static Minisat::vec<Lit> clause_vec;
		clause_vec.clear();

		auto sel_lit = mkLit(sel_var, false);
		auto gate_lit = mkLit(gate_variable(7, i, t), a);
		clause_vec.push(sel_lit);
		clause_vec.push(gate_lit);

		if (j < num_vars) { // It's a primary input, so fixed to a constant
			auto const_val = pi_value(j, t);
			if (const_val != b) {
				return;
			}
		} else {
			auto fanin1_lit = mkLit(gate_variable(7, j, t), b);
			clause_vec.push(fanin1_lit);
		}

		if (k < num_vars) {
			auto const_val = pi_value(k, t);
			if (const_val != c) {
				return;
			}
		} else {
			auto fanin2_lit = mkLit(gate_variable(7, k, t), c);
			clause_vec.push(fanin2_lit);
		}
		
		if (b || c) {
			auto function_lit = mkLit(function_variable(3, i, ((b << 1) | c) - 1, nr_gate_vars), !a);
			clause_vec.push(function_lit);
		}

		solver.addClause(clause_vec);
	}
	
	// Tries to finds a network with the specified size. Result is optional as this may
	// not be possible.
	lbool exists_fanin_2_ntk(const tt& func, Solver& solver, const unsigned nr_gates) {
		const auto num_vars = tt_num_vars(func);
		const auto tt_size = func.size() - 1;

		// Create variables that represent  the gates' truth tables
		auto nr_gate_vars = nr_gates * tt_size;
		for (auto i = 0u; i < nr_gate_vars; i++) {
			// Create gate variable x_it
			solver.newVar();
		}

		// The gate's function constraint variables
		for (auto i = 0u; i < nr_gates * 3; i++) {
			solver.newVar();
		}

		// Add selection (fanin) constraints
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < num_vars + i; j++) {
				for (auto k = j + 1; k < num_vars + i; k++) {
					for (auto t = 0u; t < tt_size; t++) {
						auto sel_var = solver.newVar();
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 0, 0, 1);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 0);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 1);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 0);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 1);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 0);
						add_selection_clause(solver, num_vars, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 1);
					}
				}
			}
		}

		// The final gate's truth table should match the one from the specification
		static Minisat::vec<Lit> spec_clause;
		spec_clause.clear();
		for (auto t = 0u; t < tt_size; t++) {
			spec_clause.push(mkLit(gate_variable(7, nr_gates - 1, t), !func.test(t + 1)));
		}

		if (solver.solve(spec_clause)) {
			return l_True;
		} else {
			return l_False;
		}
	}

	logic_ntk extract_fanin_2_ntk(const tt& func, const Solver& solver, unsigned nr_gates) {
		logic_ntk ntk;

		const auto num_vars = tt_num_vars(func);
		const auto tt_size = 7;
		
		for (auto i = 0u; i < num_vars; i++) {
			ntk.create_input();
		}

		// The last node is the output node
		ntk.create_output(ntk.nnodes() - 1);

		ntk.create_dummy_names();

		return ntk;
	}

}