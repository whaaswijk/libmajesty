#include <exact.h>
#include <convert.h>
#include <vector>
#include <iostream>

using Minisat::Solver;
using Minisat::Var;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::lbool;

using namespace cirkit;
using namespace boost;
using namespace std;

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

		return std::move(*ntk);
	}

	static inline void add_selection_clause(Solver& solver, int num_vars, int tt_size, Var sel_var, int nr_gate_vars,
		int t, int i, int j, int k, bool a, bool b, bool c) {
		static Minisat::vec<Lit> clause_vec;
		clause_vec.clear();

		if (j < num_vars) { // It's a primary input, so fixed to a constant
			auto const_val = pi_value(j, t+1);
			if (const_val != b) {
				return;
			}
		} else {
			auto fanin1_lit = mkLit(gate_variable(tt_size, j - num_vars, t), b);
			clause_vec.push(fanin1_lit);
		}

		if (k < num_vars) {
			auto const_val = pi_value(k, t+1);
			if (const_val != c) {
				return;
			}
		} else {
			auto fanin2_lit = mkLit(gate_variable(tt_size, k - num_vars, t), c);
			clause_vec.push(fanin2_lit);
		}

		auto sel_lit = mkLit(sel_var, true);
		auto gate_lit = mkLit(gate_variable(tt_size, i, t), a);
		clause_vec.push(sel_lit);
		clause_vec.push(gate_lit);
		
		if (b | c) {
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
			//cout << "adding gate tt var " << i << endl;
			solver.newVar();
		}

		// The gate's function constraint variables
		for (auto i = 0u; i < nr_gates * 3; i++) {
			//cout << "adding function var " << i << endl;
			solver.newVar();
		}

		// Add selection (fanin) constraints
		static Minisat::vec<Lit> sel_vars_clause;
		for (auto i = 0u; i < nr_gates; i++) {
			sel_vars_clause.clear();
			for (auto j = 0u; j < num_vars + i; j++) {
				for (auto k = j + 1; k < num_vars + i; k++) {
					auto sel_var = solver.newVar();
					sel_vars_clause.push(mkLit(sel_var, false));
					for (auto t = 0u; t < tt_size; t++) {
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 0, 1);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 0);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 1);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 0);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 1);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 0);
						add_selection_clause(solver, num_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 1);
					}
				}
			}
			solver.addClause(sel_vars_clause);
		}
		
		// The final gate's truth table should match the one from the specification
		static Minisat::vec<Lit> spec_clause;
		spec_clause.clear();
		for (auto t = 0u; t < tt_size; t++) {
			auto gate_var = gate_variable(tt_size, nr_gates - 1, t);
			//cout << "setting x_" << nr_gates - 1 << "_" << t << "(" << gate_var << "): " << func.test(t + 1) << endl;
			spec_clause.push(mkLit(gate_var, !func.test(t + 1)));
		}

		if (solver.solve(spec_clause)) {
			return l_True;
		} else {
			return l_False;
		}
	}
	
    logic_ntk extract_fanin_2_ntk(const tt& func, const Solver& solver, unsigned nr_gates) {
        return extract_fanin_2_ntk(func, solver, nr_gates, false);
    }

	logic_ntk extract_fanin_2_ntk(const tt& func, const Solver& solver, unsigned nr_gates, bool invert) {
		logic_ntk ntk;

		const auto num_vars = tt_num_vars(func);
		const auto tt_size = func.size() - 1;
		
		for (auto i = 0u; i < num_vars; i++) {
			ntk.create_input();
		}

		const auto nr_gate_vars = nr_gates * tt_size;
		const auto nr_function_vars = nr_gates * 3;
		auto var_offset = nr_gate_vars + nr_function_vars;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < num_vars + i; j++) {
				for (auto k = j + 1; k < num_vars + i; k++) {
					if (solver.modelValue(mkLit(var_offset++, false)) == l_True) {
						// Gate i has gates j and k as fanin
						fanin[0] = j;
						fanin[1] = k;
						nodefunc.set(0, 0);
						nodefunc.set(1, lbool_to_int(solver.modelValue(function_variable(3, i, 0, nr_gate_vars))));
						nodefunc.set(2, lbool_to_int(solver.modelValue(function_variable(3, i, 1, nr_gate_vars))));
						nodefunc.set(3, lbool_to_int(solver.modelValue(function_variable(3, i, 2, nr_gate_vars))));
						ntk.create_node(fanin, nodefunc);
					}
				}
			}
		}

		// The last node is the output node
        if (invert) {
            auto& outnode = ntk.get_node(ntk.nnodes() - 1);
            outnode.function = ~outnode.function;
        }
		ntk.create_output(ntk.nnodes() - 1);

		return ntk;
	}

	void print_fanin_2_solution(const tt& func, Solver& solver, const unsigned nr_gates) {
		const auto num_vars = tt_num_vars(func);
		const auto tt_size = func.size() - 1;

		auto nr_gate_vars = nr_gates * tt_size;
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto t = 0u; t < tt_size; t++) {
				cout << "x_" << i << "_" << t << ": " << lbool_to_int(solver.modelValue(gate_variable(tt_size, i, t))) << endl;
			}
		}

		auto nr_function_vars = nr_gates * 3;
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < 3; j++) {
				cout << "f_" << i << "_" << j << ": " << lbool_to_int(solver.modelValue(function_variable(3, i, j, nr_gate_vars))) << endl;
			}
		}

		auto var_offset = nr_gate_vars + nr_function_vars;
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < num_vars + i; j++) {
				for (auto k = j + 1; k < num_vars + i; k++) {
					cout << "s_" << i << "_" << j << "_" << k << ": " << lbool_to_int(solver.modelValue(mkLit(var_offset++, false))) << endl;
				}
			}
		}
	}

}
