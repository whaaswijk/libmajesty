#include <exact.h>
#include <convert.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

/*
using Minisat::Solver;
using Minisat::Var;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::lbool;
*/

using namespace cirkit;
using namespace boost;
using namespace std;

namespace majesty {

	static inline void create_variables(sat_solver* solver, synth_spec* spec) {
		auto var_ctr = 0u;

		spec->selection_var_offset = var_ctr;
		unsigned nr_selection_vars = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					++var_ctr;
					++nr_selection_vars;
				}
			}
		}
		spec->nr_selection_vars = nr_selection_vars;

		spec->gate_var_offset = var_ctr;
		const auto nr_gate_vars = spec->nr_gates * 3;
		spec->nr_gate_vars = nr_gate_vars;
		var_ctr += nr_gate_vars;

		spec->simulation_var_offset = var_ctr;
		const auto nr_simulation_vars = spec->nr_gates * spec->tt_size;
		spec->nr_simulation_vars = nr_simulation_vars;
		var_ctr += nr_simulation_vars;
		
		sat_solver_setnvars(solver, var_ctr);
	}

	static inline void create_variables_ns(sat_solver* solver, synth_spec* spec) {
		auto var_ctr = 0u;

		spec->selection_var_offset = var_ctr;
		const auto nr_selection_vars = spec->nr_gates * spec->nr_vars + (((spec->nr_gates - 1) * spec->nr_gates) / 2);
		var_ctr += nr_selection_vars;
		spec->nr_selection_vars = nr_selection_vars;

		spec->gate_var_offset = var_ctr;
		const auto nr_gate_vars = spec->nr_gates * 3;
		spec->nr_gate_vars = nr_gate_vars;
		var_ctr += nr_gate_vars;

		spec->simulation_var_offset = var_ctr;
		const auto nr_simulation_vars = spec->nr_gates * spec->tt_size;
		spec->nr_simulation_vars = nr_simulation_vars;
		var_ctr += nr_simulation_vars;
		
		sat_solver_setnvars(solver, var_ctr);
	}

	unsigned optimum_ntk_size(uint64_t func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.

		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		spec->nr_gates = 1u;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		auto solver = sat_solver_new();
		while (true) {
			sat_solver_restart(solver);
			auto network_exists = exists_fanin_2_ntk(func, solver, spec);
			if (network_exists == l_True) {
				if (spec->verbose) {
					print_fanin_2_solution(func, solver, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		sat_solver_delete(solver);
		
		return spec->nr_gates;
	}

	unsigned optimum_ntk_size_ns(uint64_t func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.

		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		spec->nr_gates = 1u;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		auto solver = sat_solver_new();
		while (true) {
			sat_solver_restart(solver);
			auto network_exists = exists_fanin_2_ntk_ns(func, solver, spec);
			if (network_exists == l_True) {
				if (spec->verbose) {
					print_fanin_2_solution_ns(func, solver, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		sat_solver_delete(solver);
		
		return spec->nr_gates;
	}

	logic_ntk size_optimum_ntk(uint64_t func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.


		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		logic_ntk ntk;
		spec->nr_gates = 1u;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		auto solver = sat_solver_new();
		while (true) {
			sat_solver_restart(solver);
			auto network_exists = exists_fanin_2_ntk(func, solver, spec);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk(func, solver, spec, invert);
				if (spec->verbose) {
					print_fanin_2_solution(func, solver, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		sat_solver_delete(solver);
		
		return std::move(ntk);
	}

	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.


		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		logic_ntk ntk;
		spec->nr_gates = 1u;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		auto solver = sat_solver_new();
		while (true) {
			sat_solver_restart(solver);
			auto network_exists = exists_fanin_2_ntk_ns(func, solver, spec);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk_ns(func, solver, spec, invert);
				if (spec->verbose) {
					print_fanin_2_solution(func, solver, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		sat_solver_delete(solver);
		
		return std::move(ntk);
	}
	

	static inline void add_selection_clause(sat_solver* solver, synth_spec* spec,
		int t, int i, int j, int k, bool a, bool b, bool c) {
		static lit plits[5];
		unsigned ctr = 0u;

		if (j < spec->nr_vars) { // It's a primary input, so fixed to a constant
			auto const_val = pi_value(j, t+1);
			if (const_val != b) {
				return;
			}
		} else {
			plits[ctr++] = Abc_Var2Lit(simulation_variable(j - spec->nr_vars, t, spec), b);
		}

		if (k < spec->nr_vars) {
			auto const_val = pi_value(k, t+1);
			if (const_val != c) {
				return;
			}
		} else {
			plits[ctr++] = Abc_Var2Lit(simulation_variable(k - spec->nr_vars, t, spec), c);
		}

		auto sel_var = selection_variable(spec, i, j, k);
		auto sel_lit = Abc_Var2Lit(sel_var, true);
		auto gate_lit = Abc_Var2Lit(simulation_variable(i, t, spec), a);
		plits[ctr++] = sel_lit;
		plits[ctr++] = gate_lit;
		
		if (b | c) {
			plits[ctr++] = Abc_Var2Lit(gate_variable(i, ((c << 1) | b) - 1, spec), !a);
		}
		sat_solver_addclause(solver, plits, plits + ctr);
	}

	static inline void add_selection_clause_ns(sat_solver* solver, synth_spec* spec, 
		int t, int i, int j, int k, bool a, bool b, bool c) {
		static lit plits[6];
		unsigned ctr = 0u;

		if (j < spec->nr_vars) { // It's a primary input, so fixed to a constant
			auto const_val = pi_value(j, t+1);
			if (const_val != b) {
				return;
			}
		} else {
			plits[ctr++] = Abc_Var2Lit(simulation_variable(j - spec->nr_vars, t, spec), b);
		}

		if (k < spec->nr_vars) {
			auto const_val = pi_value(k, t+1);
			if (const_val != c) {
				return;
			}
		} else {
			plits[ctr++] = Abc_Var2Lit(simulation_variable(k - spec->nr_vars, t, spec), c);
		}

		auto sel_var1 = selection_variable_ns(spec, i, j);
		auto sel_var2 = selection_variable_ns(spec, i, k);
		auto sel_lit1 = Abc_Var2Lit(sel_var1, true);
		auto sel_lit2 = Abc_Var2Lit(sel_var2, true);
		auto gate_lit = Abc_Var2Lit(simulation_variable(i, t, spec), a);
		plits[ctr++] = sel_lit1;
		plits[ctr++] = sel_lit2;
		plits[ctr++] = gate_lit;
		
		if (b | c) {
			plits[ctr++] = Abc_Var2Lit(gate_variable(i, ((c << 1) | b) - 1, spec), !a);
		}
		sat_solver_addclause(solver, plits, plits + ctr);
	}

	static inline lbool cegar_solve(const uint64_t func, sat_solver* solver, synth_spec* spec, Vec_Int_t* vlits) {
		Vec_IntClear(vlits);

		while (true) {
			auto res = sat_solver_solve(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits), 0, 0, 0, 0);
			if (res == l_False) {
				return res;
			}

			// Extract the network, simulate it, and find the first bit where it's different from the specification
			auto ntk = extract_fanin_2_ntk(func, solver, spec);
			auto ntk_func = ntk.simulate();

			bool found_solution = true;
			// Check if the solution found matches the specification
			for (auto t = 0u; t < spec->tt_size; t++) {
				const auto ntk_val = ntk_func[t + 1];
				const bool spec_val = (func >> (t + 1)) & 1;
				if (ntk_val != spec_val) {
					// Constrain the solver further by adding an additional constraint for this truth table row
					const auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
					Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !spec_val));
					for (auto i = 0u; i < spec->nr_gates; i++) {
						for (auto j = 0u; j < spec->nr_vars + i; j++) {
							for (auto k = j + 1; k < spec->nr_vars + i; k++) {
								add_selection_clause(solver, spec, t, i, j, k, 0, 0, 1);
								add_selection_clause(solver, spec, t, i, j, k, 0, 1, 0);
								add_selection_clause(solver, spec, t, i, j, k, 0, 1, 1);
								add_selection_clause(solver, spec, t, i, j, k, 1, 0, 0);
								add_selection_clause(solver, spec, t, i, j, k, 1, 0, 1);
								add_selection_clause(solver, spec, t, i, j, k, 1, 1, 0);
								add_selection_clause(solver, spec, t, i, j, k, 1, 1, 1);
							}
						}
					}
					found_solution = false;
					break;
				}
			}
			if (found_solution) {
				break;
			}
		}
		return l_True;
	}

	static inline lbool cegar_solve_ns(const uint64_t func, sat_solver* solver, synth_spec* spec, Vec_Int_t* vlits) {
		Vec_IntClear(vlits);

		while (true) {
			auto res = sat_solver_solve(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits), 0, 0, 0, 0);
			if (res == l_False) {
				return res;
			}

			// Extract the network, simulate it, and find the first bit where it's different from the specification
			auto ntk = extract_fanin_2_ntk_ns(func, solver, spec);
			auto ntk_func = ntk.simulate();

			bool found_solution = true;
			// Check if the solution found matches the specification
			for (auto t = 0u; t < spec->tt_size; t++) {
				const auto ntk_val = ntk_func[t + 1];
				const bool spec_val = (func >> (t + 1)) & 1;
				if (ntk_val != spec_val) {
					// Constrain the solver further by adding an additional constraint for this truth table row
					const auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
					Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !spec_val));
					for (auto i = 0u; i < spec->nr_gates; i++) {
						for (auto j = 0u; j < spec->nr_vars + i; j++) {
							for (auto k = j + 1; k < spec->nr_vars + i; k++) {
								add_selection_clause_ns(solver, spec, t, i, j, k, 0, 0, 1);
								add_selection_clause_ns(solver, spec, t, i, j, k, 0, 1, 0);
								add_selection_clause_ns(solver, spec, t, i, j, k, 0, 1, 1);
								add_selection_clause_ns(solver, spec, t, i, j, k, 1, 0, 0);
								add_selection_clause_ns(solver, spec, t, i, j, k, 1, 0, 1);
								add_selection_clause_ns(solver, spec, t, i, j, k, 1, 1, 0);
								add_selection_clause_ns(solver, spec, t, i, j, k, 1, 1, 1);
							}
						}
					}
					found_solution = false;
					break;
				}
			}
			if (found_solution) {
				break;
			}
		}
		return l_True;
	}

	lbool exists_fanin_2_ntk(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		static lit plits[3];
		
		create_variables(solver, spec);

		// The gate's function constraint variables
		if (spec->use_no_triv_ops) {
			for (auto i = 0u; i < spec->nr_gates; i++) {
				const auto func_var0 = gate_variable(i, 0, spec);
				const auto func_var1 = gate_variable(i, 1, spec);
				const auto func_var2 = gate_variable(i, 2, spec);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 0);
				sat_solver_addclause(solver, plits, plits + 3);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 1);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				sat_solver_addclause(solver, plits, plits + 3);
				
				plits[0] = Abc_Var2Lit(func_var0, 1);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				sat_solver_addclause(solver, plits, plits + 3);
			}
		}

		// Add selection (fanin) constraints
		Vec_Int_t * vlits = Vec_IntAlloc(spec->nr_selection_vars);
		for (auto i = 0u; i < spec->nr_gates; i++) {
			Vec_IntClear(vlits);
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					auto sel_var = selection_variable(spec, i, j, k);
					Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					if (!spec->use_cegar) {
						for (auto t = 0u; t < spec->tt_size; t++) {
							add_selection_clause(solver, spec, t, i, j, k, 0, 0, 1);
							add_selection_clause(solver, spec, t, i, j, k, 0, 1, 0);
							add_selection_clause(solver, spec, t, i, j, k, 0, 1, 1);
							add_selection_clause(solver, spec, t, i, j, k, 1, 0, 0);
							add_selection_clause(solver, spec, t, i, j, k, 1, 0, 1);
							add_selection_clause(solver, spec, t, i, j, k, 1, 1, 0);
							add_selection_clause(solver, spec, t, i, j, k, 1, 1, 1);
						}
					}
				}
			}
			sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
			// Allow at most one selection variable to be true at a time
			if (spec->use_exact_nr_svars || spec->use_cegar) {
				for (auto u = 0u; u < Vec_IntSize(vlits) - 1; u++) {
					auto svar1 = Abc_Lit2Var(Vec_IntEntry(vlits, u));
					for (auto v = u + 1; v < Vec_IntSize(vlits); v++) {
						auto svar2 = Abc_Lit2Var(Vec_IntEntry(vlits, v));
						plits[0] = Abc_Var2Lit(svar1, true);
						plits[1] = Abc_Var2Lit(svar2, true);
						sat_solver_addclause(solver, plits, plits + 2);
					}
				}
			}
		}

		// Check for co-lexicographic order: given gates i and i+1 with fanin (x, y) and (x', y') respectively,
		// we require y < y' OR y = y' AND x <= x'. This only if i is not a fanin of i+1.
		if (spec->use_colex_order) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto j = 0u; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto jp = 0u; jp < spec->nr_vars + i; jp++) {
							for (auto kp = jp + 1; kp < spec->nr_vars + i; kp++) {
								if (k == kp && j > jp || k > kp) {
									auto sivar = selection_variable(spec, i, j, k);
									auto sipvar = selection_variable(spec, i+1, jp, kp);
									plits[0] = Abc_Var2Lit(sivar, true);
									plits[1] = Abc_Var2Lit(sipvar, true);
									sat_solver_addclause(solver, plits, plits + 2);
								}
							}
						}
					}
				}
			}
		}

		// Enforce that a solution uses all gates
		if (spec->use_all_gates) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				Vec_IntClear(vlits);
				for (auto ip = i + 1; ip < spec->nr_gates; ip++) {
					for (auto j = 0u; j < i + spec->nr_vars; j++) {
						auto sel_var = selection_variable(spec, ip, j, i + spec->nr_vars);
						Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					}
					for (auto j = i + spec->nr_vars + 1; j < ip + spec->nr_vars; j++) {
						auto sel_var = selection_variable(spec, ip, i + spec->nr_vars, j);
						Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					}
				}
				sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
			}
		}

		// Do not allow reapplication of operators
		if (spec->use_no_reapplication) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto ip = i + 1; ip < spec->nr_gates; ip++) {
					for (auto j = 0u; j < spec->nr_vars + i; j++) {
						for (auto k = j + 1; k < spec->nr_vars + i; k++) {
							auto sivar = selection_variable(spec, i, j, k);
							plits[0] = Abc_Var2Lit(sivar, true);
							auto sipvar = selection_variable(spec, i + 1, j, i + spec->nr_vars);
							plits[1] = Abc_Var2Lit(sipvar, true);
							sat_solver_addclause(solver, plits, plits + 2);
							sipvar = selection_variable(spec, i + 1, k, i + spec->nr_vars);
							plits[1] = Abc_Var2Lit(sipvar, true);
							sat_solver_addclause(solver, plits, plits + 2);
						}
					}
				}
			}
		}


		if (spec->use_cegar) {
			auto res = cegar_solve(func, solver, spec, vlits);
			Vec_IntFree(vlits);
			return res;
		} else {
			// The final gate's truth table should match the one from the specification
			Vec_IntClear(vlits);
			for (auto t = 0u; t < spec->tt_size; t++) {
				auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
				Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !((func >> (t + 1)) & 1)));
			}

			auto res = sat_solver_solve(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits), 0, 0, 0, 0);
			Vec_IntFree(vlits);
			return res;
		}
	}
	
	// Tries to find a network using the new selection variable implementation
	lbool exists_fanin_2_ntk_ns(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		static lit plits[4];
		
		create_variables_ns(solver, spec);

		// The gate's function constraint variables
		if (spec->use_no_triv_ops) {
			for (auto i = 0u; i < spec->nr_gates; i++) {
				const auto func_var0 = gate_variable(i, 0, spec);
				const auto func_var1 = gate_variable(i, 1, spec);
				const auto func_var2 = gate_variable(i, 2, spec);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 0);
				sat_solver_addclause(solver, plits, plits + 3);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 1);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				sat_solver_addclause(solver, plits, plits + 3);
				
				plits[0] = Abc_Var2Lit(func_var0, 1);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				sat_solver_addclause(solver, plits, plits + 3);
			}
		}

		// Add selection (fanin) constraints
		Vec_Int_t * vlits = Vec_IntAlloc(spec->nr_selection_vars);
		const auto k = 2;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			const auto n = spec->nr_vars + i;
			// Exactly 2 selection vars should be true, since every gate has fanin 2
			// We do this by first preventing more than k vars from being true by selecting n choose k+1 subsets
			if (n > 2 && (spec->use_exact_nr_svars || spec->use_cegar)) {
				vector<bool> v(n);
				std::fill(v.end() - (k + 1), v.end(), true);
				//cout << "nr spec vars: " << spec->nr_vars << ", i: " << i << endl;
				do {
					Vec_IntClear(vlits);
					for (int j = 0; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), true));
						}
					}
					sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}
			// Next we prevent more than n-k vars from being false by selecting n choose n-k+1 subsets
			{
				vector<bool> v(n);
				std::fill(v.end() - (n - k + 1), v.end(), true);
				do {
					Vec_IntClear(vlits);
					for (int j = 0; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), false));
						}
					}
					sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}

			// Finally, we add the constaint clauses
			if (!spec->use_cegar) {
				for (auto j = 0; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto t = 0u; t < spec->tt_size; t++) {
							add_selection_clause_ns(solver, spec, t, i, j, k, 0, 0, 1);
							add_selection_clause_ns(solver, spec, t, i, j, k, 0, 1, 0);
							add_selection_clause_ns(solver, spec, t, i, j, k, 0, 1, 1);
							add_selection_clause_ns(solver, spec, t, i, j, k, 1, 0, 0);
							add_selection_clause_ns(solver, spec, t, i, j, k, 1, 0, 1);
							add_selection_clause_ns(solver, spec, t, i, j, k, 1, 1, 0);
							add_selection_clause_ns(solver, spec, t, i, j, k, 1, 1, 1);
						}
					}
				}
			}
		}

		// Enforce that a solution uses all gates
		if (spec->use_all_gates) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				Vec_IntClear(vlits);
				for (auto ip = i + 1; ip < spec->nr_gates; ip++) {
					auto sel_var = selection_variable_ns(spec, ip, i + spec->nr_vars);
					Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
				}
				sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
			}
		}

		// Check for co-lexicographic order: given gates i and i+1 with fanin (x, y) and (x', y') respectively,
		// we require y < y' OR y = y' AND x <= x'. This only if i is not a fanin of i+1.
		if (spec->use_colex_order) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto j = 0u; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto jp = 0u; jp < spec->nr_vars + i; jp++) {
							for (auto kp = jp + 1; kp < spec->nr_vars + i; kp++) {
								if (k == kp && j > jp || k > kp) {
									auto sivar1 = selection_variable_ns(spec, i, j);
									auto sivar2 = selection_variable_ns(spec, i, k);
									auto sipvar1 = selection_variable_ns(spec, i+1, jp);
									auto sipvar2 = selection_variable_ns(spec, i+1, kp);
									plits[0] = Abc_Var2Lit(sivar1, true);
									plits[1] = Abc_Var2Lit(sivar2, true);
									plits[2] = Abc_Var2Lit(sipvar1, true);
									plits[3] = Abc_Var2Lit(sipvar2, true);
									sat_solver_addclause(solver, plits, plits + 4);
								}
							}
						}
					}
				}
			}
		}

		// Do not allow reapplication of operators
		if (spec->use_no_reapplication) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto ip = i + 1; ip < spec->nr_gates; ip++) {
					for (auto j = 0u; j < spec->nr_vars + i; j++) {
						for (auto k = j + 1; k < spec->nr_vars + i; k++) {
							auto sivar1 = selection_variable_ns(spec, i, j);
							auto sivar2 = selection_variable_ns(spec, i, k);
							plits[0] = Abc_Var2Lit(sivar1, true);
							plits[1] = Abc_Var2Lit(sivar2, true);
							auto sipvar1 = selection_variable_ns(spec, ip, i + spec->nr_vars);
							auto sipvar2 = selection_variable_ns(spec, ip, j);
							plits[2] = Abc_Var2Lit(sipvar1, true);
							plits[3] = Abc_Var2Lit(sipvar2, true);
							sat_solver_addclause(solver, plits, plits + 4);
							sipvar2 = selection_variable_ns(spec, ip, k);
							plits[3] = Abc_Var2Lit(sipvar2, true);
							sat_solver_addclause(solver, plits, plits + 4);
						}
					}
				}
			}
		}
		
		if (spec->use_cegar) {
			auto res = cegar_solve_ns(func, solver, spec, vlits);
			Vec_IntFree(vlits);
			return res;
		} else {
			// The final gate's truth table should match the one from the specification
			Vec_IntClear(vlits);
			for (auto t = 0u; t < spec->tt_size; t++) {
				auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
				Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !((func >> (t + 1)) & 1)));
			}
			auto res = sat_solver_solve(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits), 0, 0, 0, 0);
			Vec_IntFree(vlits);
			return res;
		}
	}
	
    logic_ntk extract_fanin_2_ntk(const tt& func, sat_solver* solver, synth_spec* spec) {
        return extract_fanin_2_ntk(func, solver, spec, false);
    }

	logic_ntk extract_fanin_2_ntk(const tt& func, sat_solver* solver, synth_spec* spec, bool invert) {
		logic_ntk ntk;

		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		auto sel_var_ctr = 0u;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					if (sat_solver_var_value(solver, spec->selection_var_offset + sel_var_ctr)) {
						// Gate i has gates j and k as fanin
						fanin[0] = j;
						fanin[1] = k;
						nodefunc.set(0, 0);
						nodefunc.set(1, sat_solver_var_value(solver, gate_variable(i, 0, spec)));
						nodefunc.set(2, sat_solver_var_value(solver, gate_variable(i, 1, spec)));
						nodefunc.set(3, sat_solver_var_value(solver, gate_variable(i, 2, spec)));
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
	
	logic_ntk extract_fanin_2_ntk(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		return extract_fanin_2_ntk(func, solver, spec, false);
	}

	logic_ntk extract_fanin_2_ntk(const uint64_t func, sat_solver* solver, synth_spec* spec, bool invert) {
		logic_ntk ntk;
		
		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					if (sat_solver_var_value(solver, spec->selection_var_offset + selection_var_ctr++)) {
						// Gate i has gates j and k as fanin
						fanin[0] = j;
						fanin[1] = k;
						nodefunc.set(0, 0);
						nodefunc.set(1, sat_solver_var_value(solver, gate_variable(i, 0, spec)));
						nodefunc.set(2, sat_solver_var_value(solver, gate_variable(i, 1, spec)));
						nodefunc.set(3, sat_solver_var_value(solver, gate_variable(i, 2, spec)));
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

	logic_ntk extract_fanin_2_ntk_ns(const tt& func, sat_solver* solver, synth_spec* spec) {
		return extract_fanin_2_ntk_ns(func, solver, spec, false);
	}

	logic_ntk extract_fanin_2_ntk_ns(const tt& func, sat_solver* solver, synth_spec* spec, bool invert) {
		logic_ntk ntk;

		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				if (sat_solver_var_value(solver, spec->selection_var_offset + selection_var_ctr++)) {
					if (inputs_found == 0u) {
						in1 = j;
					} else {
						in2 = j;
					}
					++inputs_found;
					if (inputs_found == 2) {
						fanin[0] = in1;
						fanin[1] = in2;
						nodefunc.set(0, 0);
						nodefunc.set(1, sat_solver_var_value(solver, gate_variable(i, 0, spec)));
						nodefunc.set(2, sat_solver_var_value(solver, gate_variable(i, 1, spec)));
						nodefunc.set(3, sat_solver_var_value(solver, gate_variable(i, 2, spec)));
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

	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		return extract_fanin_2_ntk_ns(func, solver, spec, false);
	}

	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, sat_solver* solver, synth_spec* spec, bool invert) {
		logic_ntk ntk;
		
		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				if (sat_solver_var_value(solver, spec->selection_var_offset + selection_var_ctr++)) {
					if (inputs_found == 0u) {
						in1 = j;
					} else {
						in2 = j;
					}
					++inputs_found;
					if (inputs_found == 2) {
						fanin[0] = in1;
						fanin[1] = in2;
						nodefunc.set(0, 0);
						nodefunc.set(1, sat_solver_var_value(solver, gate_variable(i, 0, spec)));
						nodefunc.set(2, sat_solver_var_value(solver, gate_variable(i, 1, spec)));
						nodefunc.set(3, sat_solver_var_value(solver, gate_variable(i, 2, spec)));
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

	void print_fanin_2_solution(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		tt functt((1u << spec->nr_vars), func);
		print_fanin_2_solution(functt, solver, spec);
	}

	void print_fanin_2_solution(const tt& func, sat_solver* solver, synth_spec* spec) {
		cout << "Solution for function " << to_string(func) << endl;

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto t = 0u; t < spec->tt_size; t++) {
				cout << "x_" << i << "_" << t << ": " << sat_solver_var_value(solver, simulation_variable(i, t, spec)) << endl;
			}
		}

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < 3; j++) {
				cout << "f_" << i << "_" << j << ": " << sat_solver_var_value(solver, gate_variable(i, j, spec)) << endl;
			}
		}

		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					cout << "s_" << i << "_" << j << "_" << k << ": " << 
						sat_solver_var_value(solver, spec->selection_var_offset + selection_var_ctr++) << endl;
				}
			}
		}
	}
	
	void print_fanin_2_solution_ns(const uint64_t func, sat_solver* solver, synth_spec* spec) {
		tt functt((1u << spec->nr_vars), func);
		print_fanin_2_solution_ns(functt, solver, spec);
	}

	void print_fanin_2_solution_ns(const tt& func, sat_solver* solver, synth_spec* spec) {
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto t = 0u; t < spec->tt_size; t++) {
				cout << "x_" << i << "_" << t << ": " << sat_solver_var_value(solver, simulation_variable(i, t, spec)) << endl;
			}
		}

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < 3; j++) {
				cout << "f_" << i << "_" << j << ": " << sat_solver_var_value(solver, gate_variable(i, j, spec)) << endl;
			}
		}

		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				cout << "s_" << i << "_" << j << ": " << sat_solver_var_value(solver, spec->selection_var_offset + selection_var_ctr++) << endl;
			}
		}
	}

}
