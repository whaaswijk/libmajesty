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

// function has to live in the std namespace 
// so that it is picked up by argument-dependent name lookup (ADL).
namespace std {
	namespace
	{

		// Code from boost
		// Reciprocal of the golden ratio helps spread entropy
		//     and handles duplicates.
		// See Mike Seymour in magic-numbers-in-boosthash-combine:
		//     http://stackoverflow.com/questions/4948780

		template <class T>
		inline void hash_combine(std::size_t& seed, T const& v)
		{
			seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}

		// Recursive template code derived from Matthieu M.
		template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
		struct HashValueImpl
		{
			static void apply(size_t& seed, Tuple const& tuple)
			{
				HashValueImpl<Tuple, Index - 1>::apply(seed, tuple);
				hash_combine(seed, get<Index>(tuple));
			}
		};

		template <class Tuple>
		struct HashValueImpl<Tuple, 0>
		{
			static void apply(size_t& seed, Tuple const& tuple)
			{
				hash_combine(seed, get<0>(tuple));
			}
		};
	}

	template <typename ... TT>
	struct hash<std::tuple<TT...>>
	{
		size_t
			operator()(std::tuple<TT...> const& tt) const
		{
			size_t seed = 0;
			HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
			return seed;
		}

	};
}

namespace majesty {

	//static Minisat::vec<Lit> spec_clause;

	static inline svar_map create_variables(sat_solver* solver, synth_spec* spec) {
		auto var_ctr = 0u;

		spec->selection_var_offset = var_ctr;
		unordered_map<tuple<unsigned, unsigned, unsigned>, int> svar_map;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					svar_map[make_tuple(i, j, k)] = spec->selection_var_offset + var_ctr++;
				}
			}
		}

		spec->gate_var_offset = var_ctr;
		const auto nr_gate_vars = spec->nr_gates * 3;
		var_ctr += nr_gate_vars;

		spec->simulation_var_offset = var_ctr;
		const auto nr_simulation_vars = spec->nr_gates * spec->tt_size;
		var_ctr += nr_simulation_vars;
		
		sat_solver_setnvars(solver, var_ctr);

		return svar_map;
	}

	logic_ntk size_optimum_ntk(uint64_t func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.


		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		optional<logic_ntk> ntk;
		spec->nr_gates = 1u;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		auto solver = sat_solver_new();
		while (true) {
			sat_solver_restart(solver);
			auto svar_map = create_variables(solver, spec);
			auto network_exists = exists_fanin_2_ntk(func, solver, spec, svar_map);
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
		
		return std::move(*ntk);
	}
	

	/*
	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_spec* opt, const unsigned nr_vars, const unsigned gate_size) {
		// TODO: Check if function is constant or variable. If so, return early.


		// Make the function normal if it isn't already
		bool invert = func & 1;
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		optional<logic_ntk> ntk;
		auto nr_gates = 1u;
		while (true) {
			Solver solver;
			auto network_exists = exists_fanin_2_ntk_ns(func, solver, opt, nr_vars, nr_gates);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk_ns(func, solver, nr_vars, nr_gates, invert);
				if (opt->verbose) {
					print_fanin_2_solution_ns(func, solver, nr_vars, nr_gates);
				}
				break;
			}
			++nr_gates;
		}
		
		return std::move(*ntk);
	}
	*/

	static inline void add_selection_clause(sat_solver* solver, synth_spec* spec, 
		const unordered_map<tuple<unsigned,unsigned,unsigned>,int>& svar_map,
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

		auto sel_var = svar_map.at(make_tuple(i, j, k));
		auto sel_lit = Abc_Var2Lit(sel_var, true);
		auto gate_lit = Abc_Var2Lit(simulation_variable(i, t, spec), a);
		plits[ctr++] = sel_lit;
		plits[ctr++] = gate_lit;
		
		if (b | c) {
			plits[ctr++] = Abc_Var2Lit(gate_variable(i, ((c << 1) | b) - 1, spec), !a);
		}
		sat_solver_addclause(solver, plits, plits + ctr);
	}

	static inline lbool cegar_solve(const uint64_t func, sat_solver* solver, synth_spec* spec, const svar_map& svar_map, 
		Vec_Int_t* vlits) {
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
								auto sel_var = svar_map.at(make_tuple(i, j, k));
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 0, 1);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 1, 0);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 1, 1);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 0, 0);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 0, 1);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 1, 0);
								add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 1, 1);
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

	lbool exists_fanin_2_ntk(const uint64_t func, sat_solver* solver, synth_spec* spec, const svar_map& svar_map) {
		static lit plits[3];

		// The gate's function constraint variables
		if (spec->no_triv_ops) {
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
		Vec_Int_t * vlits = Vec_IntAlloc(svar_map.size());
		for (auto i = 0u; i < spec->nr_gates; i++) {
			Vec_IntClear(vlits);
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					auto sel_var = svar_map.at(make_tuple(i, j, k));
					Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					if (!spec->use_cegar) {
						for (auto t = 0u; t < spec->tt_size; t++) {
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 0, 1);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 1, 0);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 0, 1, 1);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 0, 0);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 0, 1);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 1, 0);
							add_selection_clause(solver, spec, svar_map, t, i, j, k, 1, 1, 1);
						}
					}
				}
			}
			sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
			// Allow at most one selection variable to be true at a time
			if (spec->exact_nr_svars || spec->use_cegar) {
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
		if (spec->colex_order) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto j = 0u; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto jp = 0u; jp < spec->nr_vars + i; jp++) {
							for (auto kp = jp + 1; kp < spec->nr_vars + i; kp++) {
								if (k == kp && j > jp || k > kp) {
									auto sivar = svar_map.at(make_tuple(i, j, k));
									auto sipvar = svar_map.at(make_tuple(i+1, jp, kp));
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
						auto sel_var = svar_map.at(make_tuple(ip, j, i + spec->nr_vars));
						Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					}
					for (auto j = i + spec->nr_vars + 1; j < ip + spec->nr_vars; j++) {
						auto sel_var = svar_map.at(make_tuple(ip, i + spec->nr_vars, j));
						Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
					}
				}
				sat_solver_addclause(solver, Vec_IntArray(vlits), Vec_IntLimit(vlits));
			}
		}

		// Do not allow reapplication of operators
		if (spec->no_reapplication) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto ip = i + 1; ip < spec->nr_gates; ip++) {
					for (auto j = 0u; j < spec->nr_vars + i; j++) {
						for (auto k = j + 1; k < spec->nr_vars + i; k++) {
							auto sivar = svar_map.at(make_tuple(i, j, k));
							plits[0] = Abc_Var2Lit(sivar, true);
							auto sipvar = svar_map.at(make_tuple(i + 1, j, i + spec->nr_vars));
							plits[1] = Abc_Var2Lit(sipvar, true);
							sat_solver_addclause(solver, plits, plits + 2);
							sipvar = svar_map.at(make_tuple(i + 1, k, i + spec->nr_vars));
							plits[1] = Abc_Var2Lit(sipvar, true);
							sat_solver_addclause(solver, plits, plits + 2);
						}
					}
				}
			}
		}


		if (spec->use_cegar) {
			auto res = cegar_solve(func, solver, spec, svar_map, vlits);
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
	/*
	lbool exists_fanin_2_ntk_ns(const uint64_t func, Solver& solver, synth_spec* opt, const unsigned nr_vars, const unsigned nr_gates) {
		const auto tt_size = (1 << nr_vars) - 1;
		
		// Create variables that represent  the gates' truth tables
		auto nr_gate_vars = nr_gates * tt_size;
		for (auto i = 0u; i < nr_gate_vars; i++) {
			// Create gate variable x_it
			//cout << "adding gate tt var " << i << endl;
			solver.newVar();
		}

		// The gate's function constraint variables
		for (auto i = 0u; i < nr_gates; i++) {
			auto fi01 = solver.newVar();
			auto fi10 = solver.newVar();
			auto fi11 = solver.newVar();
			if (opt->no_triv_ops) {
				solver.addClause(mkLit(fi01, false), mkLit(fi10, false), mkLit(fi11, false));
				solver.addClause(mkLit(fi01, false), mkLit(fi10, true), mkLit(fi11, true));
				solver.addClause(mkLit(fi01, true), mkLit(fi10, false), mkLit(fi11, true));
			}
		}

		// Add selection (fanin) constraints
		static vector<Var> selection_vars;
		static Minisat::vec<Lit> selection_vars_clause;
		static unordered_map<tuple<unsigned, unsigned>, Var> svar_map;
		svar_map.clear();
		for (auto i = 0u; i < nr_gates; i++) {
			selection_vars.clear();
			const auto n = nr_vars + i;
			const auto k = 2;
			for (auto j = 0u; j < n; j++) {
				auto sel_var = solver.newVar();
				svar_map[make_tuple(i, j)] = sel_var;
				selection_vars.push_back(sel_var);
			}
			// Exactly 2 selection vars should be true, since every gate has fanin 2
			/* We do this by first preventing more than k vars from being true by selecting n choose k+1 subsets
			{
				vector<bool> v(n);
				std::fill(v.end() - (k + 1), v.end(), true);
				do {
					selection_vars_clause.clear();
					for (int j = 0; j < n; j++) {
						if (v[j]) {
							selection_vars_clause.push(mkLit(selection_vars[j], true));
						}
					}
					solver.addClause(selection_vars_clause);
				} while (std::next_permutation(v.begin(), v.end()));
			}
			// Next we prevent more than n-k vars from being false by selecting n choose n-k+1 subsets
			{
				vector<bool> v(n);
				std::fill(v.end() - (n - k + 1), v.end(), true);
				do {
					selection_vars_clause.clear();
					for (int j = 0; j < n; j++) {
						if (v[j]) {
							selection_vars_clause.push(mkLit(selection_vars[j], false));
						}
					}
					solver.addClause(selection_vars_clause);
				} while (std::next_permutation(v.begin(), v.end()));
			}

			// Finally, we add the constaint clauses
			for (auto j = 0; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
					for (auto t = 0u; t < tt_size; t++) {
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 0, 0, 1);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 0, 1, 0);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 0, 1, 1);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 1, 0, 0);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 1, 0, 1);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 1, 1, 0);
						add_selection_clause(solver, nr_vars, tt_size, selection_vars, nr_gate_vars, t, i, j, k, 1, 1, 1);
					}
				}
			}
		}

		// Enforce that a solution uses all gates
		if (opt->use_all_gates) {
			for (auto i = 0u; i < nr_gates - 1; i++) {
				selection_vars_clause.clear();
				for (auto ip = i + 1; ip < nr_gates; ip++) {
					auto sel_var = svar_map.at(make_tuple(ip, i + nr_vars));
					selection_vars_clause.push(mkLit(sel_var, false));
				}
				solver.addClause(selection_vars_clause);
			}
		}
		
		// The final gate's truth table should match the one from the specification
		static Minisat::vec<Lit> spec_clause;
		spec_clause.clear();
		for (auto t = 0u; t < tt_size; t++) {
			auto gate_var = gate_variable(tt_size, nr_gates - 1, t);
			spec_clause.push(mkLit(gate_var, !((func >> (t + 1)) & 1)));
		}

		if (solver.solve(spec_clause)) {
			return l_True;
		} else {
			return l_False;
		}
	}
	*/
	
	// Tries to finds a network with the specified size. Result is optional as this may
	// not be possible.
	/*
	lbool exists_fanin_2_ntk(const tt& func, Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		const auto tt_size = (1 << nr_vars) - 1;

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
			for (auto j = 0u; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
					auto sel_var = solver.newVar();
					sel_vars_clause.push(mkLit(sel_var, false));
					for (auto t = 0u; t < tt_size; t++) {
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 0, 1);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 0);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 0, 1, 1);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 0);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 0, 1);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 0);
						add_selection_clause(solver, nr_vars, tt_size, sel_var, nr_gate_vars, t, i, j, k, 1, 1, 1);
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
	*/
	
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
		return extract_fanin_2_ntk(func, solver, spec, false);
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
