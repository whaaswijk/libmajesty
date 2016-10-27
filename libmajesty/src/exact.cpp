#include <exact.h>
#include <convert.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <tuple>

/*
extern "C" {
#include <sat/bsat/satSolver.h>
}
*/

using Minisat::Solver;
using Minisat::Var;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::lbool;

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

	static Minisat::vec<Lit> spec_clause;

	logic_ntk size_optimum_ntk(uint64_t func, synth_options* opt, const unsigned nr_vars, const unsigned gate_size) {
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
			auto network_exists = exists_fanin_2_ntk(func, solver, opt, nr_vars, nr_gates);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk(func, solver, nr_vars, nr_gates, invert);
				if (opt->verbose) {
					print_fanin_2_solution(func, solver, nr_vars, nr_gates);
				}
				break;
			}
			++nr_gates;
		}
		
		return std::move(*ntk);
	}

	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_options* opt, const unsigned nr_vars, const unsigned gate_size) {
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
			auto function_lit = mkLit(function_variable(3, i, ((c << 1) | b) - 1, nr_gate_vars), !a);
			clause_vec.push(function_lit);
		}

		solver.addClause(clause_vec);
	}

	static inline void add_selection_clause(Solver& solver, int num_vars, int tt_size, vector<Var>& selection_vars, int nr_gate_vars,
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

		auto sel_lit1 = mkLit(selection_vars[j], true);
		auto sel_lit2 = mkLit(selection_vars[k], true);
		auto gate_lit = mkLit(gate_variable(tt_size, i, t), a);
		clause_vec.push(sel_lit1);
		clause_vec.push(sel_lit2);
		clause_vec.push(gate_lit);
		
		if (b | c) {
			auto function_lit = mkLit(function_variable(3, i, ((c << 1) | b) - 1, nr_gate_vars), !a);
			clause_vec.push(function_lit);
		}

		solver.addClause(clause_vec);
	}

	static inline Var get_selection_var(unsigned var_offset, unsigned i, unsigned j, unsigned k) {
		return var_offset;
	}

	static inline lbool cegar_solve(const uint64_t func, Solver& solver, 
		unordered_map<tuple<unsigned, unsigned, unsigned>, Var>& svar_map, 
		unsigned tt_size, unsigned nr_gates, unsigned nr_vars, unsigned nr_gate_vars) {
		spec_clause.clear();

		while (true) {
			assert(spec_clause.size() <= tt_size);
			auto res = solver.solve(spec_clause);
			if (res == false) {
				return l_False;
			}

			// Extract the network, simulate it, and find the first bit where it's different from the specification
			auto ntk = extract_fanin_2_ntk(func, solver, nr_vars, nr_gates);
			auto ntk_func = ntk.simulate();

			bool found_solution = true;
			// Check if the solution found matches the specification
			for (auto t = 0u; t < tt_size; t++) {
				const auto ntk_val = ntk_func[t + 1];
				const bool spec_val = (func >> (t + 1)) & 1;
				if (ntk_val != spec_val) {
					// Constrain the solver further by adding an additional constraint for this truth table row
					const auto gate_var = gate_variable(tt_size, nr_gates - 1, t);
					spec_clause.push(mkLit(gate_var, !spec_val));
					for (auto i = 0u; i < nr_gates; i++) {
						for (auto j = 0u; j < nr_vars + i; j++) {
							for (auto k = j + 1; k < nr_vars + i; k++) {
								auto sel_var = svar_map.at(make_tuple(i, j, k));
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

	lbool exists_fanin_2_ntk(const uint64_t func, Solver& solver, synth_options* opt, const unsigned nr_vars, const unsigned nr_gates) {
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
		static Minisat::vec<Lit> sel_vars_clause;
		static unordered_map<tuple<unsigned, unsigned, unsigned>, Var> svar_map;
		svar_map.clear();
		for (auto i = 0u; i < nr_gates; i++) {
			sel_vars_clause.clear();
			for (auto j = 0u; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
					auto sel_var = solver.newVar();
					svar_map[make_tuple(i, j, k)] = sel_var;
					sel_vars_clause.push(mkLit(sel_var, false));
					if (!opt->use_cegar) {
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
			}
			solver.addClause(sel_vars_clause);
			// Allow at most one selection variable to be true at a time
			if (opt->exact_nr_svars || opt->use_cegar) {
				for (auto u = 0u; u < sel_vars_clause.size() - 1; u++) {
					auto svar1 = var(sel_vars_clause[u]);
					for (auto v = u + 1; v < sel_vars_clause.size(); v++) {
						auto svar2 = var(sel_vars_clause[v]);
						solver.addClause(mkLit(svar1, true), mkLit(svar2, true));
					}
				}
			}
		}

		// Check for co-lexicographic order: given gates i and i+1 with fanin (x, y) and (x', y') respectively,
		// we require y < y' OR y = y' AND x <= x'. This only if i is not a fanin of i+1.
		if (opt->colex_order) {
			for (auto i = 0u; i < nr_gates - 1; i++) {
				for (auto j = 0u; j < nr_vars + i; j++) {
					for (auto k = j + 1; k < nr_vars + i; k++) {
						for (auto jp = 0u; jp < nr_vars + i; jp++) {
							for (auto kp = jp + 1; kp < nr_vars + i; kp++) {
								if (k == kp && j > jp || k > kp) {
									auto sivar = svar_map.at(make_tuple(i, j, k));
									auto sipvar = svar_map.at(make_tuple(i+1, jp, kp));
									solver.addClause(mkLit(sivar, true), mkLit(sipvar, true));
								}
							}
						}
					}
				}
			}
		}

		// Enforce that a solution uses all gates
		if (opt->use_all_gates) {
			for (auto i = 0u; i < nr_gates - 1; i++) {
				sel_vars_clause.clear();
				for (auto ip = i + 1; ip < nr_gates; ip++) {
					for (auto j = 0u; j < i + nr_vars; j++) {
						auto sel_var = svar_map.at(make_tuple(ip, j, i + nr_vars));
						sel_vars_clause.push(mkLit(sel_var, false));
					}
					for (auto j = i + nr_vars + 1; j < ip + nr_vars; j++) {
						auto sel_var = svar_map.at(make_tuple(ip, i + nr_vars, j));
						sel_vars_clause.push(mkLit(sel_var, false));
					}
				}
				solver.addClause(sel_vars_clause);
			}
		}

		// Do not allow reapplication of operators
		if (opt->no_reapplication) {
			for (auto i = 0u; i < nr_gates - 1; i++) {
				for (auto ip = i + 1; ip < nr_gates; ip++) {
					for (auto j = 0u; j < nr_vars + i; j++) {
						for (auto k = j + 1; k < nr_vars + i; k++) {
							auto sivar = svar_map.at(make_tuple(i, j, k));
							auto sipvar = svar_map.at(make_tuple(i + 1, j, i + nr_vars));
							solver.addClause(mkLit(sivar, true), mkLit(sipvar, true));
							sipvar = svar_map.at(make_tuple(i + 1, k, i + nr_vars));
							solver.addClause(mkLit(sivar, true), mkLit(sipvar, true));
						}
					}
				}
			}
		}

		if (opt->use_cegar) {
			return cegar_solve(func, solver, svar_map, tt_size, nr_gates, nr_vars, nr_gate_vars);
		} else {
			// The final gate's truth table should match the one from the specification
			spec_clause.clear();
			for (auto t = 0u; t < tt_size; t++) {
				auto gate_var = gate_variable(tt_size, nr_gates - 1, t);
				//cout << "setting x_" << nr_gates - 1 << "_" << t << "(" << gate_var << "): " << func.test(t + 1) << endl;
				spec_clause.push(mkLit(gate_var, !((func >> (t + 1)) & 1)));
			}

			if (solver.solve(spec_clause)) {
				return l_True;
			} else {
				return l_False;
			}
		}
	}
	
	// Tries to find a network using the new selection variable implementation
	lbool exists_fanin_2_ntk_ns(const uint64_t func, Solver& solver, synth_options* opt, const unsigned nr_vars, const unsigned nr_gates) {
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
			}*/
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
	
	// Tries to finds a network with the specified size. Result is optional as this may
	// not be possible.
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
	
    logic_ntk extract_fanin_2_ntk(const tt& func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
        return extract_fanin_2_ntk(func, solver, nr_gates, nr_vars, false);
    }

	logic_ntk extract_fanin_2_ntk(const tt& func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates, bool invert) {
		logic_ntk ntk;

		const auto tt_size = func.size() - 1;
		
		for (auto i = 0u; i < nr_vars; i++) {
			ntk.create_input();
		}

		const auto nr_gate_vars = nr_gates * tt_size;
		const auto nr_function_vars = nr_gates * 3;
		auto var_offset = nr_gate_vars + nr_function_vars;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
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
	
	logic_ntk extract_fanin_2_ntk(const uint64_t func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		return extract_fanin_2_ntk(func, solver, nr_vars, nr_gates, false);
	}

	logic_ntk extract_fanin_2_ntk(const uint64_t func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates, bool invert) {
		logic_ntk ntk;

		const auto tt_size = (1u << nr_vars) - 1;
		
		for (auto i = 0u; i < nr_vars; i++) {
			ntk.create_input();
		}

		const auto nr_gate_vars = nr_gates * tt_size;
		const auto nr_function_vars = nr_gates * 3;
		auto var_offset = nr_gate_vars + nr_function_vars;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < nr_gates; i++) {
			for (auto j = 0u; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
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

	logic_ntk extract_fanin_2_ntk_ns(const tt& func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		return extract_fanin_2_ntk_ns(func, solver, nr_vars, nr_gates, false);
	}

	logic_ntk extract_fanin_2_ntk_ns(const tt& func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates, bool invert) {
		logic_ntk ntk;

		const auto tt_size = func.size() - 1;
		
		for (auto i = 0u; i < nr_vars; i++) {
			ntk.create_input();
		}

		const auto nr_gate_vars = nr_gates * tt_size;
		const auto nr_function_vars = nr_gates * 3;
		auto var_offset = nr_gate_vars + nr_function_vars;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			for (auto j = 0u; j < nr_vars + i; j++) {
				if (solver.modelValue(var_offset++) == l_True) {
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

	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		return extract_fanin_2_ntk(func, solver, nr_vars, nr_gates, false);
	}

	logic_ntk extract_fanin_2_ntk_ns(const uint64_t func, const Solver& solver, const unsigned nr_vars, const unsigned nr_gates, bool invert) {
		logic_ntk ntk;

		const auto tt_size = (1u << nr_vars) - 1;
		
		for (auto i = 0u; i < nr_vars; i++) {
			ntk.create_input();
		}

		const auto nr_gate_vars = nr_gates * tt_size;
		const auto nr_function_vars = nr_gates * 3;
		auto var_offset = nr_gate_vars + nr_function_vars;
		static vector<nodeid> fanin(2);
		static tt nodefunc(4);
		for (auto i = 0u; i < nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			for (auto j = 0u; j < nr_vars + i; j++) {
				if (solver.modelValue(var_offset++) == l_True) {
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

	void print_fanin_2_solution(const uint64_t func, Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		tt functt((1u << nr_vars), func);
		print_fanin_2_solution(functt, solver, nr_vars, nr_gates);
	}

	void print_fanin_2_solution(const tt& func, Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		const auto tt_size = (1u << nr_vars) - 1;

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
			for (auto j = 0u; j < nr_vars + i; j++) {
				for (auto k = j + 1; k < nr_vars + i; k++) {
					cout << "s_" << i << "_" << j << "_" << k << ": " << lbool_to_int(solver.modelValue(mkLit(var_offset++, false))) << endl;
				}
			}
		}
	}
	
	void print_fanin_2_solution_ns(const uint64_t func, Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		tt functt((1u << nr_vars), func);
		print_fanin_2_solution_ns(functt, solver, nr_vars, nr_gates);
	}

	void print_fanin_2_solution_ns(const tt& func, Solver& solver, const unsigned nr_vars, const unsigned nr_gates) {
		const auto tt_size = (1u << nr_vars) - 1;

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
			for (auto j = 0u; j < nr_vars + i; j++) {
				cout << "s_" << i << "_" << j << ": " << lbool_to_int(solver.modelValue(mkLit(var_offset++, false))) << endl;
			}
		}
	}

}
