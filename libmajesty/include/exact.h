#pragma once

#include <logic_network.h>
#include <truth_table_utils.hpp>
#include <boost/optional.hpp>
#include <tuple>
#include <sat_interface.h>
#include <vector>

using namespace cirkit;
	
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
		unsigned gate_size;
		unsigned gate_tt_size;
		unsigned selection_var_offset;
		unsigned nr_selection_vars;
		unsigned gate_var_offset;
		unsigned nr_gate_vars;
		unsigned simulation_var_offset;
		unsigned nr_simulation_vars;

		synth_spec() : verbose(false), use_cegar(false), use_colex_order(true),
			use_no_triv_ops(true), use_all_gates(true), use_exact_nr_svars(true),
			use_no_reapplication(true), use_colex_functions(true), gate_size(2u), gate_tt_size(3u) {}
	};

	static inline int selection_variable(const synth_spec* spec, unsigned gate_i, unsigned fanin_j, unsigned fanin_k) {
		auto ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					if (i == gate_i && j == fanin_j && k == fanin_k) {
						return spec->selection_var_offset + ctr;
					}
					++ctr;
				}
			}
		}
		return 0;
	}

	static inline int selection_variable_ns(const synth_spec* spec, unsigned gate_i, unsigned fanin_j) {
		return spec->selection_var_offset + fanin_j + gate_i * spec->nr_vars + (((gate_i - 1) * gate_i) / 2);
	}
	
	static inline int simulation_variable(int gate_i, int t, const synth_spec* spec) {
		return spec->tt_size * gate_i + t + spec->simulation_var_offset;
	}

	static inline int gate_variable(int gate, int idx, const synth_spec* spec) {
		return gate * spec->gate_tt_size + idx + spec->gate_var_offset;
	}
	
	// Tests a primary input's truth table (nonzero) at the specified index
	static inline bool pi_value(const int pi_num, const int tt_idx) {
		const auto num_zeros = (1u << pi_num);
		const auto proj_idx = tt_idx % (2 * num_zeros);
		return proj_idx >= num_zeros;
	}

	template<typename S>
	void print_fanin_2_solution(const uint64_t func, synth_spec* spec) {
		tt functt((1u << spec->nr_vars), func);
		print_fanin_2_solution<S>(functt, spec);
	}

	template<typename S>
	void print_fanin_2_solution(const tt& func, synth_spec* spec) {
		std::cout << "Solution for function " << to_string(func) << std::endl;

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto t = 0u; t < spec->tt_size; t++) {
				std::cout << "x_" << i << "_" << t << ": " << var_value<S>(simulation_variable(i, t, spec)) << std::endl;
			}
		}

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->gate_tt_size; j++) {
				std::cout << "f_" << i << "_" << j << ": " << var_value<S>(gate_variable(i, j, spec)) << std::endl;
			}
		}

		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					std::cout << "s_" << i << "_" << j << "_" << k << ": " << 
						var_value<S>(spec->selection_var_offset + selection_var_ctr++) << std::endl;
				}
			}
		}
	}
	
	template<typename S>
	void print_fanin_2_solution_ns(const uint64_t func, synth_spec* spec) {
		tt functt((1u << spec->nr_vars), func);
		print_fanin_2_solution_ns<S>(functt, spec);
	}

	template<typename S>
	void print_fanin_2_solution_ns(const tt& func, synth_spec* spec) {
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto t = 0u; t < spec->tt_size; t++) {
				std::cout << "x_" << i << "_" << t << ": " << var_value<S>(simulation_variable(i, t, spec)) << std::endl;
			}
		}

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->gate_tt_size; j++) {
				std::cout << "f_" << i << "_" << j << ": " << var_value<S>(gate_variable(i, j, spec)) << std::endl;
			}
		}

		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				std::cout << "s_" << i << "_" << j << ": " << var_value<S>(spec->selection_var_offset + selection_var_ctr++) << std::endl;
			}
		}
	}

	template<typename S>
	void print_fanin_3_solution_ns(const tt& func, synth_spec* spec) {
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto t = 0u; t < spec->tt_size; t++) {
				std::cout << "x_" << i << "_" << t << ": " << var_value<S>(simulation_variable(i, t, spec)) << std::endl;
			}
		}

		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->gate_tt_size; j++) {
				std::cout << "f_" << i << "_" << j << ": " << var_value<S>(gate_variable(i, j, spec)) << std::endl;
			}
		}

		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				std::cout << "s_" << i << "_" << j << ": " << var_value<S>(spec->selection_var_offset + selection_var_ctr++) << std::endl;
			}
		}
	}

	template<typename S>
	void print_fanin_3_solution_ns(const uint64_t func, synth_spec* spec) {
		tt functt((1u << spec->nr_vars), func);
		print_fanin_3_solution_ns<S>(functt, spec);
	}
	
	template<typename S>
	logic_ntk extract_fanin_2_ntk(synth_spec* spec) {
		return extract_fanin_2_ntk<S>(spec, false);
	}

	template<typename S>
	logic_ntk extract_fanin_2_ntk(synth_spec* spec, bool invert) {
		logic_ntk ntk;
		
		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static std::vector<nodeid> fanin(2);
		static tt nodefunc(4);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				for (auto k = j + 1; k < spec->nr_vars + i; k++) {
					if (var_value<S>(spec->selection_var_offset + selection_var_ctr++)) {
						// Gate i has gates j and k as fanin
						fanin[0] = j;
						fanin[1] = k;
						nodefunc.set(0, 0);
						nodefunc.set(1, var_value<S>(gate_variable(i, 0, spec)));
						nodefunc.set(2, var_value<S>(gate_variable(i, 1, spec)));
						nodefunc.set(3, var_value<S>(gate_variable(i, 2, spec)));
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

	template<typename S>
	logic_ntk extract_fanin_2_ntk_ns(synth_spec* spec) {
		return extract_fanin_2_ntk_ns<S>(spec, false);
	}

	template<typename S>
	logic_ntk extract_fanin_2_ntk_ns(synth_spec* spec, bool invert) {
		logic_ntk ntk;

		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static std::vector<nodeid> fanin(2);
		static tt nodefunc(4);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				if (var_value<S>(spec->selection_var_offset + selection_var_ctr++)) {
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
						nodefunc.set(1, var_value<S>(gate_variable(i, 0, spec)));
						nodefunc.set(2, var_value<S>(gate_variable(i, 1, spec)));
						nodefunc.set(3, var_value<S>(gate_variable(i, 2, spec)));
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

	template<typename S>
	logic_ntk extract_fanin_3_ntk_ns(synth_spec* spec) {
		return extract_fanin_3_ntk_ns<S>(spec, false);
	}
	
	template<typename S>
	logic_ntk extract_fanin_3_ntk_ns(synth_spec* spec, bool invert) {
		logic_ntk ntk;

		for (auto i = 0u; i < spec->nr_vars; i++) {
			ntk.create_input();
		}

		static std::vector<nodeid> fanin(3);
		static tt nodefunc(8);
		auto selection_var_ctr = 0u;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			auto inputs_found = 0u;
			auto in1 = 0u;
			auto in2 = 0u;
			auto in3 = 0u;
			for (auto j = 0u; j < spec->nr_vars + i; j++) {
				if (var_value<S>(spec->selection_var_offset + selection_var_ctr++)) {
					if (inputs_found == 0u) {
						in1 = j;
					} else if (inputs_found == 1u) {
						in2 = j;
					} else {
						in3 = j;
					}
					++inputs_found;
					if (inputs_found == 3) {
						fanin[0] = in1;
						fanin[1] = in2;
						fanin[2] = in3;
						nodefunc.set(0, 0);
						nodefunc.set(1, var_value<S>(gate_variable(i, 0, spec)));
						nodefunc.set(2, var_value<S>(gate_variable(i, 1, spec)));
						nodefunc.set(3, var_value<S>(gate_variable(i, 2, spec)));
						nodefunc.set(4, var_value<S>(gate_variable(i, 3, spec)));
						nodefunc.set(5, var_value<S>(gate_variable(i, 4, spec)));
						nodefunc.set(6, var_value<S>(gate_variable(i, 5, spec)));
						nodefunc.set(7, var_value<S>(gate_variable(i, 6, spec)));
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

	template<typename S>
	static inline void create_variables(synth_spec* spec) {
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
		const auto nr_gate_vars = spec->nr_gates * spec->gate_tt_size;
		spec->nr_gate_vars = nr_gate_vars;
		var_ctr += nr_gate_vars;

		spec->simulation_var_offset = var_ctr;
		const auto nr_simulation_vars = spec->nr_gates * spec->tt_size;
		spec->nr_simulation_vars = nr_simulation_vars;
		var_ctr += nr_simulation_vars;
		
		set_nr_vars<S>(var_ctr);
	}
	
	template<typename S>
	static inline void create_variables_ns(synth_spec* spec) {
		auto var_ctr = 0u;

		spec->selection_var_offset = var_ctr;
		const auto nr_selection_vars = spec->nr_gates * spec->nr_vars + (((spec->nr_gates - 1) * spec->nr_gates) / 2);
		var_ctr += nr_selection_vars;
		spec->nr_selection_vars = nr_selection_vars;

		spec->gate_var_offset = var_ctr;
		const auto nr_gate_vars = spec->nr_gates * spec->gate_tt_size;
		spec->nr_gate_vars = nr_gate_vars;
		var_ctr += nr_gate_vars;

		spec->simulation_var_offset = var_ctr;
		const auto nr_simulation_vars = spec->nr_gates * spec->tt_size;
		spec->nr_simulation_vars = nr_simulation_vars;
		var_ctr += nr_simulation_vars;
		
		set_nr_vars<S>(var_ctr);
	}

	template<typename S>
	static inline void add_selection_clause(synth_spec* spec,
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
		add_clause<S>(plits, plits + ctr);
	}

	template<typename S>
	static inline void add_selection_clause_ns(synth_spec* spec, 
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
		add_clause<S>(plits, plits + ctr);
	}

	template<typename S>
	static inline void add_selection_clause_ns(synth_spec* spec, 
		int t, int i, int j, int k, int l, bool a, bool b, bool c, bool d) {
		static lit plits[8];
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

		if (l < spec->nr_vars) {
			auto const_val = pi_value(l, t + 1);
			if (const_val != d) {
				return;
			}
		} else {
			plits[ctr++] = Abc_Var2Lit(simulation_variable(l - spec->nr_vars, t, spec), d);
		}

		auto sel_var1 = selection_variable_ns(spec, i, j);
		auto sel_var2 = selection_variable_ns(spec, i, k);
		auto sel_var3 = selection_variable_ns(spec, i, l);
		auto sel_lit1 = Abc_Var2Lit(sel_var1, true);
		auto sel_lit2 = Abc_Var2Lit(sel_var2, true);
		auto sel_lit3 = Abc_Var2Lit(sel_var3, true);
		auto gate_lit = Abc_Var2Lit(simulation_variable(i, t, spec), a);
		plits[ctr++] = sel_lit1;
		plits[ctr++] = sel_lit2;
		plits[ctr++] = sel_lit3;
		plits[ctr++] = gate_lit;
		
		if (b | c | d) {
			plits[ctr++] = Abc_Var2Lit(gate_variable(i, ((d << 2) | (c << 1) | b) - 1, spec), !a);
		}
		add_clause<S>(plits, plits + ctr);
	}
	
	template<typename S>
	static inline lbool cegar_solve(const uint64_t func, synth_spec* spec, Vec_Int_t* vlits) {
		tt functt(1 << spec->nr_vars, func);
		return cegar_solve<S>(functt, spec, vlits);
	}

	template<typename S>
	static inline lbool cegar_solve(const tt& func, synth_spec* spec, Vec_Int_t* vlits) {
		Vec_IntClear(vlits);

		while (true) {
			auto res = solve<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			if (res == l_False) {
				return res;
			}

			// Extract the network, simulate it, and find the first bit where it's different from the specification
			auto ntk = extract_fanin_2_ntk<S>(spec);
			auto ntk_func = ntk.simulate();

			bool found_solution = true;
			// Check if the solution found matches the specification
			for (auto t = 0u; t < spec->tt_size; t++) {
				const auto ntk_val = ntk_func[t + 1];
				const bool spec_val = func.test(t + 1);
				if (ntk_val != spec_val) {
					// Constrain the solver further by adding an additional constraint for this truth table row
					const auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
					Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !spec_val));
					for (auto i = 0u; i < spec->nr_gates; i++) {
						for (auto j = 0u; j < spec->nr_vars + i; j++) {
							for (auto k = j + 1; k < spec->nr_vars + i; k++) {
								add_selection_clause<S>(spec, t, i, j, k, 0, 0, 1);
								add_selection_clause<S>(spec, t, i, j, k, 0, 1, 0);
								add_selection_clause<S>(spec, t, i, j, k, 0, 1, 1);
								add_selection_clause<S>(spec, t, i, j, k, 1, 0, 0);
								add_selection_clause<S>(spec, t, i, j, k, 1, 0, 1);
								add_selection_clause<S>(spec, t, i, j, k, 1, 1, 0);
								add_selection_clause<S>(spec, t, i, j, k, 1, 1, 1);
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
	
	template<typename S>
	static inline lbool cegar_solve_ns(const uint64_t func, synth_spec* spec, Vec_Int_t* vlits) {
		tt functt(1 << spec->nr_vars, func);
		return cegar_solve_ns<S>(functt, spec, vlits);
	}

	template<typename S>
	static inline lbool cegar_solve_ns(const tt& func, synth_spec* spec, Vec_Int_t* vlits) {
		Vec_IntClear(vlits);

		while (true) {
			auto res = solve<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			if (res == l_False) {
				return res;
			}

			// Extract the network, simulate it, and find the first bit where it's different from the specification
			auto ntk = extract_fanin_2_ntk_ns<S>(spec);
			auto ntk_func = ntk.simulate();

			bool found_solution = true;
			// Check if the solution found matches the specification
			for (auto t = 0u; t < spec->tt_size; t++) {
				const auto ntk_val = ntk_func[t + 1];
				const bool spec_val = func.test(t + 1);
				if (ntk_val != spec_val) {
					// Constrain the solver further by adding an additional constraint for this truth table row
					const auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
					Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !spec_val));
					for (auto i = 0u; i < spec->nr_gates; i++) {
						for (auto j = 0u; j < spec->nr_vars + i; j++) {
							for (auto k = j + 1; k < spec->nr_vars + i; k++) {
								add_selection_clause_ns<S>(spec, t, i, j, k, 0, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, 0, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, 0, 1, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, 1, 0, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, 1, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, 1, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, 1, 1, 1);
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
	
	template<typename S>
	lbool exists_fanin_2_ntk(const uint64_t func, synth_spec* spec) {
		tt functt(1 << spec->nr_vars, func);
		return exists_fanin_2_ntk<S>(functt, spec);
	}

	template<typename S>
	lbool exists_fanin_2_ntk(const tt& func, synth_spec* spec) {
		static lit plits[3];
		
		create_variables<S>(spec);

		// The gate's function constraint variables
		if (spec->use_no_triv_ops) {
			for (auto i = 0u; i < spec->nr_gates; i++) {
				const auto func_var0 = gate_variable(i, 0, spec);
				const auto func_var1 = gate_variable(i, 1, spec);
				const auto func_var2 = gate_variable(i, 2, spec);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 0);
				add_clause<S>(plits, plits + 3);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 1);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
				
				plits[0] = Abc_Var2Lit(func_var0, 1);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
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
							add_selection_clause<S>(spec, t, i, j, k, 0, 0, 1);
							add_selection_clause<S>(spec, t, i, j, k, 0, 1, 0);
							add_selection_clause<S>(spec, t, i, j, k, 0, 1, 1);
							add_selection_clause<S>(spec, t, i, j, k, 1, 0, 0);
							add_selection_clause<S>(spec, t, i, j, k, 1, 0, 1);
							add_selection_clause<S>(spec, t, i, j, k, 1, 1, 0);
							add_selection_clause<S>(spec, t, i, j, k, 1, 1, 1);
						}
					}
				}
			}
			add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			// Allow at most one selection variable to be true at a time
			if (spec->use_exact_nr_svars || spec->use_cegar) {
				for (auto u = 0u; u < Vec_IntSize(vlits) - 1; u++) {
					auto svar1 = Abc_Lit2Var(Vec_IntEntry(vlits, u));
					for (auto v = u + 1; v < Vec_IntSize(vlits); v++) {
						auto svar2 = Abc_Lit2Var(Vec_IntEntry(vlits, v));
						plits[0] = Abc_Var2Lit(svar1, true);
						plits[1] = Abc_Var2Lit(svar2, true);
						add_clause<S>(plits, plits + 2);
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
									add_clause<S>(plits, plits + 2);
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
				add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
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
							add_clause<S>(plits, plits + 2);
							sipvar = selection_variable(spec, i + 1, k, i + spec->nr_vars);
							plits[1] = Abc_Var2Lit(sipvar, true);
							add_clause<S>(plits, plits + 2);
						}
					}
				}
			}
		}


		if (spec->use_cegar) {
			auto res = cegar_solve<S>(func, spec, vlits);
			Vec_IntFree(vlits);
			return res;
		} else {
			// The final gate's truth table should match the one from the specification
			Vec_IntClear(vlits);
			for (auto t = 0u; t < spec->tt_size; t++) {
				auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
				Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !func.test(t + 1)));
			}

			auto res = solve<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			Vec_IntFree(vlits);
			return res;
		}
	}
	
	// Tries to find a network using the new selection variable implementation
	template<typename S>
	lbool exists_fanin_2_ntk_ns(const uint64_t func, synth_spec* spec) {
		tt functt(1 << spec->nr_vars, func);
		return exists_fanin_2_ntk_ns<S>(functt, spec);
	}
	
	// Tries to find a network using the new selection variable implementation
	template<typename S>
	lbool exists_fanin_2_ntk_ns(const tt& func, synth_spec* spec) {
		static lit plits[4];
		
		create_variables_ns<S>(spec);

		// The gate's function constraint variables
		if (spec->use_no_triv_ops) {
			for (auto i = 0u; i < spec->nr_gates; i++) {
				const auto func_var0 = gate_variable(i, 0, spec);
				const auto func_var1 = gate_variable(i, 1, spec);
				const auto func_var2 = gate_variable(i, 2, spec);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 0);
				add_clause<S>(plits, plits + 3);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 1);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
				
				plits[0] = Abc_Var2Lit(func_var0, 1);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
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
				std::vector<bool> v(n);
				std::fill(v.end() - (k + 1), v.end(), true);
				do {
					Vec_IntClear(vlits);
					for (int j = 0u; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), true));
						}
					}
					add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}
			// Next we prevent more than n-k vars from being false by selecting n choose n-k+1 subsets
			{
				std::vector<bool> v(n);
				std::fill(v.end() - (n - k + 1), v.end(), true);
				do {
					Vec_IntClear(vlits);
					for (int j = 0u; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), false));
						}
					}
					add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}

			// Finally, we add the constaint clauses
			if (!spec->use_cegar) {
				for (auto j = 0; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto t = 0u; t < spec->tt_size; t++) {
							add_selection_clause_ns<S>(spec, t, i, j, k, 0, 0, 1);
							add_selection_clause_ns<S>(spec, t, i, j, k, 0, 1, 0);
							add_selection_clause_ns<S>(spec, t, i, j, k, 0, 1, 1);
							add_selection_clause_ns<S>(spec, t, i, j, k, 1, 0, 0);
							add_selection_clause_ns<S>(spec, t, i, j, k, 1, 0, 1);
							add_selection_clause_ns<S>(spec, t, i, j, k, 1, 1, 0);
							add_selection_clause_ns<S>(spec, t, i, j, k, 1, 1, 1);
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
				add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
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
								if ((k == kp && j > jp) || k > kp) {
									auto sivar1 = selection_variable_ns(spec, i, j);
									auto sivar2 = selection_variable_ns(spec, i, k);
									auto sipvar1 = selection_variable_ns(spec, i+1, jp);
									auto sipvar2 = selection_variable_ns(spec, i+1, kp);
									plits[0] = Abc_Var2Lit(sivar1, true);
									plits[1] = Abc_Var2Lit(sivar2, true);
									plits[2] = Abc_Var2Lit(sipvar1, true);
									plits[3] = Abc_Var2Lit(sipvar2, true);
									add_clause<S>(plits, plits + 4);
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
							add_clause<S>(plits, plits + 4);
							sipvar2 = selection_variable_ns(spec, ip, k);
							plits[3] = Abc_Var2Lit(sipvar2, true);
							add_clause<S>(plits, plits + 4);
						}
					}
				}
			}
		}
		
		if (spec->use_cegar) {
			auto res = cegar_solve_ns<S>(func, spec, vlits);
			Vec_IntFree(vlits);
			return res;
		} else {
			// The final gate's truth table should match the one from the specification
			Vec_IntClear(vlits);
			for (auto t = 0u; t < spec->tt_size; t++) {
				auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
				Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !func.test(t + 1)));
			}
			auto res = solve<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			Vec_IntFree(vlits);
			return res;
		}
	}
	
	template<typename S>
	lbool exists_fanin_3_ntk_ns(const uint64_t func, synth_spec* spec) {
		tt functt(1 << spec->nr_vars, func);
		return exists_fanin_3_ntk_ns<S>(functt, spec);
	}

	// Tries to find a network using the new selection variable implementation
	template<typename S>
	lbool exists_fanin_3_ntk_ns(const tt& func, synth_spec* spec) {
		static lit plits[4];
		
		create_variables_ns<S>(spec);

		// The gate's function constraint variables
		if (spec->use_no_triv_ops) {
			for (auto i = 0u; i < spec->nr_gates; i++) {
				/*
				const auto func_var0 = gate_variable(i, 0, spec);
				const auto func_var1 = gate_variable(i, 1, spec);
				const auto func_var2 = gate_variable(i, 2, spec);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 0);
				add_clause<S>(plits, plits + 3);

				plits[0] = Abc_Var2Lit(func_var0, 0);
				plits[1] = Abc_Var2Lit(func_var1, 1);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
				
				plits[0] = Abc_Var2Lit(func_var0, 1);
				plits[1] = Abc_Var2Lit(func_var1, 0);
				plits[2] = Abc_Var2Lit(func_var2, 1);
				add_clause<S>(plits, plits + 3);
				*/
			}
		}

		// Add selection (fanin) constraints
		Vec_Int_t * vlits = Vec_IntAlloc(spec->nr_selection_vars);
		const auto k = 3;
		for (auto i = 0u; i < spec->nr_gates; i++) {
			const auto n = spec->nr_vars + i;
			// Exactly k selection vars should be true, since every gate has fanin k
			// We do this by first preventing more than k vars from being true by selecting n choose k+1 subsets
			if (n > k && (spec->use_exact_nr_svars || spec->use_cegar)) {
				std::vector<bool> v(n);
				std::fill(v.end() - (k + 1), v.end(), true);
				do {
					Vec_IntClear(vlits);
					for (int j = 0u; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), true));
						}
					}
					add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}
			// Next we prevent more than n-k vars from being false by selecting n choose n-k+1 subsets
			{
				std::vector<bool> v(n);
				std::fill(v.end() - (n - k + 1), v.end(), true);
				do {
					Vec_IntClear(vlits);
					for (int j = 0u; j < n; j++) {
						if (v[j]) {
							Vec_IntPush(vlits,  Abc_Var2Lit(selection_variable_ns(spec, i, j), false));
						}
					}
					add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
				} while (std::next_permutation(v.begin(), v.end()));
			}

			// Finally, we add the constaint clauses
			if (!spec->use_cegar) {
				for (auto j = 0; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto l = k + 1; l < spec->nr_vars + i; l++) {
							for (auto t = 0u; t < spec->tt_size; t++) {
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 0, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 0, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 0, 1, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 1, 0, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 1, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 1, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 0, 1, 1, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 0, 0, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 0, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 0, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 0, 1, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 1, 0, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 1, 0, 1);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 1, 1, 0);
								add_selection_clause_ns<S>(spec, t, i, j, k, l, 1, 1, 1, 1);
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
					auto sel_var = selection_variable_ns(spec, ip, i + spec->nr_vars);
					Vec_IntPush(vlits, Abc_Var2Lit(sel_var, false));
				}
				add_clause<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			}
		}

		// Check for co-lexicographic order: given gates i and i+1 with fanin (x, y) and (x', y') respectively,
		// we require y < y' OR y = y' AND x <= x'. This only if i is not a fanin of i+1.
		if (false && spec->use_colex_order) {
			for (auto i = 0u; i < spec->nr_gates - 1; i++) {
				for (auto j = 0u; j < spec->nr_vars + i; j++) {
					for (auto k = j + 1; k < spec->nr_vars + i; k++) {
						for (auto jp = 0u; jp < spec->nr_vars + i; jp++) {
							for (auto kp = jp + 1; kp < spec->nr_vars + i; kp++) {
								if ((k == kp && j > jp) || k > kp) {
									auto sivar1 = selection_variable_ns(spec, i, j);
									auto sivar2 = selection_variable_ns(spec, i, k);
									auto sipvar1 = selection_variable_ns(spec, i+1, jp);
									auto sipvar2 = selection_variable_ns(spec, i+1, kp);
									plits[0] = Abc_Var2Lit(sivar1, true);
									plits[1] = Abc_Var2Lit(sivar2, true);
									plits[2] = Abc_Var2Lit(sipvar1, true);
									plits[3] = Abc_Var2Lit(sipvar2, true);
									add_clause<S>(plits, plits + 4);
								}
							}
						}
					}
				}
			}
		}

		// Do not allow reapplication of operators
		if (spec->use_no_reapplication) {
			/*
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
							add_clause<S>(plits, plits + 4);
							sipvar2 = selection_variable_ns(spec, ip, k);
							plits[3] = Abc_Var2Lit(sipvar2, true);
							add_clause<S>(plits, plits + 4);
						}
					}
				}
			}
			*/
		}
		
		if (false && spec->use_cegar) {
			auto res = cegar_solve_ns<S>(func, spec, vlits);
			Vec_IntFree(vlits);
			return res;
		} else {
			// The final gate's truth table should match the one from the specification
			Vec_IntClear(vlits);
			for (auto t = 0u; t < spec->tt_size; t++) {
				auto gate_var = simulation_variable(spec->nr_gates - 1, t, spec);
				Vec_IntPush(vlits, Abc_Var2Lit(gate_var, !(func.test(t + 1))));
			}
			auto res = solve<S>(Vec_IntArray(vlits), Vec_IntLimit(vlits));
			Vec_IntFree(vlits);
			return res;
		}
	}


	template<typename S>
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
		init_solver<S>();
		while (true) {
			restart_solver<S>();
			auto network_exists = exists_fanin_2_ntk<S>(func, spec);
			if (network_exists == l_True) {
				if (spec->verbose) {
					print_fanin_2_solution<S>(func, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		destroy_solver<S>();
		
		return spec->nr_gates;
	}

	template<typename S>
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
		init_solver<S>();
		while (true) {
			restart_solver<S>();
			auto network_exists = exists_fanin_2_ntk_ns<S>(func, spec);
			if (network_exists == l_True) {
				if (spec->verbose) {
					print_fanin_2_solution_ns<S>(func, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		destroy_solver<S>();
		
		return spec->nr_gates;
	}

	template<typename S>
	logic_ntk size_optimum_ntk(uint64_t func, synth_spec* spec) {
		tt functt(1 << spec->nr_vars, func);
		return size_optimum_ntk<S>(functt, spec);
	}
	
	template<typename S>
	logic_ntk size_optimum_ntk(const tt& spec_func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.
		tt func = spec_func;


		// Make the function normal if it isn't already
		auto invert = func.test(0);
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		logic_ntk ntk;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		spec->gate_tt_size = (1 << spec->gate_size) - 1;
		
		spec->nr_gates = 1u;
		init_solver<S>();
		while (true) {
			if (spec->verbose) {
				std::cout << "trying with " << spec->nr_gates << " gates\n";
			}
			restart_solver<S>();
			auto network_exists = exists_fanin_2_ntk<S>(func, spec);
			if (network_exists == l_True) {
				ntk = extract_fanin_2_ntk<S>(spec, invert);
				if (spec->verbose) {
					print_fanin_2_solution<S>(func, spec);
				}
				break;
			}
			++spec->nr_gates;
		}
		destroy_solver<S>();
		
		return std::move(ntk);
	}

	template<typename S>
	logic_ntk size_optimum_ntk_ns(uint64_t func, synth_spec* spec) {
		tt functt(1 << spec->nr_vars, func);
		return size_optimum_ntk_ns<S>(functt, spec);
	}
	
	template<typename S>
	logic_ntk size_optimum_ntk_ns(const tt& spec_func, synth_spec* spec) {
		// TODO: Check if function is constant or variable. If so, return early.
		tt func = spec_func;


		// Make the function normal if it isn't already
		auto invert = func.test(0);
		if (invert) {
			func = ~func;
		}

		// Else start by checking for networks with increasing nrs. of gates.
		logic_ntk ntk;
		spec->tt_size = (1 << spec->nr_vars) - 1;
		spec->gate_tt_size = (1 << spec->gate_size) - 1;

		spec->nr_gates = 1u;
		init_solver<S>();
		if (spec->verbose) {
			std::cout << "\nnr vars: " << spec->nr_vars << "\n";
		}
		while (true) {
			if (spec->verbose) {
				std::cout << "trying with " << spec->nr_gates << " gates\n";
			}
			restart_solver<S>();
			if (spec->gate_size == 3) {
				auto network_exists = exists_fanin_3_ntk_ns<S>(func, spec);
				if (network_exists == l_True) {
					ntk = extract_fanin_3_ntk_ns<S>(spec, invert);
					if (spec->verbose) {
						//print_fanin_3_solution_ns<S>(func, spec);
					}
					break;
				}
			} else {
				auto network_exists = exists_fanin_2_ntk_ns<S>(func, spec);
				if (network_exists == l_True) {
					ntk = extract_fanin_2_ntk_ns<S>(spec, invert);
					if (spec->verbose) {
						//print_fanin_2_solution_ns<S>(func, spec);
					}
					break;
				}
			}
			++spec->nr_gates;
		}
		destroy_solver<S>();
		
		return std::move(ntk);
	}
	

		
		
}
