#pragma once

extern "C" {
#include <base/abc/abc.h>
#include <misc/vec/vecInt.h>
#include <misc/vec/vecPtr.h>
#include <sat/bsat/satSolver.h>
}

namespace majesty {
	struct synth_spec;

	template<typename S>
	void init_solver();

	template<typename S>
	void restart_solver();

	template<typename S>
	void set_nr_vars(unsigned nr_vars);

	template<typename S>
	void add_clause(lit* begin, lit* end);

	template<typename S>
	int solve(lit* begin, lit* end);

	template<typename S>
	int var_value(int var);

	template<typename S>
	void destroy_solver();
	
	static sat_solver* abc_solver = NULL;

	template<>
	inline void init_solver<sat_solver>() {
		assert(abc_solver == NULL);
		abc_solver = sat_solver_new();
	}

	template<>
	inline void restart_solver<sat_solver>() {
		sat_solver_restart(abc_solver);
	}

	template<>
	inline void set_nr_vars<sat_solver>(unsigned nr_vars) {
		sat_solver_setnvars(abc_solver, nr_vars);
	}

	template<>
	inline void add_clause<sat_solver>(lit* begin, lit* end) {
		sat_solver_addclause(abc_solver, begin, end);
	}

	template<>
	inline int var_value<sat_solver>(int var) {
		return sat_solver_var_value(abc_solver, var);
	}

	template<>
	inline int solve<sat_solver>(lit* begin, lit* end) {
		return sat_solver_solve(abc_solver, begin, end, 0, 0, 0, 0);
	}

	template<>
	inline void destroy_solver<sat_solver>() {
		assert(abc_solver != NULL);
		sat_solver_delete(abc_solver);
		abc_solver = NULL;
	}
}

