#include <sat_interface.h>

namespace majesty {
	sat_solver* abc_solver = NULL;

	template<>
	void init_solver<sat_solver>() {
		assert(abc_solver == NULL);
		abc_solver = sat_solver_new();
	}

	template<>
	void restart_solver<sat_solver>() {
		sat_solver_restart(abc_solver);
	}

	template<>
	void set_nr_vars<sat_solver>(unsigned nr_vars) {
		sat_solver_setnvars(abc_solver, nr_vars);
	}

	template<>
	void add_clause<sat_solver>(lit* begin, lit* end) {
		sat_solver_addclause(abc_solver, begin, end);
	}

	template<>
	int var_value<sat_solver>(int var) {
		return sat_solver_var_value(abc_solver, var);
	}

	template<>
	int solve<sat_solver>(lit* begin, lit* end) {
		return sat_solver_solve(abc_solver, begin, end, 0, 0, 0, 0);
	}

	template<>
	void destroy_solver<sat_solver>() {
		assert(abc_solver != NULL);
		sat_solver_delete(abc_solver);
		abc_solver = NULL;
	}
}