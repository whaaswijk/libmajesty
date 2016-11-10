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
}
