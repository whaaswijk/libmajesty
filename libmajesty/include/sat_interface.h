#pragma once

extern "C" {
#include <base/abc/abc.h>
#include <misc/vec/vecInt.h>
#include <misc/vec/vecPtr.h>
#include <sat/bsat/satSolver.h>
}

#ifdef USE_CMSAT
#include <cryptominisat5/cryptominisat.h>
// A hack to undefine the CryptoMiniSat lbool definitions, 
// as they conflict with those defined by ABC.
#undef l_True
#undef l_False
#undef l_Undef
#include <thread>
#endif

namespace majesty {

	struct synth_spec;

	static unsigned conflict_limit = 0u;

	template<typename S>
	void init_solver(unsigned conflict_limit);

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

	template<typename S>
	unsigned nr_conflicts();
	
	static sat_solver* abc_solver = NULL;

	template<>
	inline void init_solver<sat_solver>(unsigned climit) {
		assert(abc_solver == NULL);
		abc_solver = sat_solver_new();
		conflict_limit = climit;
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
		return sat_solver_solve(abc_solver, begin, end, conflict_limit, 0, 0, 0);
	}

	template<>
	inline void destroy_solver<sat_solver>() {
		assert(abc_solver != NULL);
		sat_solver_delete(abc_solver);
		abc_solver = NULL;
	}

	template<>
	inline unsigned nr_conflicts<sat_solver>() {
		return sat_solver_nconflicts(abc_solver);
	}

#ifdef USE_CMSAT
	static CMSat::SATSolver* cms_solver = nullptr;
	template<>
	inline void init_solver<CMSat::SATSolver>(unsigned climit) {
		assert(cms_solver == nullptr);
		conflict_limit = climit;
		/*
		cms_solver = new CMSat::SATSolver;
		auto nr_threads = std::thread::hardware_concurrency();
		cms_solver->set_num_threads(nr_threads);
		*/
	}

	template<>
	inline void restart_solver<CMSat::SATSolver>() {
		if (cms_solver != nullptr) {
			delete cms_solver;
		}
		cms_solver = new CMSat::SATSolver;
		auto nr_threads = std::thread::hardware_concurrency();
		cms_solver->set_num_threads(nr_threads);
		if (conflict_limit > 0) {
			cms_solver->set_max_confl(conflict_limit);
		}
	}

	template<>
	inline void set_nr_vars<CMSat::SATSolver>(unsigned nr_vars) {
		cms_solver->new_vars(nr_vars);
	}

	template<>
	inline void add_clause<CMSat::SATSolver>(lit* begin, lit* end) {
		static std::vector<CMSat::Lit> clause;
		clause.clear();
		for (auto i = begin; i < end; i++) {
			clause.push_back(CMSat::Lit(Abc_Lit2Var(*i), Abc_LitIsCompl(*i)));
		}
		cms_solver->add_clause(clause);
	}

	template<>
	inline int var_value<CMSat::SATSolver>(int var) {
		return cms_solver->get_model()[var] == CMSat::boolToLBool(true);
	}

	template<>
	inline int solve<CMSat::SATSolver>(lit* begin, lit* end) {
		static std::vector<CMSat::Lit> assumps;
		assumps.clear();
		for (auto i = begin; i < end; i++) {
			assumps.push_back(CMSat::Lit(Abc_Lit2Var(*i), Abc_LitIsCompl(*i)));
		}
		auto res = cms_solver->solve(&assumps);
		if (res == CMSat::boolToLBool(true)) {
			return l_True;
		} else if (res == CMSat::boolToLBool(false)) {
			return l_False;
		} else {
			return l_Undef;
		}
	}

	template<>
	inline void destroy_solver<CMSat::SATSolver>() {
		assert(cms_solver != nullptr);
		delete cms_solver;
		cms_solver = nullptr;
	}

	template<>
	inline unsigned nr_conflicts<CMSat::SATSolver>() {
		return 0u;
	}
#endif
}

