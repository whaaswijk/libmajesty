#include <xmg.h>
#include <truth_table_utils.hpp>
#include <exact.h>
#include <mig_interface.h>
//#include <maj_io.h>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <iostream>
#include <convert.h>

using namespace majesty;
/*
using Minisat::Solver;
using Minisat::lbool;
using Minisat::Var;
*/
using cirkit::tt;


//static inline Var gate_variable(int tt_size, int gate_i, int t);

//static inline Var function_variable(int gate_func_size, int f, int idx, int nr_gate_vars);

	// Tests a primary input's truth table (nonzero) at the specified index
//static inline bool pi_value(const int pi_num, const int tt_idx);

TEST_CASE("Exact Synthesis Utility Functions Test", "[exact synthesis]") {
	REQUIRE(pi_value(0, 0) == false);
	REQUIRE(pi_value(0, 1) == true);
	REQUIRE(pi_value(0, 2) == false);
	REQUIRE(pi_value(5, 31) == false);
	REQUIRE(pi_value(5, 32) == true);
}

TEST_CASE("Trivial Exact Synthesis Test", "[exact synthesis]") {
	auto pi1 = cirkit::tt_nth_var(0);
	auto pi2 = cirkit::tt_nth_var(1);
	auto pi3 = cirkit::tt_nth_var(2);

	auto andfunc = pi1 & pi2;
	auto orfunc = andfunc | pi3;
	auto invorfunc = andfunc | ~pi3; 
	cirkit::tt_to_minbase(andfunc);
	auto andfuncstr = cirkit::to_string(andfunc);
	cirkit::tt_to_minbase(orfunc);
	auto orfuncstr = cirkit::to_string(orfunc);
	cirkit::tt_to_minbase(invorfunc);
	invorfunc = ~invorfunc;
	auto invorfuncstr = cirkit::to_string(invorfunc);

	synth_spec spec;
	spec.verbose = true;
	spec.use_cegar = true;
		
	uint64_t andfuncint = andfunc.to_ulong();
	uint64_t orfuncint = orfunc.to_ulong();
	uint64_t invorfuncint = invorfunc.to_ulong();
	{
		spec.nr_vars = 2;
		auto ntk = size_optimum_ntk(andfuncint, &spec);
		REQUIRE(ntk.ninternal() == 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk(orfuncint, &spec);
		REQUIRE(ntk.ninternal() != 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk(orfuncint, &spec);
		REQUIRE(ntk.ninternal() == 2);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk(invorfuncint, &spec);
		REQUIRE(ntk.ninternal() != 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk(invorfuncint, &spec);
		auto ninternal = ntk.ninternal();
		REQUIRE(ninternal == 2);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk(2, &spec);
		REQUIRE(ntk.ninternal() == 2);
		//print_fanin_2_solution(2, solver, 3, 2);
	}

}

TEST_CASE("New Selection Variable Test", "[exact synthesis]") {
	synth_spec spec;
	spec.nr_vars = 3;
	spec.verbose = false;
	spec.use_cegar = true;
	spec.use_all_gates = true;
	spec.use_colex_order = true;
	spec.use_no_triv_ops = true;
	spec.use_exact_nr_svars = true;
	spec.use_no_reapplication = true;
	for (auto f = 0u; f < 256; f++) {
		auto old_ntk = size_optimum_ntk(f, &spec);
		auto old_size = old_ntk.ninternal();
		auto old_simvec = old_ntk.simulate();
		auto new_ntk = size_optimum_ntk(f, &spec);
		auto new_size = old_ntk.ninternal();
		auto new_simvec = old_ntk.simulate();
		if (old_size != new_size) {
			std::cout << "hur" << std::endl;
		}
		REQUIRE(old_simvec == new_simvec);
		REQUIRE(old_size == new_size);
	}
}

#ifndef _WIN32
#include <maj_io.h>

TEST_CASE("XMG String Serialization", "[serialization]") {
    // Get a MIG and use it to test on
    auto mig = mig_int_decompose(6, 20);
    mig->create_dummy_names();
    auto mig_verilog_str = mig->to_verilog();
    auto dez_mig = verilog_to_xmg(mig_verilog_str);
    REQUIRE(mig->equals(dez_mig));
    delete mig;
}
#endif
