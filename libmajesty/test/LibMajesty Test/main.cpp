#include <xmg.h>
#include <truth_table_utils.hpp>
#include <exact.h>
//#include <maj_io.h>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace majesty;
using Minisat::Solver;
using Minisat::lbool;
using Minisat::Var;
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
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(andfunc, solver, 1);
		REQUIRE(exists == l_True);
		/*
		auto ntk = extract_fanin_2_ntk(andfunc, solver, 1);
		ntk.create_dummy_names();
		write_blif(ntk, "andfunc.blif");
		*/
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(orfunc, solver, 1);
		REQUIRE(exists == l_False);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(orfunc, solver, 2);
		/*
		auto ntk = extract_fanin_2_ntk(orfunc, solver, 2);
		ntk.create_dummy_names();
		write_blif(ntk, "orfunc.blif");
		*/
		REQUIRE(exists == l_True);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(invorfunc, solver, 1);
		REQUIRE(exists == l_False);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(invorfunc, solver, 2);
		//print_fanin_2_solution(invorfunc, solver, 2);
		REQUIRE(exists == l_True);
	}

}