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

TEST_CASE("MIG Moves Test", "[game]") {
	const auto ninputs = 3;
	const auto nfuncs = (1u << (1u << ninputs));
	for (auto i = 0u; i < nfuncs; i++) {
		auto mig = mig_int_decompose(ninputs, i);
		auto moves = compute_moves(*mig);
		for (auto& move : moves) {
			auto new_mig = apply_move(*mig, move);
			assert(new_mig->equals(*mig));
		}
		std::cout << i << std::endl;

		delete mig;
	}
}
/*
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

	synth_options opts;
	opts.verbose = false;
	opts.use_cegar = true;
		
	uint64_t andfuncint = andfunc.to_ulong();
	uint64_t orfuncint = orfunc.to_ulong();
	uint64_t invorfuncint = invorfunc.to_ulong();
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(andfuncint, solver, &opts, 2, 1);
		REQUIRE(exists == l_True);
		auto ntk = size_optimum_ntk(andfuncint, &opts, 2, 2);
		REQUIRE(ntk.ninternal() == 1);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(orfuncint, solver, &opts, 3, 1);
		REQUIRE(exists == l_False);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(orfuncint, solver, &opts, 3, 2);
		REQUIRE(exists == l_True);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(invorfuncint, solver, &opts, 3, 1);
		REQUIRE(exists == l_False);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(invorfuncint, solver, &opts, 3, 2);
		REQUIRE(exists == l_True);
		//print_fanin_2_solution(invorfunc, solver, 2);
		auto ntk = size_optimum_ntk(invorfuncint, &opts, 3, 2);
		auto ninternal = ntk.ninternal();
		REQUIRE(ninternal == 2);
	}
	{
		Solver solver;
		auto exists = exists_fanin_2_ntk(2, solver, &opts, 3, 2);
		REQUIRE(exists == l_True);
		//print_fanin_2_solution(2, solver, 3, 2);
	}

}

TEST_CASE("New Selection Variable Test", "[exact synthesis]") {
	for (auto f = 0u; f < 256; f++) {
		auto nr_gates = 1u;
		while (true) {
			Solver solver_old, solver_new;
			synth_options opts;
			opts.colex_order = false;
			opts.no_triv_ops = false;
			opts.use_all_gates = false;
			opts.verbose = false;
			auto old_exists = exists_fanin_2_ntk(f, solver_old, &opts, 3, nr_gates);
			auto new_exists = exists_fanin_2_ntk_ns(f, solver_new, &opts, 3, nr_gates);
			REQUIRE(old_exists == new_exists);
			if (old_exists == l_True) {
				break;
			}
			++nr_gates;
		}
	}
}
*/

#ifndef _WIN32
#include <maj_io.h>
TEST_CASE("Exact Synthesis Extraction Test", "[exact synthesis]") {
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
		exists_fanin_2_ntk(andfunc, solver, 1);
		auto ntk = extract_fanin_2_ntk(andfunc, solver, 1);
		ntk.create_dummy_names();
		write_blif(ntk, "andfunc.blif");
	}
	{
		Solver solver;
		exists_fanin_2_ntk(orfunc, solver, 2);
		auto ntk = extract_fanin_2_ntk(orfunc, solver, 2);
		ntk.create_dummy_names();
		write_blif(ntk, "orfunc.blif");
	}
    {
		Solver solver;
		exists_fanin_2_ntk(invorfunc, solver, 2);
		auto ntk = extract_fanin_2_ntk(invorfunc, solver, 2, true);
		ntk.create_dummy_names();
		write_blif(ntk, "invorfunc.blif");
	}
}

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
