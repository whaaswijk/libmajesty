#include <xmg.h>
#include <truth_table_utils.hpp>
#include <exact.h>
#include <mig_interface.h>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <iostream>
#include <convert.h>

using namespace majesty;
using cirkit::tt;

TEST_CASE("Exact Synthesis Utility Functions", "[exact synthesis]") {
	REQUIRE(pi_value(0, 0) == false);
	REQUIRE(pi_value(0, 1) == true);
	REQUIRE(pi_value(0, 2) == false);
	REQUIRE(pi_value(5, 31) == false);
	REQUIRE(pi_value(5, 32) == true);
}

/*
TEST_CASE("Trivial Exact Synthesis", "[exact synthesis]") {
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
		auto ntk = size_optimum_ntk<sat_solver>(andfuncint, &spec);
		REQUIRE(ntk.ninternal() == 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk<sat_solver>(orfuncint, &spec);
		REQUIRE(ntk.ninternal() != 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk<sat_solver>(orfuncint, &spec);
		REQUIRE(ntk.ninternal() == 2);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk<sat_solver>(invorfuncint, &spec);
		REQUIRE(ntk.ninternal() != 1);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk<sat_solver>(invorfuncint, &spec);
		auto ninternal = ntk.ninternal();
		REQUIRE(ninternal == 2);
	}
	{
		spec.nr_vars = 3;
		auto ntk = size_optimum_ntk<sat_solver>(2, &spec);
		REQUIRE(ntk.ninternal() == 2);
	}
}
*/



TEST_CASE("New Selection Variable", "[exact synthesis]") {
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
		auto old_ntk = size_optimum_ntk<sat_solver>(f, &spec);
		auto old_size = old_ntk.ninternal();
		auto old_simvec = old_ntk.simulate();
		auto new_ntk = size_optimum_ntk_ns<sat_solver>(f, &spec);
		auto new_size = new_ntk.ninternal();
		auto new_simvec = new_ntk.simulate();
		if (old_size != new_size) {
			std::cout << "hur" << std::endl;
		}
		REQUIRE(old_simvec == new_simvec);
		REQUIRE(old_size == new_size);
	}
}

TEST_CASE("Logic Ntk String Serialization", "[serialization]") {
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
		auto old_ntk = size_optimum_ntk_ns<sat_solver>(f, &spec);
		auto old_size = old_ntk.ninternal();
		auto old_simvec = old_ntk.simulate();

		auto ntk_str = logic_ntk_to_string(old_ntk);

		auto str_ntk = string_to_logic_ntk(ntk_str);
		auto str_size = str_ntk.ninternal();
		auto str_simvec = str_ntk.simulate();
		if (old_size != str_size) {
			std::cout << "hur" << std::endl;
		}
		if (old_simvec != str_simvec) {
			std::cout << "dur" << std::endl;
		}
		REQUIRE(old_simvec == str_simvec);
		REQUIRE(old_size == str_size);
	}
}

#ifndef _WIN32
#include <maj_io.h>
#include <convert.h>
#include <lut_cover.h>
#include <lut_optimize.h>

TEST_CASE("CryptoMiniSaat Trivial Exact Synthesis", "[exact synthesis]") {
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
		auto ntk = size_optimum_ntk<CMSat::SATSolver>(andfuncint, &spec);
		REQUIRE(ntk.ninternal() == 1);
	}
}

TEST_CASE("CryptoMiniSat", "[exact synthesis]") {
	synth_spec spec;
	spec.nr_vars = 3;
	spec.verbose = false;
	spec.use_cegar = false;//true;
	spec.use_all_gates = true;
	spec.use_colex_order = true;
	spec.use_no_triv_ops = true;
	spec.use_exact_nr_svars = true;
	spec.use_no_reapplication = true;
	for (auto f = 0u; f < 256; f++) {
		auto old_ntk = size_optimum_ntk_ns<sat_solver>(f, &spec);
		auto old_size = old_ntk.ninternal();
		auto old_simvec = old_ntk.simulate();
		auto new_ntk = size_optimum_ntk_ns<CMSat::SATSolver>(f, &spec);
		auto new_size = new_ntk.ninternal();
		auto new_simvec = new_ntk.simulate();
		if (old_size != new_size) {
			std::cout << "hur" << std::endl;
		}
		REQUIRE(old_simvec == new_simvec);
		REQUIRE(old_size == new_size);
	}
}

TEST_CASE("XMG to Logic Network", "[conversion]") {
    auto xmg = read_verilog("../assets/adder.v");
    auto ntk = xmg_to_logic_ntk(xmg);
    write_blif(ntk, "adder_converted.blif");
}

/*
TEST_CASE("LUT Mapping", "[techmapping]") {
    auto xmg = read_verilog("../assets/adder.v");
    auto ntk = xmg_to_logic_ntk(xmg);

    auto cut_params = default_cut_params();
    cut_params->klut_size = 4;
    auto lut_ntk = lut_map_area(ntk, cut_params.get());
    write_blif(lut_ntk, "adder.blif");

    xmg = read_verilog("../assets/multiplier.v");
    ntk = xmg_to_logic_ntk(xmg);

    lut_ntk = lut_map_area(ntk, cut_params.get());
    write_blif(lut_ntk, "multiplier.blif");

    xmg = read_verilog("../assets/div.v");
    ntk = xmg_to_logic_ntk(xmg);

    lut_ntk = lut_map_area(ntk, cut_params.get());
    write_blif(lut_ntk, "div.blif");
}
*/

TEST_CASE("Logic Network LUT Decomposition", "[optimization]") {
    auto xmg = read_verilog("../assets/adder.v");
    auto ntk = xmg_to_logic_ntk(xmg);
    auto cut_params = default_cut_params();
    cut_params->klut_size = 5;
    auto lut_ntk = lut_map_area(ntk, cut_params.get());
    auto opt_ntk = logic_ntk_from_luts(lut_ntk);
    write_blif(opt_ntk, "decomp_adder.blif");
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
