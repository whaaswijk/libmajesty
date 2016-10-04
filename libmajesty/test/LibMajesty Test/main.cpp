#include <xmg.h>
#include <truth_table_utils.hpp>
#include <exact.h>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace majesty;
using Minisat::Solver;
using Minisat::lbool;
using cirkit::tt;

TEST_CASE("Trivial Exact Synthesis Tests", "[exact synthesis]") {
	Solver solver;
	tt andfunc(std::string("0001"));
	REQUIRE(exists_fanin_2_ntk(andfunc, solver, 1) == l_True);
}