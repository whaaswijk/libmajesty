#pragma once

#include <logic_network.h>
#include <truth_table_utils.hpp>
#include <boost/optional.hpp>
#include "minisat/Solver.h"
#include "minisat/SolverTypes.h"

namespace majesty {

	logic_ntk size_optimum_ntk(cirkit::tt& function, unsigned gate_size);
	Minisat::lbool exists_fanin_2_ntk(const cirkit::tt& func, Minisat::Solver&, const unsigned nr_gates);
	logic_ntk extract_fanin_2_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_gates);
	logic_ntk extract_fanin_3_ntk(const cirkit::tt& func, const Minisat::Solver&, const unsigned nr_gates);
	boost::optional<logic_ntk> find_fanin_2_ntk(const cirkit::tt& function, const unsigned nr_gates);
	boost::optional<logic_ntk> find_fanin_3_ntk(const cirkit::tt& function, const unsigned nr_gates);
	
	inline unsigned nr_ordered_tuples(unsigned tuple_size, unsigned bound) {
		if (bound <= tuple_size) {
			return 0;
		} else if (tuple_size == bound + 1) {
			return 1;
		}
		return nr_ordered_tuples(tuple_size, bound - 1) + nr_ordered_tuples(tuple_size - 1, bound - 1);
	}
}