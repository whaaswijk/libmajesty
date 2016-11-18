#ifndef FUNCTION_STORE_H
#define FUNCTION_STORE_H

#ifndef _WIN32
extern "C" {
	#include <hiredis.h>
}
#endif
#include "truth_table_utils.hpp"
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/format.hpp>
#include <exact.h>
#include <map>

namespace majesty {

	struct store_stats {
		std::map<unsigned, unsigned> vars_size_map;
	};

	class function_store_entry {
	public:
		std::string expr;
		bool is_exact;
		unsigned conflict_limit;

		function_store_entry() {
		}

		function_store_entry(std::string expr, bool is_exact, unsigned conflict_limit) : 
			expr(expr), is_exact(is_exact), conflict_limit(conflict_limit) {
		}

		std::string to_string() {
			return (boost::format("%s-%d-%u") % expr % is_exact % conflict_limit).str();
		}
	};

	function_store_entry to_function_store_entry(const std::string&);

	class function_store {
		private:
			boost::filesystem::path::string_type _logfile;
#ifndef _WIN32
			redisContext* _rcontext = NULL;
#endif
			bool _remove_log = false;

		public:
			function_store();
			function_store(const std::string&, const unsigned port);
			~function_store();
			std::string min_size_xmg(const cirkit::tt&);
			boost::optional<std::string> min_size_xmg(const cirkit::tt&, unsigned);
			boost::optional<std::string> heuristic_size_xmg(const cirkit::tt&, unsigned);
			void set_heuristic_size_xmg(const cirkit::tt&, const std::string&, unsigned);
			cirkit::tt npn_canon(const cirkit::tt&, 
					cirkit::tt& phase, std::vector<unsigned>& perm);
			boost::optional<unsigned> get_last_size(const cirkit::tt& func);

			void set_entry(const cirkit::tt& function, const function_store_entry& entry);
			void set_entry(const cirkit::tt& function, const std::string& expr, const bool is_exact, const unsigned conflict_limit);
			boost::optional<function_store_entry> get_entry(const cirkit::tt& function);
			
			//boost::optional<std::string> get_size_optimum_ntk_ns(const cirkit::tt&, synth_spec* spec, unsigned conflict_limit);

			//logic_ntk get_logic_ntk(const std::string& key);

			store_stats get_store_stats();
	};
}

#endif
