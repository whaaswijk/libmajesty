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

namespace majesty {
	
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
	};
}

#endif
