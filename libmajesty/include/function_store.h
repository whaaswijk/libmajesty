#ifndef FUNCTION_STORE_H
#define FUNCTION_STORE_H

#ifndef _WIN32
extern "C" {
	#include <hiredis.h>
}
#endif
#include "truth_table_utils.hpp"
#include <boost/filesystem.hpp>

namespace majesty {
	
	enum NODE_TYPE { ALL, MAJ };

	class function_store {
		private:
			std::string compute_min_size_depth_xmg(const cirkit::tt&);
			std::string compute_min_size_depth_xmg(const cirkit::tt&, NODE_TYPE);
			boost::filesystem::path::string_type _logfile;
#ifndef _WIN32
			redisContext* _rcontext = NULL;
#endif
			bool _remove_log = false;

		public:
			function_store();
			function_store(const std::string&, const unsigned port);
			~function_store();
			std::string min_size_depth_xmg(const cirkit::tt&);
			std::string min_size_depth_xmg(const cirkit::tt&, NODE_TYPE);
			cirkit::tt npn_canon(const cirkit::tt&, 
					cirkit::tt& phase, std::vector<unsigned>& perm);
	};
}

#endif
