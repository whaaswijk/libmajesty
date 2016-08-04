#ifndef MIG_MAP_H
#define MIG_MAP_H

#include <vector>
#include "xmg.h"


namespace majesty {
	class strashmap {
		private:
			unsigned int _size = 0u;
			std::vector<std::vector<int32_t>> _table;
			xmg_stats& _stats;

		public:
			strashmap(size_t, xmg_stats&);
			int32_t find(maj3signature, xmg&) const;
			int32_t find_or_add(maj3signature, xmg&);
			int32_t find_or_add(maj3signature, xmg&, 
					varmap&, Minisat::Solver&, fanoutmap&);
			unsigned hash(maj3signature) const;

			int32_t find(xorsignature, xmg&) const;
			int32_t find_or_add(xorsignature, xmg&);
			int32_t find_or_add(xorsignature, xmg&, 
					varmap&, Minisat::Solver&, fanoutmap&);
			unsigned hash(xorsignature) const;

	};
}

#endif
