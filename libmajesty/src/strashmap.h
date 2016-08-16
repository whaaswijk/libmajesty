#ifndef MIG_MAP_H
#define MIG_MAP_H

#include <vector>
#include "xmg.h"


namespace majesty {
	class strashmap {
		private:
			unsigned int _size = 0u;
			std::vector<std::vector<nodeid>> _table;
			xmg_stats& _stats;

		public:
			strashmap(size_t, xmg_stats&);
			nodeid find(maj3signature, xmg&) const;
			nodeid find_or_add(maj3signature, xmg&);
			nodeid find_or_add(maj3signature, xmg&, 
					varmap&, Minisat::Solver&, fanoutmap&);
			unsigned hash(maj3signature) const;

			nodeid find(xorsignature, xmg&) const;
			nodeid find_or_add(xorsignature, xmg&);
			nodeid find_or_add(xorsignature, xmg&, 
					varmap&, Minisat::Solver&, fanoutmap&);
			unsigned hash(xorsignature) const;

	};
}

#endif
