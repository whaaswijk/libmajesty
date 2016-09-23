#ifndef CUT_H
#define CUT_H

#include <xmg.h>
#include <logic_network.h>
#include "truth_table_utils.hpp"
#include <cstdint>
#include <vector>
#include <map>
#include <unordered_map>
#include <limits>

#define DEFAULT_KLUT_SIZE 6
#define DEFAULT_CUT_LIMIT 8

const unsigned CHASH_BUF = sizeof(unsigned)*8;

namespace majesty {
	struct cut_params {
		unsigned int klut_size;
		unsigned int num_cuts;
	};

	std::unique_ptr<cut_params> default_cut_params();
	class cut;

	class cut {
		private:
			std::vector<nodeid> _nodes;
			unsigned _sig = 0;
			cut *_c1 = nullptr, *_c2 = nullptr, *_c3 = nullptr;
			bool _is_required = false;
			void compute_maj3_function(const node&,
				std::unordered_map<const cut*,std::unique_ptr<cirkit::tt>>&);
			void compute_xor_function(const node&,
				std::unordered_map<const cut*,std::unique_ptr<cirkit::tt>>&);

		public:
			cut();
			cut(const cut&);
			cut(nodeid);
			cut(cut*, cut*);
			cut(cut*, cut*, cut*);
			unsigned sig() const { return _sig; }
			unsigned int size() const { return _nodes.size(); }
			bool is_required() const { return _is_required; }
			void set_required();
			bool equal(const cut*) const;
			bool dominates(const cut*) const;
			const std::vector<nodeid>& nodes() const { return _nodes; }
			void computesignature();
			void computefunction(const node&,
					std::unordered_map<const cut*,std::unique_ptr<cirkit::tt>>&);

	};

	using cutvec = std::vector<std::unique_ptr<cut>>;
	using cutmap = std::vector<cutvec>;
	using reqmap = std::vector<unsigned int>;
	using funcmap = std::unordered_map<const cut*,std::unique_ptr<cirkit::tt>>;

	cutvec eqclass_cuts(const std::vector<node>&, const node&, const cutmap&, const cut_params*);
	cutvec eqclass_cuts_eval_funcs(const std::vector<node>&, const node&, const cutmap&, 
		const cut_params*, funcmap&, const std::vector<cirkit::tt>&);

	cutvec node_cuts(const node&, const cutmap&, const cut_params*);
	cutvec node_cuts(const ln_node&, const cutmap&, const cut_params*);
	
	cutvec crossproduct(const cutvec&, const cutvec&, const cut_params*);

	cutvec crossproduct(const cutvec&, const cutvec&, const cutvec&, 
			const cut_params*);

	void safeinsert(cutvec& v, std::unique_ptr<cut> c);

	cutmap enumerate_cuts(const xmg&, const cut_params*);
	cutmap enumerate_cuts(const logic_ntk&, const cut_params*);
	cutmap enumerate_cuts_eval_funcs(const xmg&, const cut_params*, funcmap&);

	cirkit::tt compute_function(const ln_node&, cutmap&);
}


#endif

