#include "shell_env.h"
#include <xmg.h>
#include <mig_interface.h>

using namespace std;

namespace majesty
{
	bool dls(const xmg*m, unsigned orig_nnodes, unsigned depth)
	{
		if (depth == 0)
		{
			if (m->nnodes() < orig_nnodes)
				return true;
			else
				return false;
		}
		auto moves = compute_moves(*m);
		for (auto move : moves)
		{
			auto move_xmg = apply_move(*m, move);
			if (dls(move_xmg, orig_nnodes, depth-1))
				return true;
		}
		return false;
	}

	bool find_improvement(const xmg* m, unsigned orig_size, unsigned depth)
	{
		for (auto i = 1u; i <= depth; i++) {
			if (dls(m, orig_size, i))
				return true;
		}
		return false;
	}

	// Attempts a brute-force search on the MIG action space
	// to find a possible improvement within the specified
	// number of steps.
	command_code
	search_improvement(shell_env* env, const vector<string>& argv)
	{
		if (env->current_ntk == nullptr)
		{
			env->error("no MIG available\n");
			return cmd_error;
		}
		if (argv.size() != 2)
		{
			env->error("specify the depth of the search tree\n");
			return arg_error;
		}
		int depth = 0;
		try
		{
			depth = stoi(argv[1]);
			if (depth <= 0)
			{
				env->error("invalid depth\n");
				return arg_error;
			}
		}
		catch (invalid_argument e)
		{
			env->error("invalid depth\n");
			return arg_error;
		}

		auto xmg = env->current_ntk.get();
		auto found_improvement = find_improvement(xmg, xmg->nnodes(), depth);
		if (found_improvement)
			env->print("Found improvement\n");
		else
			env->print("No improvement found\n");

		return success;
	}
}
