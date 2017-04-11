#include "shell_env.h"
#include <xmg.h>
#include <mig_interface.h>
#include <unordered_set>

using namespace std;

namespace majesty
{
	bool
	dfs(const xmg*m, unsigned orig_nnodes, unsigned depth)
	{
		if (depth == 0)
		{
			return m->nnodes() < orig_nnodes;
		}
		auto moves = compute_moves(*m);
		for (auto& move : moves)
		{
			auto move_xmg = apply_move(*m, move);
			auto strashed_xmg = strash_xmg(*move_xmg, true);
			auto found_improvement = dfs(strashed_xmg, orig_nnodes, depth-1);
			delete move_xmg;
			delete strashed_xmg;

			if (found_improvement) {
				printf("move(type=%d,id1=%u,id2=%u,id3=%u)\n",
					   move.type, move.nodeid1, move.nodeid2, move.nodeid3);
				return true;
			}
		}
		return false;
	}

	bool
	depth_dfs(const xmg*m, int orig_depth, unsigned search_depth)
	{
		if (search_depth == 0)
		{
			return m->depth() < orig_depth;
		}
		auto moves = compute_moves(*m);
		for (auto& move : moves)
		{
			auto move_xmg = apply_move(*m, move);
			auto strashed_xmg = strash_xmg(*move_xmg, true);
			auto found_improvement = depth_dfs(strashed_xmg, orig_depth, search_depth-1);
			delete move_xmg;
			delete strashed_xmg;
			if (found_improvement) {
				printf("move(type=%d,id1=%u,id2=%u,id3=%u)\n",
					   move.type, move.nodeid1, move.nodeid2, move.nodeid3);
				return true;
			}
		}
		return false;
	}

	bool
	find_improvement(const xmg* m, unsigned orig_size, unsigned depth)
	{
		return dfs(m, orig_size, depth);
	}

	bool
	find_depth_improvement(const xmg* m, unsigned orig_size, unsigned search_depth)
	{
		return depth_dfs(m, orig_size, search_depth);
	}
	

	// Attempts a brute-force search on the MIG action space to find
	// an improvement within the specified number of steps.
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

	// Attempts a brute-force search on the MIG action space to find
	// an improvement within the specified number of steps.
	command_code
	search_depth_improvement(shell_env* env, const vector<string>& argv)
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
		auto found_improvement = find_depth_improvement(xmg, xmg->depth(), depth);
		if (found_improvement)
			env->print("Found improvement\n");
		else
			env->print("No improvement found\n");

		return success;
	}

	command_code
	mig_apply_move(shell_env* env, const std::vector<std::string>& argv) 
	{
		if (env->current_ntk == nullptr)
		{
			env->error("no MIG available\n");
			return cmd_error;
		}		
		if (argv.size() < 2)
		{
			env->error("specify move type\n");
			return arg_error;
		}
		int move_type = 0;
		int nodeid1 = -1;
		int nodeid2 = -1;
		int nodeid3 = -1;
		
		try
		{
			move_type = stoi(argv[1]);
			if (move_type < 0)
			{
				env->error("invalid move type\n");
				return arg_error;
			}
		}
		catch (invalid_argument e)
		{
			env->error("invalid move type\n");
			return arg_error;
		}
		if (argv.size() > 2) {
			try
			{			
				nodeid1 = stoi(argv[2]);
				if (nodeid1 < 0)
				{
					env->error("invalid node id\n");
					return arg_error;
				}
			}
			catch (invalid_argument e)
			{
				env->error("invalid node id\n");
				return arg_error;
			}		
		}
		if (argv.size() > 3) {
			try
			{			
				nodeid2 = stoi(argv[3]);
				if (nodeid2 < 0)
				{
					env->error("invalid node id\n");
					return arg_error;
				}
			}
			catch (invalid_argument e)
			{
				env->error("invalid node id\n");
				return arg_error;
			}		
		}
		if (argv.size() > 4) {
			try
			{			
				nodeid3 = stoi(argv[4]);
				if (nodeid3 < 0)
				{
					env->error("invalid node id\n");
					return arg_error;
				}
			}
			catch (invalid_argument e)
			{
				env->error("invalid node id\n");
				return arg_error;
			}		
		}
		
		move m;
		m.type = (MoveType)move_type;
		m.nodeid1 = nodeid1;
		m.nodeid2 = nodeid2;
		m.nodeid3 = nodeid3;

		auto res_mig = apply_move(*env->current_ntk, m);
		if (res_mig == NULL)
		{
			env->error("invalid move\n");
			return arg_error;
		}
		
		auto strash_mig = strash_xmg(*res_mig, true);
		delete res_mig;

		for (auto inname : env->current_ntk->innames())
		{
			strash_mig->add_inname(inname);
		}
		for (auto outname : env->current_ntk->outnames())
		{
			strash_mig->add_inname(outname);
		}
		
		
		env->current_ntk.reset(strash_mig);
		
		return success;
	}

	command_code
	compute_nr_moves(shell_env* env, const vector<string>& argv)
	{
		if (env->current_ntk == nullptr)
		{
			env->error("no MIG available\n");
			return cmd_error;
		}

		auto moves = compute_moves(*env->current_ntk);
		env->print("%d moves available\n", moves.size());

		return success;
	}
	
	command_code
	compute_nr_moves_fast(shell_env* env, const vector<string>& argv)
	{
		if (env->current_ntk == nullptr)
		{
			env->error("no MIG available\n");
			return cmd_error;
		}

		auto moves = compute_moves_fast(*env->current_ntk);
		env->print("%d moves available\n", moves.size());

		return success;
	}
	
}
