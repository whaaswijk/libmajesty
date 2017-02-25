#ifndef GAME_COMMANDS
#define GAME_COMMANDS

#include "shell_env.h"

namespace majesty
{
	command_code
		search_improvement(shell_env* env, const std::vector<std::string>& argv);
	command_code
		search_depth_improvement(shell_env* env, const std::vector<std::string>& argv);
	command_code
		mig_apply_move(shell_env* env, const std::vector<std::string>& argv);
	command_code
		compute_nr_moves(shell_env* env, const std::vector<std::string>& argv);
	command_code
		compute_nr_moves_fast(shell_env* env, const std::vector<std::string>& argv);
}

#endif
