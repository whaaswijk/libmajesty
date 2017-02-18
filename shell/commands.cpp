#include "shell_env.h"
#include <maj_io.h>
#include <xmg.h>
#include "game_commands.h"
#include <convert.h>

using namespace std;

namespace majesty
{
	command_code
	test_command(shell_env* env, const vector<string>& argv)
	{
		env->warning("this is a test\n");
		return success;
	}

	command_code
	test_command2(shell_env* env, const vector<string>& argv)
	{
		string str(argv[0]);
		for (auto i = 1u; i < argv.size(); i++)
		{
			str = str + " " + argv[i];
		}
		
		env->message("this is a test two \"%s\"\n", str.c_str());
		return success;
	}
	
	command_code
	test_command3(shell_env* env, const vector<string>& argv)
	{
		env->error("this is a test\n");
		return success;
	}

	command_code
	quit_command(shell_env* env, const vector<string>& argv)
	{
		return quit;
	}

	command_code
	help_command(shell_env* env, const vector<string>& argv)
	{
		printf("\n");
		
		printf("Welcome to Majesty\n\n");
		printf("Unfortunately this help is still very incomplete...\n\n");
		
		printf("Prepend '!' to commands to interact with the system shell\n");
		printf("E.g. type \"!ls\" to list the files in the current directory\n");

		printf("\n");
		
		return success;
	}

	command_code
	xmg_read_verilog(shell_env* env, const vector<string>& argv)
	{
		if (argv.size() != 2)
		{
			env->error("input file not specified\n");
			return arg_error;
		}

		auto input_filename = argv[1];
		auto in_xmg = ptr_read_verilog(input_filename);
		env->current_ntk.reset(in_xmg);
		
		return success;
	}

	command_code
	xmg_write_verilog(shell_env* env, const vector<string>& argv)
	{
		if (env->current_ntk == nullptr)
		{
			env->error("no current network available\n");
			return cmd_error;
		}
		if (argv.size() != 2)
		{
			env->error("output file not specified\n");
			return arg_error;
		}

		auto output_filename = argv[1];
		write_verilog(*env->current_ntk, output_filename);
		
		return success;
	}

	command_code
	print_stats(shell_env* env, const vector<string>& argv)
	{
		if (env->current_ntk == nullptr)
		{
			env->print("No network available\n");
			return cmd_error;
		}
		else
		{
			auto xmg_ptr = env->current_ntk.get();
			env->print("XMG\ti/o = %u/%u\tnd = %u\tlevels = %d\n",
					   xmg_ptr->nin(), xmg_ptr->nout(),
					   xmg_ptr->nnodes_proper(),
					   xmg_ptr->depth());
		}
		return success;
	}

	command_code
	read_truth(shell_env* env, const vector<string>& argv) 
	{
		if (argv.size() != 2)
		{
			env->error("no truth table provided\n");
			return cmd_error;
		}

		cirkit::tt table(argv[1]);
		env->current_tt = table;

		return success;
	}

	command_code
	tt_print_stats(shell_env* env, const vector<string>& argv)
	{
		if (!env->current_tt)
		{
			env->error("no truth table has been loaded\n");
			return cmd_error;	
		}
		auto num_vars = cirkit::tt_num_vars(*env->current_tt);
		
		env->print("Truth Table \tvars = %u\trepr = %s\n",
				   num_vars, cirkit::to_string(*env->current_tt).c_str());
		return success;
	}

	command_code
	tt_to_mig(shell_env* env, const vector<string>& argv)
	{
		if (!env->current_tt)
		{
			env->error("no truth table has been loaded\n");
			return cmd_error;	
		}
		auto num_vars = cirkit::tt_num_vars(*env->current_tt);
		auto decomp_mig = mig_shannon_decompose(num_vars, *env->current_tt);
		auto copy_mig = new xmg(decomp_mig);
		env->current_ntk.reset(copy_mig);		
		return success;
	}
	
 	void
	register_commands(shell_env& env) 
	{
		env.register_command("help", help_command);
		env.register_command("quit", quit_command);
		env.register_command("xmg_read_verilog", xmg_read_verilog);
		env.register_command("xmg_write_verilog", xmg_write_verilog);
		env.register_command("xmg_print_stats", print_stats);
		env.register_command("search_improvement", search_improvement);
		env.register_command("search_depth_improvement", search_depth_improvement);
		env.register_command("apply_move", apply_move);
		env.register_command("tt", read_truth);
		env.register_command("tt_print_stats", tt_print_stats);
		env.register_command("tt_to_mig", tt_to_mig);
	}
}
