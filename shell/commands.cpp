#include "shell_env.h"
#include <maj_io.h>
#include <xmg.h>
#include "game_commands.h"

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
	read_verilog_xmg(shell_env* env, const vector<string>& argv)
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
	write_verilog_xmg(shell_env* env, const vector<string>& argv)
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
		}
		else
		{
			auto xmg_ptr = env->current_ntk.get();
			env->print("XMG\ti/o = %u/%u\tnd = %u\n",
					  xmg_ptr->nin(), xmg_ptr->nout(),
					  xmg_ptr->nnodes_proper());
		}
		return success;
	}

  	void
	register_commands(shell_env& env) 
	{
		env.register_command("help", help_command);
		env.register_command("quit", quit_command);
		env.register_command("read_verilog", read_verilog_xmg);
		env.register_command("write_verilog", write_verilog_xmg);
		env.register_command("print_stats", print_stats);
		env.register_command("search_improvement", search_improvement);
		env.register_command("apply_move", apply_move);
	}
}
