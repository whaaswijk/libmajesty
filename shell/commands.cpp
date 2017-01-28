#include "shell_env.h"

using namespace std;

namespace majesty
{
	int test_command(shell_env* env, const vector<string>& argv)
	{
		env->warning("this is a test\n");
		return 0;
	}

	int test_command2(shell_env* env, const vector<string>& argv)
	{
		string str(argv[0]);
		for (auto i = 1u; i < argv.size(); i++)
		{
			str = str + " " + argv[i];
		}
		
		env->message("this is a test two \"%s\"\n", str.c_str());
		return 0;
	}
	
	int test_command3(shell_env* env, const vector<string>& argv)
	{
		env->error("this is a test\n");
		return 0;
	}

	int help_command(shell_env* env, const vector<string>& argv)
	{
		printf("\n");
		
		printf("Welcome to Majesty\n\n");
		printf("Unfortunately this help is still very incomplete...\n\n");
		
		printf("Prepend '!' to commands to interact with the system shell\n");
		printf("E.g. type \"!ls\" to list the files in the current directory\n");

		printf("\n");
		
		return 0;
	}
	
	void register_commands(shell_env& env) 
	{
		env.register_command("help", help_command);
		
	}
}
