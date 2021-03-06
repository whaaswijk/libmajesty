#include <cstdio>
#include <cstdlib>
#include <readline/readline.h>
#include <readline/history.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <boost/filesystem.hpp>
#include "majesty_shell.h"
#include "shell_env.h"
#include "commands.h"

using namespace std;
using namespace majesty;
namespace fs = boost::filesystem;

void print_header() 
{
	printf("%s%s%s\n", KWHT, header.c_str(), KNRM);
}

vector<string>* command_names_ref;

void init_shell(shell_env& env)
{
	char config_folder[256];
	
	// Check if the local configuration folder exists and try to
	// create it if not
	auto home = getenv("HOME");
	if (strlen(home) <= 0)
	{
		cout << "Warning: unable to find HOME" << endl;
		return;
	}
	snprintf(config_folder, 256, "%s/.majesty", home);
	fs::path config_path(config_folder);
	auto cfolder_exists = fs::is_directory(config_path);
	if (!cfolder_exists)
	{
		cfolder_exists =fs::create_directory(config_path);
	}
	if (!cfolder_exists)
	{
		cout << "Warning: unable to create config folder" << endl;
		env.cfolder_exists = false;
		return;
	}
	env.cfolder_exists = true;
	
	// Setup readline
	auto history_path = config_path.append("/history");
	env.history_file = history_path.string();
	if (fs::exists(history_path))
	{
		history_truncate_file(env.history_file.c_str(), 100);
		read_history(env.history_file.c_str());
	}
	else
	{
		cout << "No history file found" << endl;
	}

	register_commands(env);
	command_names_ref = &env.command_names;
}


vector<string>
parse_line(char* line)
{
	vector<string> argv;

	if (line[0] == '!') // Send command to system shell
	{
		system(line+1);
		return argv;
	}
	

	auto tokp = strtok(line, " ");
	while (tokp != NULL)
	{
		string arg(tokp);
		argv.push_back(arg);
		tokp = strtok(NULL, " ");
	}
	
	return argv;
}



char* command_name_generator(const char* text, int state)
{
	static size_t cmd_idx;
	static int len;
	
	if (!state)
	{
		cmd_idx = 0;
		len = strlen(text);

	}
	
	while (cmd_idx < command_names_ref->size()) {
		auto name = command_names_ref->at(cmd_idx).c_str();
		cmd_idx++;
		if (strncmp(name, text, len) == 0) {
			return strdup(name);
		}
	}

	return NULL;
}


char**
command_completion_function(const char* text, int start, int end)
{
	if (text[0] == '!')
	{
		return rl_completion_matches(text+1, rl_filename_completion_function);
	}
	else if (start == 0)
	{
		return rl_completion_matches(text, command_name_generator);
		return NULL;
	}
	else
	{
		return (char**)NULL;
	}
	
}

string
trim(const string& str)
{
	size_t first = str.find_first_not_of(' ');
	if (string::npos == first)
	{
		return str;
	}
	size_t last = str.find_last_not_of(' ');

	return str.substr(first, (last - first + 1));
}


int
main(void) 
{
	shell_env env;
	const char* prompt = "majesty > ";
	
	print_header();
	init_shell(env);

	rl_attempted_completion_function = command_completion_function;
	
	char *line = NULL;
	auto stop = false;
	while (!stop)
	{
		if ((line = readline(prompt)) == NULL)
		{
			break;
		}
		string line_str(line);
		auto trimmed_line_str = trim(line_str);
		auto trimmed_line = strdup(trimmed_line_str.c_str());
		add_history(trimmed_line);
		auto argv = parse_line(trimmed_line);
		free(trimmed_line);
		free(line);
		if (argv.size() > 0)
		{
			auto status = env.execute_command(argv[0], argv);
			switch (status)
			{
			case success:
				break;
			case quit:
				stop = true;
				break;
			default:
				break;
			}
		}
		
	}
	cout << endl;

	free(rl_prompt); // Hack to clean readline prompt seems to be necessary
	destroy_shell(env);
		
	return 0;
}
