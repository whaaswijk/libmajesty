#include "shell_env.h"
#include <readline/history.h>
#include <stdarg.h>

namespace majesty
{
	int
	shell_env::register_command(
		const std::string& name,
		command cmd)
	{
		auto it = command_map.find(name);
		if (it != command_map.end()) {
			warning("command \"%s\" already registered\n",
				 name.c_str());
			return 1;
		}
		
		command_map[name] = cmd;
		command_names.push_back(name);
		return 0;
	}
	
	int
	shell_env::execute_command(
		const std::string& name,
		const std::vector<std::string>& argv)
	{
		auto it = command_map.find(name);
		if (it == command_map.end()) {
			error("command \"%s\" not found\n", name.c_str());
			return 1;
		}
		return it->second(this, argv);
	}

	int
	shell_env::print(char const* fmt, ...)
	{
		va_list argptr;
		va_start(argptr, fmt);
		int ret = vprintf(fmt, argptr);
		va_end(argptr);
		return ret;
	}
	

	int
	shell_env::warning(char const* fmt, ...)
	{
		if (ll < log_level::warning)
		{
			return 0;
		}

		printf("%sWarning: %s", KMAG, KNRM);
		va_list argptr;
		va_start(argptr, fmt);
		int ret = vprintf(fmt, argptr);
		va_end(argptr);
		return ret;
	}
	
	
	int
	shell_env::error(char const* fmt, ...)
	{
		if (ll < log_level::error)
		{
			return 0;
		}
		
		printf("%sError: %s", KRED, KNRM);
		va_list argptr;
		va_start(argptr, fmt);
		int ret = vprintf(fmt, argptr);
		va_end(argptr);
		return ret;
	}
	
	int
	shell_env::message(char const* fmt, ...)
	{
		if (ll < log_level::error)
		{
			return 0;
		}
		
		printf("%sMessage: %s", KCYN, KNRM);
		va_list argptr;
		va_start(argptr, fmt);
		int ret = vprintf(fmt, argptr);
		va_end(argptr);
		return ret;
	}

	void
	destroy_shell(shell_env& env)
	{
		if (env.cfolder_exists)
		{
			write_history(env.history_file.c_str());
		}
	}
	
}
