#ifndef SHELL_ENV
#define SHELL_ENV

#include <string>
#include <vector>
#include <unordered_map>
#include <xmg.h>
#include <boost/optional.hpp>
#include <truth_table_utils.hpp>

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

namespace majesty
{
	class shell_env;
	
	enum command_code
	{
		success = 0,
		arg_error = 1, // error in command arguments
		cmd_error = 2, // error in executing command
		quit = 3
	};

	using command =
		std::function<command_code(shell_env*,
								   const std::vector<std::string>&)>;
	
	enum log_level
	{
		silent = 0,
		error = 1,
		warning = 2,
		message = 3
	};
	
	class shell_env 
	{
	private:
		std::unordered_map <std::string,command> command_map;

	public:
		std::vector<std::string> command_names;
		bool cfolder_exists;
		std::string history_file;
		log_level ll = log_level::warning;
		std::unique_ptr<xmg> current_ntk = nullptr;
		boost::optional<cirkit::tt> current_tt;
		
		int register_command(const std::string& name, command);
		int execute_command(const std::string& name,
							const std::vector<std::string>& argv);
		
		int print(char const* fmt, ...);
		int warning(char const* fmt, ...);
		int error(char const* fmt, ...);
		int message(char const* fmt, ...);
		
	};

	void destroy_shell(shell_env& env);
}

#endif
